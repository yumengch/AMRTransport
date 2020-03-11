module class_mesh
  use class_coord,     only: coord
  use class_AMRTypes,  only: FaceList, ColumnList
  use class_parameter, only: dimX, dimY
  implicit none
  private
  public :: mesh, mesh_init, mesh_sync, mesh_initOther, mesh_clean
  ! This class act as a glue among linked_list, quadtree, edges and centers
  type mesh
    ! coarse coordinate
    type(coord)                   :: coord_
    ! basic elements
    type(FaceList)                :: faces_(dimX:dimY)
    type(ColumnList)              :: columns_
  end type mesh

  interface mesh
    procedure :: mesh_init
  end interface mesh
  
  contains
    function mesh_init(nx, ny, nz, dt, beta, TestName) result(val)
    ! initialize a new mesh
      use class_AMRTypes,           only: face, FaceConnect
      use class_face,               only: face_init, FaceList_append, FaceConnect_append
      use class_column,             only: ColumnList_init
      use class_kind,               only: dp
      integer,             intent(in)  :: nx, ny, nz
      character(len =*),   intent(in)  :: TestName
      real(dp),            intent(in)  :: dt, beta  
      integer                          :: i, j
      type(mesh)                       :: val
      type(face)                       :: NewFace
      type(FaceConnect)                :: ColumnFace(0:nx -1, 0:ny-1, dimX:dimY, 2)
      type(face), pointer              :: null_face
      null_face => null()
      ! initialize the base coordinate  
      val%coord_ = coord(nx, ny, nz)
      ! initialize the faces
      do j = 0, ny - 1
        do i = 0, nx - 1
          NewFace = face_init(val%faces_(dimX)%tail_, val%coord_%x_(i), &
            [val%coord_%y_(j), val%coord_%y_(j + 1)], [i, j], val%coord_%l0_, &
            nz, TestName, 0, dt, beta, dimX, null_face)
          call FaceList_append(val%faces_(dimX), NewFace )
          call FaceConnect_append(ColumnFace(i, j, dimX, 1), val%faces_(dimX)%tail_, 2)
          call FaceConnect_append(ColumnFace(modulo(i - 1, nx), j, dimX, 2), val%faces_(dimX)%tail_, 1)
        enddo
      enddo
      ! initialize the south pole faces
      j = 0
      do i = 0, nx/2 - 1
        NewFace = face_init(val%faces_(dimY)%tail_, val%coord_%y_(j), &
          [val%coord_%x_(i), val%coord_%x_(i + 1)], [i, j], val%coord_%l0_, &
            nz, TestName, 0, dt, beta, dimY, null_face)
        call FaceList_append(val%faces_(dimY), NewFace )
        call FaceConnect_append(ColumnFace(i, j, dimY, 1), val%faces_(dimY)%tail_, 2)
        NewFace = face_init(val%faces_(dimY)%tail_, val%coord_%y_(j), [val%coord_%x_(i + nx/2), &
          val%coord_%x_(i + nx/2 + 1)], [i + nx/2, j], val%coord_%l0_, &
            nz, TestName, 0, dt, beta, dimY, null_face)
        NewFace%pole_ => val%faces_(dimY)%tail_
        call FaceList_append(val%faces_(dimY), NewFace )
        val%faces_(dimY)%tail_%prev_%pole_ => val%faces_(dimY)%tail_
        call FaceConnect_append(ColumnFace(i + nx/2, j, dimY, 1), val%faces_(dimY)%tail_, 2)
      enddo
      ! initialize non-pole faces in Y direction
      do j = 1, ny - 1
        do i = 0, nx - 1
          NewFace = face_init(val%faces_(dimY)%tail_, val%coord_%y_(j), [val%coord_%x_(i), &
            val%coord_%x_(i + 1)], [i, j], val%coord_%l0_, &
            nz, TestName, 0, dt, beta, dimY, null_face)
          call FaceList_append(val%faces_(dimY), NewFace )
          call FaceConnect_append(ColumnFace(i, j, dimY, 1), val%faces_(dimY)%tail_, 2)
          call FaceConnect_append(ColumnFace(i, j - 1, dimY, 2), val%faces_(dimY)%tail_, 1)
        enddo
      enddo
      ! initialize the north pole faces
      j = ny
      do i = 0, nx/2 - 1
        NewFace = face_init(val%faces_(dimY)%tail_, val%coord_%y_(j), &
          [val%coord_%x_(i), val%coord_%x_(i + 1)], [i, j], val%coord_%l0_, &
            nz, TestName, 0, dt, beta, dimY, null_face)
        call FaceList_append(val%faces_(dimY), NewFace )
        call FaceConnect_append(ColumnFace(i,j - 1, dimY, 2), val%faces_(dimY)%tail_, 1)
        NewFace = face_init(val%faces_(dimY)%tail_, val%coord_%y_(j), [val%coord_%x_(i + nx/2), &
          val%coord_%x_(i + nx/2 + 1)], [i + nx/2, j], val%coord_%l0_, &
            nz, TestName, 0, dt, beta, dimY, null_face)
        NewFace%pole_ => val%faces_(dimY)%tail_
        call FaceList_append(val%faces_(dimY), NewFace )
        val%faces_(dimY)%tail_%prev_%pole_ => val%faces_(dimY)%tail_
        call FaceConnect_append(ColumnFace(i + nx/2, j - 1, dimY, 2), val%faces_(dimY)%tail_, 1)
      enddo
      ! initialize the columns using the present faces
      call ColumnList_init(val%columns_, val%coord_, nx, ny, ColumnFace)
    end function mesh_init

    subroutine mesh_sync(self, nz, t, dt, beta)
    ! synchronize the old mesh with new mesh to prepare for next time step
      use class_AMRTypes,      only: column, ColumnStack
      use class_column,        only: ColumnList_Coarsen, ColumnList_refine
      use class_kind,          only: dp
      type(mesh),          intent(inout) :: self
      real(dp),            intent(in)    :: dt, beta
      integer,             intent(in)    :: t, nz
      type(column), pointer              :: MeshColumn, MeshColumnNext, StackColumn
      integer                            :: i, j
      type(ColumnStack)                  :: StackColumns
      MeshColumn => self%columns_%head_
      do while(associated(MeshColumn))

        MeshColumnNext => MeshColumn%next_
        ! if the new mesh is coarser than the current mesh coarsen it
        if (MeshColumn%OtherMesh_%l_ < MeshColumn%l_) then
          MeshColumn%prop_%NumChange_ = MeshColumn%OtherMesh_%l_ - MeshColumn%l_
          call ColumnList_Coarsen(self%columns_, self%faces_, MeshColumn, nz, "Old", t, dt, beta)
          MeshColumn%prop_%NumChange_ = 0 
          MeshColumnNext => MeshColumn%next_
        else if (associated(MeshColumn%OtherMesh_%child_)) then
        ! if the new mesh is finer than the current mesh, refine all the leaves
          StackColumns%head_ => MeshColumn
          do while (associated(StackColumns%head_))
            StackColumn => StackColumns%head_
            StackColumns%head_ => StackColumn%nextStack_
            StackColumn%nextStack_ => null()
            StackColumn%prop_%NumChange_ = 1
            call ColumnList_refine(self%columns_, StackColumn, self%faces_, nz, "Old", dt, beta, t)
            do i = 0, 1
              do j = 0, 1
                if (associated(StackColumn%child_(i, j)%OtherMesh_%child_)) then
                  StackColumn%child_(i, j)%nextStack_ => StackColumns%head_
                  StackColumns%head_ => StackColumn%child_(i, j)
                endif
              enddo
            enddo
          enddo
        endif
        MeshColumn => MeshColumnNext
      enddo
    end subroutine mesh_sync

    subroutine mesh_initOther(self, OtherMesh)
    ! initialize the other mesh property in the column
      use class_AMRTypes,      only: column
      type(mesh), intent(inout) :: self, OtherMesh
      type(column), pointer     :: MeshColumn, OtherMeshColumn
      integer                   :: i      
      MeshColumn => self%columns_%head_
      OtherMeshColumn => OtherMesh%columns_%head_
      do i = 1, self%columns_%n_
        MeshColumn%OtherMesh_ => OtherMeshColumn
        OtherMeshColumn%OtherMesh_ => MeshColumn
        MeshColumn => MeshColumn%next_
        OtherMeshColumn => OtherMeshColumn%next_
      enddo
    end subroutine mesh_initOther

    subroutine mesh_clean(self)
    ! finalize the mesh
      use class_coord, only: coord_clean
      use class_face,  only: FaceList_clean
      type(mesh), intent(inout) :: self
      call coord_clean(self%coord_)
      call FaceList_clean(self%faces_(dimX))
      call FaceList_clean(self%Faces_(dimY))
      ! column list clean
    end subroutine mesh_clean
end module class_mesh
