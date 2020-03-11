module class_column
  use class_AMRTypes, only: column, ColumnList, ColumnStack, cproperty, ColumnArray
  implicit none
  private
  public :: column_dA, column_get_leaves, column_stability
  public :: ColumnList_get_column, ColumnList_init, ColumnList_refine, ColumnList_Coarsen
  public :: ColumnList_AMR_Mark, ColumnList_q

  contains
    ! initialize a column instance
    function property_init(x, dx, i, Faces) result( val )
      use class_kind,                 only: dp
      use class_parameter,            only: dimX, dimY
      use class_AMRTypes,             only: FaceConnect, cproperty
      type(FaceConnect),  intent(in)     :: Faces(dimX:, 1:)
      real(dp),           intent(in)     :: x(dimX:), dx(dimX:)
      integer,            intent(in)     :: i(dimX:)
      type(cproperty)                    :: val

      ! assign edges to the cell
      val%x_ = x;
      val%dx_(:) = dx(:)
      val%i_(:) = i(:);
      val%cosp_ = (sin(x(dimY) + dx(dimY)) - sin(x(dimY)))/dx(dimY)
      val%faces_(1:, 1:) = Faces(1:, 1:)
      val%qmid_(:, :, :) = 0._dp
      val%div_(:) = 0._dp
    end function property_init

    ! get the area of the current column
    function column_dA(self) result(val)
      use class_kind,                 only: dp
      use class_parameter,            only: dimX, dimY
      type(column), intent(in) :: self
      real(dp)                  :: val
      val = self%prop_%dx_(dimX)*(sin(self%prop_%x_(dimY) + self%prop_%dx_(dimY)) - sin(self%prop_%x_(dimY)))
    end function column_dA

    ! get parent of the column
    function column_coarse(self, n) result(val)
      type(column),  pointer, intent(in) :: self
      integer,                intent(in) :: n
      integer                            :: i
      type(column),  pointer             :: val
      val => self
      do i = 1, n
        val => val%parent_
      enddo
    end function column_coarse

    ! coarsen the faces in the coarsend columns
    subroutine column_FaceCoarsen(CoarseColumn, Faces, NewFaces, x, q, nz, TestName, t, dt, beta)
      use class_AMRTypes,                    only: FaceList, FaceConnect, face
      use class_face,                        only: FaceConnect_appendList, FaceConnect_remove, FaceList_Remove
      use class_kind,                        only: dp
      use class_parameter,                   only: dimX, dimY
      type(column),      pointer, intent(inout) :: CoarseColumn
      type(FaceList),             intent(inout) :: Faces(dimX:)
      type(FaceConnect),          intent(inout) :: NewFaces(dimX:, 1:)
      real(dp),                   intent(inout) :: x(dimX:), q(1:, 1:)
      character(len =*),          intent(in)    :: TestName
      real(dp),                   intent(in)    :: dt, beta
      integer,                    intent(in)    :: t, nz
      type(column),      pointer                :: FineColumn, AdjColumn
      type(face),        pointer                :: MeshFace, MeshFaceNext
      type(ColumnStack)                         :: ColumnLeaves
      real(dp)                                  :: dA
      integer                                   :: Lscale
      integer                                   :: i, dim, indx
      dA = 0._dp; q(1:, 1:) = 0._dp
      ! Reorder the connectivity for faces
      ColumnLeaves = column_get_leaves(CoarseColumn)
      do while(associated(ColumnLeaves%head_))
        FineColumn => ColumnLeaves%head_
        ColumnLeaves%head_ => FineColumn%nextleaf_
        FineColumn%nextleaf_ => null()
        q(1:, 1:) = q(1:, 1:) + FineColumn%prop_%q_(1:, 1:)*column_dA(FineColumn)
        dA = dA + column_dA(FineColumn)
        ! Modify the faces for the coarsening
        ! if face at edges reorder the faces
        Lscale = 2**(FineColumn%l_ - CoarseColumn%l_)
        do dim = dimX, dimY
          do indx = 1, 2
            ! if the column is at the edge of the coarse column, we add the faces to the facelist
            if ( FineColumn%prop_%i_(dim) == Lscale*CoarseColumn%prop_%i_(dim) + (indx - 1)*(Lscale - 1) ) then
              ! point the adjacency to coarse column
              MeshFace => FineColumn%prop_%faces_(dim, indx)%head_
              do i = 1, FineColumn%prop_%faces_(dim, indx)%n_
                if (associated(MeshFace%pole_)) MeshFace%pole_%prop_%columns_(indx)%p_ => CoarseColumn
                MeshFace%prop_%columns_(3 - indx)%p_ => CoarseColumn
                MeshFace => MeshFace%nextc_(3 - indx)%p_
              enddo
              ! append a list to the list
              call FaceConnect_appendList(NewFaces(dim, indx), FineColumn%prop_%faces_(dim, indx), 3 - indx)
            else
              ! remove face
              ! pole faces do not presnt inner a coarse cell
              ! no need to consider it here
              MeshFace => FineColumn%prop_%faces_(dim, indx)%head_
              do while ( associated(MeshFace) )
                MeshFaceNext => MeshFace%nextc_(3 - indx)%p_
                ! get the adjacent column
                AdjColumn => MeshFace%prop_%columns_(indx)%p_
                ! remove the node from both columns and the face list
                call FaceConnect_remove(FineColumn%prop_%faces_(dim, indx), MeshFace, 3 - indx)
                call FaceConnect_remove(AdjColumn%prop_%faces_(dim, 3 - indx), MeshFace, indx)
                call FaceList_Remove(Faces(dim), MeshFace)
                MeshFace => MeshFaceNext
              enddo
            endif
            ! get the minimin x
            x(dimX) = min(x(dimX), FineColumn%prop_%x_(dimX))
            x(dimY) = min(x(dimY), FineColumn%prop_%x_(dimY))
          enddo
        enddo
      enddo
      q(1:, 1:) = q(1:, 1:)/dA
      call column_FaceLevelAdjust(CoarseColumn, Faces, NewFaces, nz, TestName, t, dt, beta)
    end subroutine column_FaceCoarsen

    subroutine column_FaceLevelAdjust(CoarseColumn, Faces, NewFaces, nz, TestName, t, dt, beta)
      use class_AMRTypes,                 only: FaceList, FaceConnect, face
      use class_face,                     only: FaceConnect_remove, FaceList_Remove, Face_Coarsen
      use class_parameter,                only: dimX, dimY
      use class_kind,                     only: dp
      type(column),   pointer, intent(inout) :: CoarseColumn
      type(FaceList),          intent(inout) :: Faces(dimX:)
      type(FaceConnect),       intent(inout) :: NewFaces(dimX:, 1:)
      character(len =*),       intent(in)    :: TestName
      real(dp),                intent(in)    :: dt, beta
      integer,                 intent(in)    :: t, nz
      integer                                :: dim, indx
      type(face),     pointer                :: MeshFace, MeshFaceNext
      type(column),   pointer                :: AdjColumn, OldAdjColumn
      ! Refine/Coarsen the faces
      do dim = dimX, dimY
        do indx = 1, 2
          MeshFace => NewFaces(dim, indx)%head_
          ! consider the situation around poles find the corresponding adjacent faces
          OldAdjColumn => null()
          do while( associated(MeshFace) )
            ! Get Adj Column
            MeshFaceNext => MeshFace%nextc_(3 - indx)%p_
            AdjColumn =>  MeshFace%prop_%columns_(indx)%p_
            if ( associated(OldAdjColumn, AdjColumn) ) then
              ! Remove the Face
              call FaceConnect_remove(NewFaces(dim, indx), MeshFace, 3 - indx)
              if (associated(MeshFace%pole_)) then
                call FaceConnect_remove(AdjColumn%prop_%faces_(dim, indx), MeshFace%pole_, 3 - indx)
                call FaceList_Remove(Faces(dim), MeshFace%pole_)
              else
                call FaceConnect_remove(AdjColumn%prop_%faces_(dim, 3 - indx), MeshFace, indx)
              endif
              call FaceList_Remove(Faces(dim), MeshFace)
              MeshFace => MeshFaceNext
              cycle
            endif
            ! adjust the faces
            call Face_Coarsen(MeshFace, min(0, max(CoarseColumn%l_, AdjColumn%l_) &
             - MeshFace%prop_%l_), nz, TestName, t, dt, beta, dim)
            OldAdjColumn => AdjColumn
            MeshFace => MeshFaceNext
          enddo
        enddo
      enddo
    end subroutine column_FaceLevelAdjust

    subroutine column_refinenode(self, Lscale, Children)
      type(column), pointer, intent(inout) :: self
      integer,               intent(in)    :: Lscale
      type(ColumnArray),     intent(inout) :: Children(0:, 0:)
      integer                              :: i, j, l, ii, jj
      integer                              :: LscaleTmp
      Children(0, 0)%p_ => self
      do l = 1, self%prop_%NumChange_
        LscaleTmp = 2**(self%prop_%NumChange_ - l + 1)
        do j = 0, Lscale, LscaleTmp
          do i = 0, Lscale, LscaleTmp
            if ( associated(Children(i, j)%p_%child_) ) print *, 'column_refine: program bug. Refining a non-existed cell.'
            allocate(Children(i, j)%p_%child_(0:1, 0:1))
            if (associated(Children(i, j)%p_%OtherMesh_%child_)) then
              do ii = 0, 1
                do jj = 0, 1
                  Children(i, j)%p_%OtherMesh_%child_(ii, jj)%OtherMesh_ => Children(i, j)%p_%child_(ii, jj)
                  Children(i, j)%p_%child_(ii, jj)%OtherMesh_ => Children(i, j)%p_%OtherMesh_%child_(ii, jj)
                enddo
              enddo
            else
              Children(i, j)%p_%child_(0, 0)%OtherMesh_ => Children(i, j)%p_%OtherMesh_
              Children(i, j)%p_%child_(1, 0)%OtherMesh_ => Children(i, j)%p_%OtherMesh_
              Children(i, j)%p_%child_(0, 1)%OtherMesh_ => Children(i, j)%p_%OtherMesh_
              Children(i, j)%p_%child_(1, 1)%OtherMesh_ => Children(i, j)%p_%OtherMesh_
            endif
            Children(i, j)%p_%child_(0, 0)%parent_ => Children(i, j)%p_
            Children(i, j)%p_%child_(0, 1)%parent_ => Children(i, j)%p_
            Children(i, j)%p_%child_(1, 0)%parent_ => Children(i, j)%p_
            Children(i, j)%p_%child_(1, 1)%parent_ => Children(i, j)%p_
            Children(i, j)%p_%child_(0:, 0:)%l_ = self%l_ + l

            Children(i + LscaleTmp/2, j)%p_               => Children(i, j)%p_%child_(1, 0)
            Children(i, j + LscaleTmp/2)%p_               => Children(i, j)%p_%child_(0, 1)
            Children(i + LscaleTmp/2, j + LscaleTmp/2)%p_ => Children(i, j)%p_%child_(1, 1)
            Children(i, j)%p_                             => Children(i, j)%p_%child_(0, 0)
          enddo
        enddo
      enddo

      do j = 0, Lscale
        do i = 0, Lscale
          if (Children(i, j)%p_%OtherMesh_%l_ >= Children(i, j)%p_%l_) &
          call column_refine_OtherMesh(Children(i, j)%p_%OtherMesh_, Children(i, j)%p_)
        enddo
      enddo
    end subroutine column_refinenode

    subroutine column_refine_OtherMesh(self, OtherMeshColumn)
      type(column), pointer,  intent(in)    :: self, OtherMeshColumn
      type(ColumnStack)                      :: tmp
      type(column), pointer                  :: FineColumn
      integer                                :: i, j

      tmp%head_ => self
      do while (associated(tmp%head_))
        FineColumn => tmp%head_
        tmp%head_ => FineColumn%nextStackleaf_
        FineColumn%nextStackleaf_ => null()
        FineColumn%OtherMesh_ => OtherMeshColumn
        if (associated(FineColumn%child_)) then
          do i = 0, 1
            do j = 0, 1
              FineColumn%child_(i, j)%nextStackleaf_ => tmp%head_
              tmp%head_ => FineColumn%Child_(i, j)
            enddo
          enddo
        endif
      enddo
    end subroutine column_refine_OtherMesh

    subroutine column_refineface(self, Lscale, i00, x0, dx, Faces, NewFaces, nz, TestName, dt, beta, t)
      use class_AMRTypes,                 only: face, FaceList, FaceConnect
      use class_face,                     only: Face_init, FaceList_append, FaceList_refine
      use class_face,                     only: FaceConnect_append, FaceConnect_refine, FaceConnect_remove
      use class_kind,                     only: dp
      use class_parameter,                only: dimX, dimY
      type(column),      intent(inout) :: self
      integer,           intent(in)    :: Lscale, i00(dimX:)
      real(dp),          intent(in)    :: x0(dimX:), dx(dimX:)
      type(FaceList),    intent(inout) :: Faces(dimX:)
      type(FaceConnect), intent(inout) :: NewFaces(0:, 0:, dimX:, 1:)
      character(len =*),   intent(in)  :: TestName
      real(dp),            intent(in)  :: dt, beta
      integer,             intent(in)  :: t, nz
      integer                          :: i, j, dim, indx, scale
      integer                          :: dL, nNewFaces, Ltarget, ii, jj
      type(face)                       :: NewFace
      type(face), pointer              :: MeshFace, MeshFaceNext, prev, next
      type(column), pointer            :: AdjColumn
      Ltarget = self%l_ + self%prop_%NumChange_
      ! setup the connecitivity for faces
      prev => faces(dimX)%tail_
      next => null()
      ! faces inside the columns
      do j = 0, Lscale
        do i = 0, Lscale - 1
          NewFace = face_init(prev, x0(dimX) + (i + 1)*dx(dimX), x0(dimY) + [j, j + 1]*dx(dimY), &
            i00(:) + [i + 1, j], Ltarget, nz, TestName, t, dt, beta, dimX, next)
          call FaceList_append(Faces(dimX), NewFace)
          prev => faces(dimX)%tail_
          call FaceConnect_append(Newfaces(i, j, dimX, 2), faces(dimX)%tail_, 1)
          call FaceConnect_append(Newfaces(i+1, j, dimX, 1), faces(dimX)%tail_, 2)
        enddo
      enddo
      prev => faces(dimY)%tail_
      do j = 0, Lscale - 1
        do i = 0, Lscale
          NewFace = face_init(prev, x0(dimY) + (j + 1)*dx(dimY), x0(dimX) + [i, i + 1]*dx(dimX), &
            i00(:) + [i, j + 1], Ltarget, nz, TestName, t, dt, beta, dimY, next)
          call FaceList_append(Faces(dimY), NewFace)
          prev => faces(dimY)%tail_
          call FaceConnect_append(Newfaces(i, j, dimY, 2), faces(dimY)%tail_, 1)
          call FaceConnect_append(Newfaces(i, j+1, dimY, 1), faces(dimY)%tail_, 2)
        enddo
      enddo

      ! faces at the boundaries
      do dim = dimX, dimY
        do indx = 1, 2
          ! get the head
          MeshFace => self%prop_%faces_(dim, indx)%head_
          do i = 1, self%prop_%faces_(dim, indx)%n_
            MeshFaceNext => MeshFace%nextc_(3 - indx)%p_
            ! remove the MeshFace from its current connectivity
            call FaceConnect_remove(self%prop_%faces_(dim, indx), MeshFace, 3 - indx)
            dL = Ltarget - MeshFace%prop_%l_
            nNewFaces = max(2**dL, 1)
            ! generate new faces
            call FaceList_refine(Faces(dim), MeshFace, dim, dL, nz, TestName, t, dt, beta)
            ! append new faces to the adjacent faceconnectivity list
            AdjColumn => MeshFace%prop_%columns_(indx)%p_
            if (associated(MeshFace%pole_)) then
              call FaceConnect_refine(AdjColumn%prop_%faces_(dim, indx), MeshFace%pole_, nNewFaces, 3 - indx)
            else
              call FaceConnect_refine(AdjColumn%prop_%faces_(dim, 3 - indx), MeshFace, nNewFaces, indx)
            endif
            do j = 1, nNewFaces
              ! get the corresponding columns of these faces
              MeshFace%prop_%columns_(indx)%p_ => AdjColumn
              if (associated(MeshFace%pole_)) MeshFace%pole_%prop_%columns_(3 - indx)%p_ => AdjColumn
              scale = 2**(MeshFace%prop_%l_ - Ltarget)
              ii = MeshFace%prop_%i_(dimX)/scale - i00(dimX) + (1 - indx)*(2 - dim)
              jj = MeshFace%prop_%i_(dimY)/scale - i00(dimY) + (1 - indx)*(dim - 1)
              if ( (dim == dimX) .and. (MeshFace%prop_%i_(dimX) == 0) .and. (indx == 2) ) ii = Lscale
              call FaceConnect_append(NewFaces(ii, jj, dim, indx), MeshFace, 3 - indx)
              MeshFace => MeshFace%next_
            enddo
            MeshFace => MeshFaceNext
          enddo
        enddo
      enddo
    end subroutine column_refineface

    ! get all active columns with the coarse column
    function column_get_leaves(self) result(val)
      type(column), target,  intent(in)     :: self
      type(ColumnStack)                      :: val
      type(ColumnStack)                      :: tmp
      type(column), pointer                  :: FineColumn
      integer                                :: i, j

      tmp%head_ => self
      do while (associated(tmp%head_))
        FineColumn => tmp%head_
        tmp%head_ => FineColumn%nextStackleaf_
        FineColumn%nextStackleaf_ => null()
        if (associated(FineColumn%child_)) then
          do i = 0, 1
            do j = 0, 1
              FineColumn%child_(i, j)%nextStackleaf_ => tmp%head_
              tmp%head_ => FineColumn%Child_(i, j)
            enddo
          enddo
        else
          FineColumn%nextleaf_ => val%head_
          val%head_ => FineColumn
        endif
      enddo
    end function column_get_leaves

    ! delete a column especially its properties
    subroutine column_delete(self)
      use class_parameter,          only: dimX, dimY
      use class_face,               only: FaceConnect_clean
      type(column), intent(inout)  :: self

      call FaceConnect_clean(self%prop_%faces_(dimX, 1), 2)
      call FaceConnect_clean(self%prop_%faces_(dimY, 1), 2)
      call FaceConnect_clean(self%prop_%faces_(dimX, 2), 1)
      call FaceConnect_clean(self%prop_%faces_(dimY, 2), 1)
      if (associated(self%prop_)) then
        deallocate(self%prop_)
      endif
    end subroutine column_delete

    ! get arbitrary column from root
    function ColumnList_get_column(self, i, j, l, l0) result(val)
      type(ColumnList),    intent(in)    :: self
      integer,              intent(in)    :: i, j, l, l0
      type(column), pointer               :: finer
      integer                             :: i_tmp, j_tmp, coeff

      type(column), pointer               :: Val
      coeff = 2**(l - l0)
      finer => self%TreeColumns_(i/coeff, j/coeff)%root_

      do while ( (finer%l_ < l) .and. ( associated(finer%child_) ) )
        coeff = 2**(l - finer%l_ - 1)
        i_tmp = i/coeff
        j_tmp = j/coeff

        i_tmp = i_tmp - i_tmp/2*2
        j_tmp = j_tmp - j_tmp/2*2
        finer => finer%child_(i_tmp, j_tmp)
      enddo
      Val => finer
    end function ColumnList_get_column

    subroutine ColumnList_init(self, coord_, nx, ny, ColumnFace)
      use class_coord,              only: coord
      use class_parameter,          only: dimX, dimY
      use class_AMRTypes,           only: FaceConnect
      type(ColumnList),  intent(inout) :: self
      type(FaceConnect), intent(in)    :: ColumnFace(0:, 0:, 1:, 1:)
      type(coord),       intent(in)    :: coord_
      integer,           intent(in)    :: nx, ny
      type(column)                     :: NewColumn
      integer                          :: i, j

      allocate(self%TreeColumns_(0:nx - 1, 0:ny-1))
      ! initialize the coordinate
      do j = 0, ny - 1
        do i = 0, nx - 1
          allocate(NewColumn%prop_)
          NewColumn%l_ = coord_%l0_
          NewColumn%prop_ = property_init([coord_%x_(i), coord_%y_(j)], &
            [coord_%x_(i + 1) - coord_%x_(i), coord_%y_(j + 1) - coord_%y_(j)], &
              [i, j], ColumnFace(i, j, :, :))
          allocate(self%TreeColumns_(i, j)%root_, source = NewColumn)
          call ColumnList_append(self, self%TreeColumns_(i, j)%root_)
          self%TreeColumns_(i, j)%root_%prop_%faces_(dimX, 1)%head_%prop_%columns_(2)%p_ => self%TreeColumns_(i, j)%root_
          self%TreeColumns_(i, j)%root_%prop_%faces_(dimX, 2)%head_%prop_%columns_(1)%p_ => self%TreeColumns_(i, j)%root_
          self%TreeColumns_(i, j)%root_%prop_%faces_(dimY, 1)%head_%prop_%columns_(2)%p_ => self%TreeColumns_(i, j)%root_
          self%TreeColumns_(i, j)%root_%prop_%faces_(dimY, 2)%head_%prop_%columns_(1)%p_ => self%TreeColumns_(i, j)%root_
          if (associated(self%TreeColumns_(i, j)%root_%prop_%faces_(dimY, 1)%head_%pole_)) &
            self%TreeColumns_(i, j)%root_%prop_%faces_(dimY, 1)%head_%pole_%prop_%columns_(1)%p_ => self%TreeColumns_(i, j)%root_
          if (associated(self%TreeColumns_(i, j)%root_%prop_%faces_(dimY, 2)%head_%pole_)) &
            self%TreeColumns_(i, j)%root_%prop_%faces_(dimY, 2)%head_%pole_%prop_%columns_(2)%p_ => self%TreeColumns_(i, j)%root_
        enddo
      enddo
    end subroutine ColumnList_init

    ! Append a column node into the end of the columnlist
    subroutine ColumnList_append(self, NewColumn)
      type(ColumnList),     intent(inout) :: self
      type(Column), pointer, intent(inout) :: NewColumn

      NewColumn%prev_ => self%tail_
      NewColumn%next_ => null()
      if (associated(self%head_)) then
        self%tail_%next_ => NewColumn
        self%tail_ => self%tail_%next_
      else
        self%head_ => NewColumn
        self%tail_ => self%head_
      endif
      self%n_ = self%n_ + 1
    end subroutine ColumnList_append

    ! insert a node following the previous node into the columnlist
    subroutine ColumnList_insert(self, NewColumn, PrevColumn)
      type(ColumnList),     intent(inout) :: self
      type(column), pointer, intent(inout) :: NewColumn
      type(column), pointer, intent(inout) :: PrevColumn

      if (.not. associated(PrevColumn) ) print *, "FaceList_insert: prev face is not present. programming bug"
      NewColumn%next_  => PrevColumn%next_
      NewColumn%prev_  => PrevColumn
      PrevColumn%next_ => NewColumn
      if (associated(NewColumn%next_)) then
        NewColumn%next_%prev_ => PrevColumn%next_
      else
        self%tail_ => PrevColumn%next_
      endif
      self%n_ = self%n_ + 1
    end subroutine ColumnList_insert

    ! remove column from columnlist
    subroutine ColumnList_remove(self, OldColumn)
      type(ColumnList),         intent(inout)  :: self
      type(column),     pointer, intent(inout)  :: OldColumn

      ! remove it from all FaceConnect as well
      if (.not. associated(OldColumn%prev_)) then
        self%head_ => OldColumn%next_
      else
        OldColumn%prev_%next_ => OldColumn%next_
      endif
      if (.not. associated(OldColumn%next_)) then
        self%tail_ => OldColumn%prev_
      else
        OldColumn%next_%prev_ => OldColumn%prev_
      endif

      call column_delete(OldColumn)
      nullify(OldColumn%next_)
      nullify(OldColumn%prev_)

      self%n_ = self%n_ - 1
    end subroutine ColumnList_remove

    ! delete all children of coarsecolumn from the tree and list
    subroutine ColumnList_delete_child(self, CoarseColumn)
      type(ColumnList)          , intent(inout) :: self
      type(column),      pointer, intent(inout) :: CoarseColumn
      integer                                   :: i, j
      type(column),      pointer                :: FineColumn, FinerColumn
      type(ColumnStack)                         :: stack
      logical                                   :: coarsen

      if (.not. associated(CoarseColumn%child_)) return
      stack%head_ => CoarseColumn
      do while( associated(Stack%head_) )
        ! using stack to avoid recursive procedure
        FineColumn  => Stack%head_
        if (associated(FineColumn%child_)) then
          coarsen = .true.
          do i = 0, 1
            do j = 0, 1
              FinerColumn => FineColumn%child_(i, j)
              if (associated(FinerColumn%child_)) then
                coarsen = .false.
                FinerColumn%nextStack_ => stack%head_
                stack%head_ => FinerColumn
              else
                if (associated(FinerColumn%prop_)) then
                  call columnlist_remove(self, FinerColumn)
                endif
              endif

            enddo
          enddo
          if (coarsen) then
            if (associated(FineColumn%OtherMesh_%child_)) &
              call column_Coarsen_OtherMesh(FineColumn%OtherMesh_, FineColumn) 
            deallocate(FineColumn%child_)
            stack%head_ => FineColumn%nextStack_
            FineColumn%nextStack_ => null()
          endif
        endif
      enddo
    end subroutine ColumnList_delete_child

    subroutine column_Coarsen_OtherMesh(self, OtherMeshColumn)
      type(column), pointer,  intent(in)    :: self, OtherMeshColumn
      type(ColumnStack)                      :: tmp
      type(column), pointer                  :: FineColumn
      integer                                :: i, j

      tmp%head_ => self
      do while (associated(tmp%head_))
        FineColumn => tmp%head_
        tmp%head_ => FineColumn%nextStackleaf_
        FineColumn%nextStackleaf_ => null()
        FineColumn%OtherMesh_ => OtherMeshColumn
        if (associated(FineColumn%child_)) then
          do i = 0, 1
            do j = 0, 1
              FineColumn%child_(i, j)%nextStackleaf_ => tmp%head_
              tmp%head_ => FineColumn%Child_(i, j)
            enddo
          enddo
        endif
      enddo
    end subroutine column_Coarsen_OtherMesh

    ! coarsen the FineColumn from the ColumnList
    subroutine ColumnList_Coarsen(self, Faces, FineColumn, nz, TestName, t, dt, beta)
      use class_AMRTypes,                only: FaceList
      use class_parameter,               only: dimX, dimY
      use class_kind,                    only: dp
      type(ColumnList),      intent(inout)  :: self
      type(column), pointer, intent(inout)  :: FineColumn
      type(FaceList),        intent(inout)  :: Faces(dimX:)
      integer,               intent(in)     :: nz
      character(len =*),     intent(in)     :: TestName
      real(dp),              intent(in)     :: dt, beta
      integer,               intent(in)     :: t
      type(column), pointer                 :: CoarseColumn, PrevColumn
      integer                               :: Lscale
      if (FineColumn%prop_%NumChange_ >= 0) return
      CoarseColumn => column_coarse(FineColumn, -FineColumn%prop_%NumChange_)
      allocate(CoarseColumn%prop_)
      Lscale = 2**(-FineColumn%prop_%NumChange_)
      CoarseColumn%prop_%dx_(:) = FineColumn%prop_%dx_(:)*Lscale
      CoarseColumn%prop_%i_(:)  = FineColumn%prop_%i_(:)/Lscale

      PrevColumn => FineColumn
      CoarseColumn%prop_%x_(dimX:) = FineColumn%prop_%x_(dimX:)
      CoarseColumn%prop_%qmid_(:, :, :) = 0._dp
      CoarseColumn%prop_%div_(:) = 0._dp

      call column_FaceCoarsen(CoarseColumn, Faces, CoarseColumn%prop_%Faces_, &
        CoarseColumn%prop_%x_, CoarseColumn%prop_%q_, nz, TestName, t, dt, beta)
      CoarseColumn%prop_%cosp_  = (sin(CoarseColumn%prop_%x_(dimY) + CoarseColumn%prop_%dx_(dimY)) - &
        sin(CoarseColumn%prop_%x_(dimY)))/CoarseColumn%prop_%dx_(dimY)      

      call ColumnList_Insert(self, CoarseColumn, FineColumn)
      call ColumnList_delete_child(self, CoarseColumn)
      FineColumn => CoarseColumn
    end subroutine ColumnList_Coarsen

    ! refine the coarsecolumn
    subroutine ColumnList_refine(self, CoarseColumn, Faces, nz, TestName, dt, beta, t)
      use class_kind,                          only: dp
      use class_AMRTypes,                      only: face, FaceList, FaceConnect
      use class_face,                          only: face_init
      use class_parameter,                     only: dimX, dimY
      type(ColumnList),           intent(inout)   :: self
      type(column),      pointer, intent(inout)   :: CoarseColumn
      type(FaceList),             intent(inout)   :: Faces(dimX:)
      integer,                    intent(in)      :: nz
      character(len =*),              intent(in)  :: TestName
      real(dp),                       intent(in)  :: dt, beta
      integer,                        intent(in)  :: t
      type(column),      pointer                  :: prev
      integer                                     :: i, j, k, dim, indx
      integer                                     :: Lscale, LscaleTmp
      integer                                     :: i00(dimX:dimY)
      real(dp)                                    :: x0(dimX:dimY), dx(dimX:dimY)
      type(FaceConnect) :: NewFaces(0:2**CoarseColumn%prop_%NumChange_ - 1, 0:2**CoarseColumn%prop_%NumChange_ - 1, dimX:dimY, 2)
      type(ColumnArray) :: Children(0:2**CoarseColumn%prop_%NumChange_ - 1, 0:2**CoarseColumn%prop_%NumChange_- 1)
      type(face),        pointer                  :: MeshFace

      Lscale     = 2**CoarseColumn%prop_%NumChange_ - 1
      LscaleTmp  = 2**CoarseColumn%prop_%NumChange_
      i00(dimX:) = LscaleTmp*CoarseColumn%prop_%i_(dimX:)
      x0(dimX:)  = CoarseColumn%prop_%x_(dimX:)
      dx(dimX:)  = CoarseColumn%prop_%dx_/LscaleTmp
      call column_refineface(CoarseColumn, Lscale, i00, x0, dx, Faces, NewFaces, nz, TestName, dt, beta, t)
      ! create refined high resolution tree nodes
      call column_refineNode(CoarseColumn, Lscale, Children)
      ! setup the properties of columns
      prev => CoarseColumn
      do j = 0, Lscale
        do i = 0, Lscale
          do dim = dimX, dimY
            do indx = 1, 2
              MeshFace => NewFaces(i, j, dim, indx)%head_
              do k = 1, NewFaces(i, j, dim, indx)%n_

                MeshFace%prop_%columns_(3 - indx)%p_ => Children(i, j)%p_
                if (associated(MeshFace%pole_)) MeshFace%pole_%prop_%columns_(indx)%p_ => Children(i, j)%p_
                MeshFace => MeshFace%nextc_(3 - indx)%p_
              enddo
            enddo
          enddo
          allocate(Children(i, j)%p_%prop_)
          
          Children(i, j)%p_%prop_ = property_init(x0(dimX:) + dx(dimX:)*[i, j], &
            dx, i00(dimX:) + [i, j], NewFaces(i, j, 1:, 1:))
          Children(i, j)%p_%prop_%NumChange_ = 0
          call ColumnList_insert(self, Children(i, j)%p_, prev)
          prev => prev%next_
        enddo
      enddo
      call ColumnList_remove(self, CoarseColumn)
    end subroutine ColumnList_refine

    subroutine ColumnList_AMR_Mark(self, Faces, l0, MaxRef, ntracer, nx, ny, nz, &
      theta_r, theta_c, dt, x, y, dxMax, beta, t, TestName, is_init, icycle)
      use class_kind,      only: dp
      use class_AMRTypes,  only: FaceList
      use class_parameter, only: dimX, dimY
      type(ColumnList),          intent(inout) :: self
      type(FaceList),            intent(inout) :: Faces(dimX:)
      integer,                   intent(in)    :: l0, MaxRef, ntracer, nx, ny, nz, t
      real(dp),                  intent(in)    :: theta_r, theta_c, x(0:), y(0:), dxMax(dimX:), dt, beta
      character(len =*),         intent(in)    :: TestName
      logical,                   intent(in)    :: is_init
      integer,                   intent(in)    :: icycle
      type(column),     pointer                :: MeshColumn, MeshColumnNext, CoarseColumn, LeafColumn
      type(ColumnStack)                        :: ColumnLeaves
      integer                                  :: i, l, MaxLevel, nref
      logical                                  :: coarsen
      real(dp) :: mintheta, maxtheta

      if (icycle == 0) then
        MeshColumn => self%head_
        do while(associated(MeshColumn))
          call column_criteria_grad(MeshColumn, ntracer, nz)
          MeshColumnNext => MeshColumn%next_
          if ((MeshColumn%prop_%theta_ < theta_c) .and. (MeshColumn%l_ > l0)) then
            coarsen = .true.
            CoarseColumn => MeshColumn
            do l = 1, MeshColumn%l_ - l0
              CoarseColumn => CoarseColumn%parent_
              ColumnLeaves = column_get_leaves(CoarseColumn)
              LeafColumn => ColumnLeaves%head_
              do while (associated(LeafColumn))
                if (LeafColumn%prop_%theta_ > theta_c) Coarsen = .false.
                LeafColumn => LeafColumn%nextleaf_
              enddo
              if (.not. coarsen) exit
              do while (associated(ColumnLeaves%head_ ))
                LeafColumn                 => ColumnLeaves%head_
                ColumnLeaves%head_         => LeafColumn%nextleaf_
                LeafColumn%nextleaf_       => null()
                if (Coarsen) LeafColumn%prop_%NumChange_ = -l
              enddo
            enddo
            call ColumnList_Coarsen(self, Faces, MeshColumn, nz, TestName, t, dt, beta)
            MeshColumn%prop_%NumChange_ = 0
            MeshColumnNext => MeshColumn%next_
          endif
          MeshColumn => MeshColumnNext
        enddo
      endif
      mintheta = 99999999.
      maxtheta = 0._dp
      nref = 0
      MaxLevel = l0 + MaxRef
      MeshColumn => self%head_
      do i = 1, self%n_
        call column_criteria_grad(MeshColumn, ntracer, nz)
          mintheta = min(mintheta, MeshColumn%prop_%theta_)
          maxtheta = max(maxtheta, MeshColumn%prop_%theta_)        
        if (MeshColumn%prop_%theta_ > theta_r) then
          MeshColumn%prop_%NumChange_ = &
            max(MeshColumn%prop_%NumChange_, MaxLevel - MeshColumn%l_)
          nref = nref + 4
          if (.not. is_init) then
            ! refine the departure cells (and change its theta)
            call ColumnList_depart(self, MeshColumn, x, nx, dxMax(dimX), dt, l0, MaxRef, dimX)
            call ColumnList_depart(self, MeshColumn, y, ny, dxMax(dimY), dt, l0, MaxRef, dimY)
          endif
        endif
        call column_ratio(MeshColumn)
        MeshColumn => MeshColumn%next_
      enddo
    end subroutine ColumnList_AMR_Mark

    subroutine ColumnList_depart(columns, MeshColumn, x, nx, dxMax, dt, l0, MaxRef, dim)
      use class_kind,      only: dp
      use class_AMRTypes,  only: face
      use class_parameter, only: dimX, dimY, ae, pi
      type(ColumnList),           intent(in) :: columns
      type(column),      pointer, intent(in) :: MeshColumn
      real(dp)                  , intent(in) :: x(0:), dxMax, dt
      integer                   , intent(in) :: l0, MaxRef, dim, nx
      type(face),        pointer             :: MeshFace
      type(column),      pointer             :: DepartColumn, LeafColumn
      type(ColumnStack)                      :: ColumnLeaves
      real(dp)                               :: xd, xdtmp, dx, cosp
      integer                                :: i, n, idx(2, dimX:dimY), Lscale, indx
      do indx = 1, 2
        xd = merge(2*pi, -2._dp*pi, indx == 1)

        MeshFace => MeshColumn%prop_%faces_(dim, indx)%head_
        do i = 1, MeshColumn%prop_%faces_(dim, indx)%n_
          if (dim == dimX) then
            cosp = (sin(MeshFace%prop_%y_(2)) - sin(MeshFace%prop_%y_(1)))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
            if (indx == 1) xdtmp = minval(MeshFace%prop_%x_ - MeshFace%prop_%u_(:)*dt/cosp/ae)
            if (indx == 2) xdtmp = maxval(MeshFace%prop_%x_ - MeshFace%prop_%u_(:)*dt/cosp/ae)
            xdtmp = modulo(xdtmp, 2._dp*pi)
            if (xdtmp >= x(nx)) xdtmp = xdtmp - 2._dp*pi
          else
            if (indx == 1) xdtmp = minval(sin(MeshFace%prop_%x_) - MeshFace%prop_%u_(:)*cos(MeshFace%prop_%x_)*dt/ae)
            if (indx == 2) xdtmp = maxval(sin(MeshFace%prop_%x_) - MeshFace%prop_%u_(:)*cos(MeshFace%prop_%x_)*dt/ae)
            xdtmp = asin(max(min(1._dp, xdtmp), -1._dp))
          endif
          xd = merge(min(xd, xdtmp), max(xd, xdtmp), indx == 1)
          MeshFace => MeshFace%nextc_(3 - indx)%p_
        enddo
        ! find the idx
        idx(indx, dim) = INT((xd - x(0))/dxMax)
        do while (xd > x(idx(indx, dim) + 1) )
          idx(indx, dim) = idx(indx, dim) + 1
        end do

        Lscale = 2**(MeshColumn%l_ - l0)
        dx = (x(idx(indx, dim) + 1) - x(idx(indx, dim)))/Lscale
        n = 1
        do while (xd > x(idx(indx, dim)) + n*dx)
          n = n + 1
        end do
        idx(indx, dim) = idx(indx, dim)*Lscale + n - 1
        idx(indx, 3 - dim) = MeshColumn%prop_%i_(3 - dim)
      enddo

      do i = 1, 2
        DepartColumn => ColumnList_get_column(columns, idx(i, dimX), idx(i, dimY), MeshColumn%l_, l0)
        if (.not. associated(DepartColumn%child_)) then
          DepartColumn%prop_%NumChange_ = l0 + MaxRef - DepartColumn%l_
          call column_ratio(DepartColumn)
        else
          ColumnLeaves = column_get_leaves(DepartColumn)
          do while (associated(ColumnLeaves%head_))
            LeafColumn => ColumnLeaves%head_
            ColumnLeaves%head_ => LeafColumn%nextleaf_
            LeafColumn%nextleaf_ => null()
            LeafColumn%prop_%NumChange_ = l0 + MaxRef - LeafColumn%l_
            call column_ratio(LeafColumn)
          enddo
        endif
      enddo
    end subroutine ColumnList_depart

    subroutine Column_stability(MeshColumn, dt)
      use class_AMRTypes,                 only: face
      use class_kind,                     only: dp
      use class_parameter,                only: dimX, dimY, ae
      type(column), pointer, intent(inout) :: MeshColumn
      real(dp)             , intent(in)    :: dt
      type(face),   pointer                :: MeshFace
      integer                              :: dim, j
      real(dp)                             :: u(dimX:dimY, 2), y(2), cosp_, xd1, xd2
      ! use leaf as temporary stack
      u(:, 1) = 99999._dp
      u(:, 2) = -99999._dp
      do dim = dimX, dimY
        if (MeshColumn%prop_%faces_(dim, 1)%n_ == &
          MeshColumn%prop_%faces_(dim, 2)%n_) then
          u(dim, :) = 0._dp
          cycle
        endif
        MeshFace => MeshColumn%prop_%faces_(dim, 1)%head_
        do j = 1, MeshColumn%prop_%faces_(dim, 1)%n_
          cosp_ = 1._dp
          if (dim == dimX) cosp_ = (sin(MeshFace%prop_%y_(2)) - &
            sin(MeshFace%prop_%y_(1)))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
          u(dim, 1) = min(u(dim, 1), minval(MeshFace%prop_%u_(1:)/cosp_))
          MeshFace => MeshFace%nextc_(2)%p_
        enddo
        MeshFace => MeshColumn%prop_%faces_(dim, 2)%head_
        do j = 1, MeshColumn%prop_%faces_(dim, 2)%n_
          cosp_ = 1._dp
          if (dim == dimX) cosp_ = (sin(MeshFace%prop_%y_(2)) - &
            sin(MeshFace%prop_%y_(1)))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
          u(dim, 2) = max(u(dim, 2), maxval(MeshFace%prop_%u_(1:)/cosp_))
          MeshFace => MeshFace%nextc_(1)%p_
        enddo
      enddo
      y(1) = MeshColumn%prop_%x_(dimY)
      y(2) = MeshColumn%prop_%dx_(dimY) + MeshColumn%prop_%x_(dimY)
      xd2  = asin(max(min(1._dp, sin(y(2)) - u(dimY, 2)*dt*cos(y(2))), -1._dp))
      xd1  = asin(max(min(1._dp, sin(y(1)) - u(dimY, 1)*dt*cos(y(1))), -1._dp))
      if ( (MeshColumn%prop_%dx_(dimX) < (u(dimX, 2) - u(dimX, 1))*dt/ae) .or. &
        (xd2 < xd1) ) then
        MeshColumn%prop_%NumChange_ = max(1, MeshColumn%prop_%NumChange_)
        call column_ratio(MeshColumn)
      endif
    end subroutine Column_stability

    ! get refinement criteria
    subroutine column_criteria_grad(self, ntracer, nz)
      use class_kind,                 only: dp
      use class_parameter,            only: dimX, dimY
      type(column), pointer, intent(inout) :: self
      integer,               intent(in)    :: ntracer, nz
      type(column), pointer                :: AdjColumn, LeafColumn
      type(ColumnStack)                    :: ColumnLeaves
      integer                              :: dim, indx
      integer                              :: Lscale, coeff
      real(dp)                             :: dA, dx, q(ntracer + 1, nz)
      real(dp)                             :: grad(dimX:dimY, 2, ntracer + 1, nz)

      ! get the gradient
      do dim = dimX, dimY
        do indx = 1, 2
          AdjColumn => self%prop_%faces_(dim, indx)%head_%prop_%columns_(indx)%p_
          if (AdjColumn%l_ < self%l_) then
            ColumnLeaves = column_get_leaves(column_coarse(self, self%l_ - AdjColumn%l_))
          else
            ColumnLeaves = column_get_leaves(column_coarse(AdjColumn, AdjColumn%l_ - self%l_))
          endif
          dA = 0._dp; q(1:, 1:) = 0._dp;
          do while (associated(ColumnLeaves%head_))
            LeafColumn => ColumnLeaves%head_
            ColumnLeaves%head_ => LeafColumn%nextleaf_
            LeafColumn%nextleaf_ => null()
            q(1:, 1:) = q(1:, 1:) + LeafColumn%prop_%q_(1:, 1:)*column_dA(LeafColumn)
            dA = dA + column_dA(LeafColumn)
          enddo
          q(1:, 1:) = q(1:, 1:)/dA
          Lscale = 2**abs(self%l_ - AdjColumn%l_)
          if (AdjColumn%l_ < self%l_) then
            dx = (self%prop_%dx_(dim)*Lscale + AdjColumn%prop_%dx_(dim))/2._dp
            if (dim == dimX) dx = dx/AdjColumn%prop_%cosp_
            coeff = 3 - 2*indx
            grad(dim, indx, 1:, 1:) = coeff*(q(1:, 1:) - AdjColumn%prop_%q_(1:, 1:))/dx
          else
            dx = (self%prop_%dx_(dim) + AdjColumn%prop_%dx_(dim)*Lscale)/2.
            if (dim == dimX) dx = dx/Self%prop_%cosp_
            coeff = 3 - 2*indx
            grad(dim, indx, 1:, 1:) = coeff*(self%prop_%q_(1:, 1:) - q(1:, 1:))/dx
          endif
        enddo
      enddo
      self%prop_%theta_ = maxval(abs(grad))
      ! if (self%prop_%i_(dimX) == 10) self%prop_%theta_ = 1.
    end subroutine column_criteria_grad

    subroutine column_ratio(MeshColumnIn)
      use class_AMRTypes,               only: face
      use class_parameter,              only: dimX, dimY
      type(column), pointer, intent(inout) :: MeshColumnIn
      type(column), pointer                :: MeshColumn, AdjColumn
      type(face),   pointer                :: MeshFace
      type(ColumnStack)                    :: UnChecked
      integer                              :: TargetLevel, AdjTargetLevel
      integer                              :: dim, indx, i
      UnChecked%head_ => MeshColumnIn
      do while(associated(UnChecked%head_))
        MeshColumn => UnChecked%head_
        UnChecked%head_ => MeshColumn%nextStack_
        MeshColumn%nextStack_ => null()
        ! target level of the current column
        TargetLevel = MeshColumn%l_ + MeshColumn%prop_%NumChange_

        do dim = dimX, dimY
          do indx = 1, 2
            ! if the
            MeshFace => MeshColumn%prop_%faces_(dim, indx)%head_
            do i = 1, MeshColumn%prop_%faces_(dim, indx)%n_
              AdjColumn => MeshFace%prop_%columns_(indx)%p_
              AdjTargetLevel = AdjColumn%l_ + AdjColumn%prop_%NumChange_
              if (TargetLevel - AdjTargetLevel > 1) then
                ! get all leaves
                AdjColumn%prop_%NumChange_ = TargetLevel - 1 - AdjColumn%l_
                AdjColumn%nextStack_ => UnChecked%head_
                UnChecked%head_ => AdjColumn
              endif
              MeshFace => MeshFace%nextc_(3 - indx)%p_
            enddo
          enddo
        enddo
      enddo
    end subroutine column_ratio

    subroutine ColumnList_q(self)
      use class_kind, only: dp
      type(ColumnList),          intent(inout) :: self
      type(column),     pointer                :: MeshColumn, OtherMeshColumn
      type(ColumnStack)                        :: ColumnLeaves
      real(dp)                                 :: dA
      integer                                  :: i
      MeshColumn => self%head_
      do i = 1, self%n_
        ColumnLeaves = column_get_leaves(MeshColumn%OtherMesh_)
        dA = 0._dp
        MeshColumn%prop_%q0_(1:, 1:) = 0._dp
        do while (associated(ColumnLeaves%head_))
          OtherMeshColumn           => ColumnLeaves%head_
          ColumnLeaves%head_        => OtherMeshColumn%nextleaf_
          OtherMeshColumn%nextleaf_ => null()
          MeshColumn%prop_%q0_(1:, 1:) = MeshColumn%prop_%q0_(1:, 1:) + OtherMeshColumn%prop_%q0_(1:, 1:)*column_dA(OtherMeshColumn)
          dA = dA + column_dA(OtherMeshColumn)
        enddo
        MeshColumn%prop_%q0_(1:, 1:) = MeshColumn%prop_%q0_(1:, 1:)/dA

        MeshColumn => MeshColumn%next_
      enddo
    end subroutine ColumnList_q
end module class_column
