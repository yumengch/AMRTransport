module class_AMRTypes
  use class_kind,      only: dp
  use class_parameter, only: dimX, dimY, nz, ntracer, Inner, Outer
  implicit none

  ! array of pointers for column object
  type :: ColumnArray
    type(column), pointer :: p_ => null()
  end type ColumnArray

  ! properties of faces
  type :: fproperty
    real(dp)              :: x_           ! position of the face
    real(dp)              :: y_(2)        ! position of the two vertiex of the face
    integer               :: i_(2), l_    ! index of the face
    type(ColumnArray)     :: columns_(2)  ! adjacent columns of the face
    real(dp)              :: u_(nz)        ! velocity at the faces with size u(nz)
    type(ColumnArray)     :: DepartColumns_(Inner:Outer)
    real(dp)              :: Departx0_(dimX:dimY), Departdx_(dimX:dimY)
    real(dp)              :: xd_
  end type fproperty

  ! array of pointers for column object
  type FaceArray
    type(Face), pointer   :: p_ => null()
  end type FaceArray

  ! class face
  type face
    ! face list node
    type(face),      pointer :: prev_ => null()
    type(face),      pointer :: next_ => null()
    ! face connect node
    type(FaceArray)          :: prevc_(2)
    type(FaceArray)          :: nextc_(2)
    ! face pole node
    type(face),      pointer :: pole_ => null()
    ! property of the face
    type(fproperty), pointer :: prop_ => null()
  end type face

  ! list of faces
  type FaceList
    type(face), pointer :: head_ => null()
    type(face), pointer :: tail_ => null()
    integer             :: n_ = 0
  end type FaceList

  ! pointer to the faces in a column
  type FaceConnect
    type(face), pointer :: head_ => null()
    type(face), pointer :: tail_ => null()
    integer             :: n_ = 0
  end type FaceConnect

  ! class for lagrangian volume
  type LV
    real(dp)              :: x_(2)
    integer               :: jLimit, lLimit
    type(column), pointer :: ArrivalColumn => null()
  end type LV

  type cproperty
    ! index, connectivity, geometry
    integer                     :: i_(dimX:dimY)                ! index of the column
    real(dp)                    :: x_(dimX:dimY)                ! position of the minimum vertex
    real(dp)                    :: dx_(dimX:dimY)               ! size of the column
    real(dp)                    :: cosp_                        ! scale factor for latitude size
    type(FaceConnect)           :: faces_(dimX:dimY, 2)         ! Faces adajcent to the column
    ! old and new tracer concentration
    real(dp)                    :: q0_(ntracer + 1, nz)                    ! tracer distribution in previous time step q0(ntracer + 1, nz)
    real(dp)                    :: qinter_(ntracer + 1, nz, dimX:dimY)             ! tracer distribution for the intermediate steps
    real(dp)                    :: q_(ntracer + 1, nz)                     ! tracer distribution in the new time step q(ntracer + 1, nz)
    ! properties for the scheme
    real(dp)                    :: div_(dimX:dimY)              ! divergence of the cell in each dimension
    real(dp)                    :: qmid_(ntracer + 1, dimX:dimY, 0:Outer)               ! tracer distribution of intermiediate step qmid(ntracer + 1, dimX:dimY, Inner:Outer)
    ! values for cell integrated semi-Lagrangian scheme
    type(LV)                    :: LV_(dimX:dimY, 2)                       ! Lagrangian volume for the computation in CISL
    real(dp)                    :: qx_(dimX:dimY, ntracer + 1, -2:2)   ! value for interpolations size: q0(dimX:dimY, ntracer + 1, -2:2)
    real(dp)                    :: qdx_(dimX:dimY, -2:2)               ! distance for interpolations
    integer                     :: nGhost_   = 0                       ! Number of Ghost Cells ! if nGhost > 1 there should be warnings
    integer                     :: dimStart_ = 0                       ! start dim
    integer                     :: dimEnd_   = 0                       ! end dim
    integer                     :: StencilStart_(dimX:dimY) = 0        ! start dim and stencil
    integer                     :: StencilEnd_(dimX:dimY) = 0          ! stencil , determines the stencil of the column (PPM2D, 1D or no PPM)
    real(dp)                    :: a_(dimX:dimY, ntracer + 1), b_(dimX:dimY, ntracer + 1)  ! coefficient of PPM scheme size: a(dimX:dimY, ntracer + 1), b(dimX:dimY, ntracer + 1)
    real(dp), allocatable       :: qghost_(:, :, :)                    ! values for the ghost cell size of (0:2**nGhost, 0:2**nGhost, ntracer + 1)
    logical                     :: has_dim(dimX:dimY) = .false.
    ! refinement
    integer                     :: NumChange_ = 0               ! Number of Refinement
    real(dp)                    :: theta_     = 0
    real(dp)                    :: area(2)       = 1._dp
  end type cproperty

  type :: column
    type(column),    pointer  :: child_(:, :) => null()         ! children of the column
    type(column),    pointer  :: parent_      => null()         ! coarser column
    type(column),    pointer  :: next_        => null()         ! next in the columnlist
    type(column),    pointer  :: prev_        => null()         ! prev in the columnList
    type(column),    pointer  :: nextStack_   => null()         ! for temporary stack of columns
    type(column),    pointer  :: nextStackLeaf_ => null()       ! for temporary stack for getting leaves
    type(column),    pointer  :: nextleaf_    => null()         ! for temporary leaves of columns
    type(column),    pointer  :: OtherMesh_   => null()         ! pointer to the corresponding column in the other mesh
    type(column),    pointer  :: nextLV_      => null()
    type(cproperty), pointer  :: prop_        => null()

    integer                   :: l_                             ! level of the column
  end type column

  type :: tree
    type(column), pointer   :: root_ => null()
  end type tree

  type :: ColumnStack
    type(column), pointer   :: head_ => null()
  end type ColumnStack

  type ColumnList
    type(tree), allocatable :: TreeColumns_(:, :)
    type(column), pointer   :: head_ => null()
    type(column), pointer   :: tail_ => null()
    integer                 :: n_ = 0
  end type ColumnList
end module class_AMRTypes