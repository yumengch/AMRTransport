module class_face
use class_kind,     only: dp
use class_AMRTypes, only: face, FaceList, FaceConnect
implicit none
private
public :: Face_Coarsen, face_init
public :: FaceConnect_remove, FaceConnect_append, FaceConnect_appendList, FaceConnect_clean, FaceConnect_refine
public :: FaceList_Remove, FaceList_append, FaceList_refine, FaceList_clean
contains
  function Face_init(prev, x, y, i, l, x0, y0, l0, nx, ny, nz, uInput, vInput, dim, next) result(val)
    use class_setupwind,          only: setup_interp
    real(dp),            intent(in)  :: x
    real(dp),            intent(in)  :: y(1:)
    real(dp),            intent(in)  :: x0(0:), y0(0:)
    real(dp),            intent(in)  :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)
    integer,             intent(in)  :: l0, nx, ny
    integer,             intent(in)  :: i(1:), l, nz
    type(face), pointer, intent(in)  :: prev, next
    integer,             intent(in)  :: dim
    type(face)                       :: val
    allocate(val%prop_)
    val%prop_%x_ = x
    val%prop_%y_ = y(:)
    val%prop_%i_(:) = i;
    val%prop_%l_ = l
    val%next_ => next
    val%prev_ => prev
    val%prop_%u_(1:) = 0._dp
    ! give velocity to the face
    ! call setup_ideal(val, TestName, nz, t, dt, beta, dim)
    call setup_interp(val, x0, y0, l0, nx, ny, nz, uInput, vInput, dim)
  end function Face_init

  subroutine Face_coarsen(self, l, x, y, l0, nx, ny, nz, uInput, vInput, dim)
    use class_setupwind,          only: setup_interp
    type(face), intent(inout)        :: self
    integer,    intent(in)           :: l, nz
    real(dp),            intent(in)  :: x(0:), y(0:)
    real(dp),            intent(in)  :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)
    integer,             intent(in)  :: l0, nx, ny
    integer,             intent(in)  :: dim
    integer                          :: Lscale
    real(dp)                         :: dy
    if (l == 0) return
    Lscale = 2**(-l)
    dy = Lscale*(self%prop_%y_(2) - self%prop_%y_(1))
    self%prop_%y_(2)  = self%prop_%y_(1) + dy
    self%prop_%i_(1:) = self%prop_%i_(1:)/Lscale
    self%prop_%l_     = self%prop_%l_ + l
    ! give velocity to the face
    ! call setup_ideal(self, TestName, nz, t, dt, beta, dim)
    call setup_interp(self, x, y, l0, nx, ny, nz, uInput, vInput, dim)
    if (associated(self%pole_)) then
      self%pole_%prop_%y_(2)  = self%pole_%prop_%y_(1) + dy
      self%pole_%prop_%i_(1:) = self%pole_%prop_%i_(1:)/Lscale
      self%pole_%prop_%l_     = self%pole_%prop_%l_ + l
      ! call setup_ideal(self%pole_, TestName, nz, t, dt, beta, dim)
      call setup_interp(self%pole_, x, y, l0, nx, ny, nz, uInput, vInput, dim)
    endif
  end subroutine Face_coarsen

  ! Add Face to the Face List
  subroutine FaceList_append(self, NewFace)
    class(FaceList),    intent(inout) :: self
    type(face),         intent(inout) :: NewFace

    NewFace%prev_ => self%tail_
    NewFace%next_ => null()
    if (associated(self%head_)) then
      allocate(self%tail_%next_, source = NewFace)
      self%tail_ => self%tail_%next_
    else
      allocate( self%head_, source =  NewFace )
      self%tail_ => self%head_
    endif
    self%n_ = self%n_ + 1
  end subroutine FaceList_append

  ! Insert element to this list if we know its prev position.
  subroutine FaceList_insert(self, NewFace, PrevFace)
    class(FaceList),     intent(inout) :: self
    type(face),          intent(inout) :: NewFace
    type(face), pointer, intent(inout) :: PrevFace

    if (.not. associated(PrevFace) ) print *, "FaceList_insert: prev face is not present. programming bug"
    NewFace%next_ => PrevFace%next_
    NewFace%prev_ => PrevFace
    allocate(PrevFace%next_, source = NewFace)
    if (associated(NewFace%next_)) then
      NewFace%next_%prev_ => PrevFace%next_
    else
      self%tail_ => PrevFace%next_
    endif
    self%n_ = self%n_ + 1
  end subroutine FaceList_insert

  ! refine the face based on
  subroutine FaceList_refine(self, MeshFace, dim, l, x, y, l0, nx, ny, nz, uInput, vInput)
    use class_setupwind,                 only: setup_interp
    type(FaceList),          intent(inout)  :: self
    type(Face),     pointer, intent(inout)  :: MeshFace
    integer,                 intent(in)     :: dim, l, nz
    real(dp),                intent(in)     :: x(0:), y(0:)
    real(dp),                intent(in)     :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)
    integer,                 intent(in)     :: l0, nx, ny
    integer                                 :: Lscale, i
    real(dp)                                :: dy
    type(face)                              :: NewFace
    type(face),     pointer                 :: prev, next, poleFace
    if (l < 0) return
    Lscale = 2**l
    ! get the refined levels
    dy = ( MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1) )/Lscale
    ! change the property of the original face
    MeshFace%prop_%y_(2)  = MeshFace%prop_%y_(1) + dy
    MeshFace%prop_%i_(1:) = Lscale*MeshFace%prop_%i_(1:)
    MeshFace%prop_%l_     = MeshFace%prop_%l_ + l
    !call setup_ideal(MeshFace, TestName, nz, t, dt, beta, dim)
    call setup_interp(MeshFace, x, y, l0, nx, ny, nz, uInput, vInput, dim)
    ! add the rest of the higher level faces into the FaceList
    prev => MeshFace
    next => MeshFace%next_
    do i = 1, Lscale - 1
      NewFace = face_init(prev, MeshFace%prop_%x_, MeshFace%prop_%y_(1) + [i, i + 1]*dy, &
        MeshFace%prop_%i_(1:) + [i*(dim - 1), i*(2 - dim)], MeshFace%prop_%l_, x, y, l0, nx, ny, &
        nz, uInput, vInput, dim, next)
      call FaceList_insert(self, NewFace, prev)
      prev => prev%next_
    enddo
    ! create corresponnding pole faces
    if (associated(MeshFace%pole_)) then
      MeshFace%pole_%prop_%y_(2)  = MeshFace%pole_%prop_%y_(1) + dy
      MeshFace%pole_%prop_%i_(1:) = Lscale*MeshFace%pole_%prop_%i_(1:)
      MeshFace%pole_%prop_%l_     = MeshFace%prop_%l_
      ! call setup_ideal(MeshFace%pole_, TestName, nz, t, dt, beta, dim)
      call setup_interp(MeshFace%pole_, x, y, l0, nx, ny, nz, uInput, vInput, dim)
      prev     => MeshFace%pole_
      next     => MeshFace%pole_%next_
      poleFace => MeshFace%next_
      do i = 1, Lscale - 1
        NewFace = face_init(prev, MeshFace%pole_%prop_%x_, MeshFace%pole_%prop_%y_(1) + &
          [i, i + 1]*dy, MeshFace%pole_%prop_%i_(1:) + [i*(dim - 1), i*(2 - dim)], &
          MeshFace%pole_%prop_%l_, x, y, l0, nx, ny, &
          nz, uInput, vInput, dim, next)
        NewFace%pole_  => poleFace
        call FaceList_insert(self, NewFace, prev)
        prev           => prev%next_
        poleFace%pole_ => prev
        poleFace       => poleFace%next_
      enddo
    endif
  end subroutine FaceList_refine

  ! Remove Face.
  subroutine FaceList_remove(self, OldFace)
    class(FaceList),         intent(inout)  :: self
    type(face),     pointer, intent(inout)  :: OldFace

    ! remove it from all FaceConnect as well
    if (.not. associated(OldFace%prev_)) then
      self%head_ => OldFace%next_
    else
      OldFace%prev_%next_ => OldFace%next_
    endif
    if (.not. associated(OldFace%next_)) then
      self%tail_ => OldFace%prev_
    else
      OldFace%next_%prev_ => OldFace%prev_
    endif

    if (associated(OldFace%prop_)) then
      deallocate(OldFace%prop_)
    endif
    nullify(OldFace%next_)
    nullify(OldFace%prev_)
    deallocate(OldFace)

    OldFace => null()
    self%n_= self%n_- 1
  end subroutine FaceList_remove

  ! Finalizes the components of the given list.
  subroutine FaceList_clean(self)
    class(FaceList),         intent(inout) :: self
    type(face),     pointer                :: OldFace
    if (.not. associated(self%head_)) return
    OldFace => self%head_
    do while (self%n_/= 0)
      call FaceList_remove(self, OldFace)
      OldFace => self%head_
    enddo
  end subroutine FaceList_clean

  ! append a face to the connectivity list
  subroutine FaceConnect_append(self, NewFace, indx)
    class(FaceConnect),    intent(inout) :: self
    type(face), pointer,   intent(inout) :: NewFace
    integer,               intent(in)    :: indx
    ! remove the original connectivity information
    ! Set the previous node of the new node
    NewFace%prevc_(indx)%p_ => self%tail_
    ! set the next node of the new node
    NewFace%nextc_(indx)%p_ => null()
    if (associated(self%head_)) then
      ! connect tail next to the NewFace
      self%tail_%nextc_(indx)%p_ => NewFace
      ! move tail to the current node
      self%tail_ => NewFace
    else
      self%head_ => NewFace
      self%tail_ => self%head_
    endif
    self%n_= self%n_+ 1
  end subroutine FaceConnect_append

  ! append a list of face to the connectivity list
  subroutine FaceConnect_appendList(self, NewFaceList, indx)
    class(FaceConnect), intent(inout) :: self, NewFaceList
    integer,            intent(in)    :: indx

    if (associated(self%head_)) then
      ! set the next of the tail to the head of the new list
      self%tail_%nextc_(indx)%p_ => NewFaceList%head_
      ! set the prev of the head node of the new list to the old tail
      NewFaceList%head_%prevc_(indx)%p_ => self%tail_
    else
      self%head_ => NewFaceList%head_
    endif
    self%tail_ => NewFaceList%tail_
    self%n_= self%n_+ NewFaceList%n_
    NewFaceList%head_ => null()
    NewFaceList%tail_ => null()
    NewFaceList%n_= 0
  end subroutine FaceConnect_appendList

  subroutine FaceConnect_Refine(self, MeshFaceInput, nNew, indx)
    type(FaceConnect),          intent(inout) :: self
    type(face),        pointer, intent(in)    :: MeshFaceInput
    integer,                    intent(in)    :: nNew, indx
    type(face),        pointer                :: MeshFace, MeshFaceOld
    integer                                   :: i
   
    MeshFaceOld => MeshFaceInput
    MeshFace => MeshFaceOld%next_
    do i = 1, nNew - 1
      MeshFace%prevc_(indx)%p_ => MeshFaceOld
      call FaceConnect_insert(self, MeshFace, MeshFace%prevc_(indx)%p_, indx)
      MeshFaceOld => MeshFace
      MeshFace => MeshFace%next_
    enddo
  end subroutine FaceConnect_Refine

  ! Insert element to this list
  subroutine FaceConnect_insert(self, NewFace, PrevFace, indx)
    class(FaceConnect),  intent(inout) :: self
    type(face), pointer, intent(inout) :: NewFace
    integer,             intent(in)    :: indx
    type(face), pointer, intent(inout) :: PrevFace

    if (associated(PrevFace)) then
      ! NewFace next
      NewFace%nextc_(indx)%p_ => PrevFace%nextc_(indx)%p_
      ! prev face next
      PrevFace%nextc_(indx)%p_ => NewFace
    else
      ! NewFace next
      NewFace%nextc_(indx)%p_ => null()
      self%head_              => NewFace
    endif
    ! newface prev link:
    NewFace%prevc_(indx)%p_ => PrevFace
    ! next face's prev
    if (associated(NewFace%nextc_(indx)%p_)) then
      NewFace%nextc_(indx)%p_%prevc_(indx)%p_ => NewFace
    else
      self%tail_ => NewFace
    endif
    self%n_= self%n_+ 1
  end subroutine FaceConnect_insert

  ! Remove Face.
  subroutine FaceConnect_remove(self, OldFace, indx)
    class(FaceConnect),      intent(inout)  :: self
    type(face),     pointer, intent(inout)  :: OldFace
    integer,                 intent(in)     :: indx
    if (.not. associated(OldFace%prevc_(indx)%p_)) then
      self%head_ => OldFace%nextc_(indx)%p_
    else
      OldFace%prevc_(indx)%p_%nextc_(indx)%p_ => OldFace%nextc_(indx)%p_
    endif

    if (.not. associated(OldFace%nextc_(indx)%p_)) then
      self%tail_ => OldFace%prevc_(indx)%p_
    else
      OldFace%nextc_(indx)%p_%prevc_(indx)%p_ => OldFace%prevc_(indx)%p_
    endif
    OldFace%prevc_(indx)%p_ => null()
    OldFace%nextc_(indx)%p_ => null()
    self%n_= self%n_- 1
  end subroutine FaceConnect_remove

  ! Finalizes the components of the given list.
  subroutine FaceConnect_clean(self, indx)
    class(FaceConnect), intent(inout) :: self
    integer, intent(in)                :: indx
    type(face),pointer                 :: OldFace
    if (.not. associated(self%head_)) return
    OldFace => self%head_
    do while (self%n_/= 0)
      call FaceConnect_remove(self, OldFace, indx)
      OldFace => self%head_
    enddo
  end subroutine FaceConnect_clean
end module class_face