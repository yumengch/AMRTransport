module class_depart
implicit none
private
public :: depart
contains
  subroutine Depart(m, Faces, dt, dxMax, nx, x, y, k, l0)
    ! put departure faces into departure columns
    use class_mesh,                    only: mesh
    use class_AMRTypes,                only: FaceList, face
    use class_kind,                    only: dp
    use class_parameter,               only: dimX, dimY, pi, ae, Inner, Outer
    use class_column,                  only: ColumnList_get_column
    type(mesh),     intent(in)    :: m(1:)
    type(FaceList), intent(in)    :: Faces(dimX:)
    real(dp),       intent(in)    :: dt, dxMax(dimX:)
    real(dp),       intent(in)    :: x(0:), y(0:)
    integer,        intent(in)    :: k, nx, l0
    real(dp)                      :: mu, cosp, xd
    type(face), pointer           :: MeshFace
    integer                       :: i, idx(dimX:dimY)

    MeshFace => Faces(dimX)%head_
    do i = 1, Faces(dimX)%n_
      ! get xd
      cosp = (sin(MeshFace%prop_%y_(2)) - sin(MeshFace%prop_%y_(1)))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
      xd = MeshFace%prop_%x_ - MeshFace%prop_%u_(k)*dt/cosp/ae
      xd = modulo(xd, 2._dp*pi)
      if (xd >= x(nx)) xd = xd - 2._dp*pi
      MeshFace%prop_%xd_ = xd
      ! get index
      idx = DepartIndex(MeshFace, x, dxMax(dimX), xd, l0, dimX)
      MeshFace%prop_%DepartColumns_(Inner)%p_ => ColumnList_get_column(m(1)%columns_, idx(dimX), idx(dimY), MeshFace%prop_%l_, l0)
      MeshFace%prop_%DepartColumns_(Outer)%p_ => ColumnList_get_column(m(2)%columns_, idx(dimX), idx(dimY), MeshFace%prop_%l_, l0)
      MeshFace => MeshFace%next_
    enddo


    MeshFace => Faces(dimY)%head_
    do i = 1, Faces(dimY)%n_
      ! get xd
      mu = sin(MeshFace%prop_%x_)
      cosp = cos(MeshFace%prop_%x_)
      xd = mu - MeshFace%prop_%u_(k)*cosp*dt/ae

      if (abs(xd) > 1._dp) then
        print *, "departure position is too far", MeshFace%prop_%i_(dimX), MeshFace%prop_%i_(dimY), xd
      endif
      xd = max(min(xd, 1._dp), -1._dp)
      xd = asin(xd)
      MeshFace%prop_%xd_ = xd
      ! get index
      idx = DepartIndex(MeshFace, y, dxMax(dimY), xd, l0, dimY)
      MeshFace%prop_%DepartColumns_(Inner)%p_ => ColumnList_get_column(m(1)%columns_, idx(dimX), idx(dimY), MeshFace%prop_%l_, l0)
      MeshFace%prop_%DepartColumns_(Outer)%p_ => ColumnList_get_column(m(2)%columns_, idx(dimX), idx(dimY), MeshFace%prop_%l_, l0)
      MeshFace => MeshFace%next_
    enddo
  end subroutine Depart

  function DepartIndex(MeshFace, x, dxMax, xd, l0, dim) result(val)
    ! get the index of the departure face based on departure point
    use class_AMRTypes,                only: face
    use class_kind,                    only: dp
    use class_parameter,               only: dimX, dimY
    type(face), pointer, intent(inout) :: MeshFace
    real(dp),            intent(in)    :: x(0:), dxMax, xd
    integer,             intent(in)    :: l0, dim
    integer                            :: Lscale, i
    real(dp)                           :: dx
    integer                            :: val(dimX:dimY)

    ! get the lowest possible index
    val(dim) = INT((xd - x(0))/dxMax)
    ! find the correct index on the coarsest mesh
    do while (xd > x(val(dim) + 1) )
      val(dim) = val(dim) + 1
    end do

    if (.not. ((xd >= x(val(dim))) .and. (xd <= x(val(dim) + 1)))) then
      print *,"error for index after correction"
    endif

    ! find the index on the refined mesh
    Lscale = 2**(MeshFace%prop_%l_ - l0)
    dx = (x(val(dim) + 1) - x(val(dim)))/Lscale
    MeshFace%prop_%Departdx_(dim) = dx
    i = 1
    do while (xd > x(val(dim)) + i*dx)
      i = i + 1
    end do
    MeshFace%prop_%Departx0_(dim) = x(val(dim)) + (i - 1)*dx
    val(dim) = val(dim)*Lscale + i - 1
    val(3 - dim) = MeshFace%prop_%i_(3 - dim)
  end function DepartIndex
end module class_depart