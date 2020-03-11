module class_predict
implicit none
private
public :: predict
contains
  subroutine predict(m, l0, nx, ny, nz, ntracer, dt)
    use class_mesh,        only: mesh
    use class_kind,        only: dp
    use class_AMRTypes,    only: column
    use class_parameter,   only: dimX, dimY, Inner, Outer
    type(mesh), intent(inout) :: m(1:)
    integer,    intent(in)    :: l0, nx, ny, nz, ntracer
    real(dp),   intent(in)    :: dt
    type(column), pointer     :: MeshColumn
    integer                   :: idx(m(2)%columns_%n_, dimX:dimY, dimX:dimY)
    real(dp)                  :: xd(m(2)%columns_%n_, dimX:dimY)
    integer                   :: k, i
    do k = 1, nz
      call DepartPosition(m(2)%columns_, m(2)%coord_%x_, dt, xd, k)
      CALL DepartIndex(m(2)%columns_, xd(1:, dimX), m(2)%coord_%x_, m(2)%coord_%dxMax_(dimX), idx(1:, dimX:, dimX), l0, dimX)
      call DepartIndex(m(2)%columns_, xd(1:, dimY), m(2)%coord_%y_, m(2)%coord_%dxMax_(dimY), idx(1:, dimX:, dimY), l0, dimY)

      MeshColumn => m(1)%columns_%head_
      do i = 1, m(1)%columns_%n_
        MeshColumn%prop_%qmid_(:, dimX, Inner) = MeshColumn%prop_%q0_(:, k)
        MeshColumn%prop_%qmid_(:, dimY, Inner) = MeshColumn%prop_%q0_(:, k)
       
        MeshColumn => MeshColumn%next_
      enddo

      call linear(m(2)%columns_, m(1)%columns_, xd(1:, dimX), idx(1:, dimX:, dimX), Inner, dimX, ntracer, nx, ny, l0)
      call linear(m(2)%columns_, m(1)%columns_, xd(1:, dimY), idx(1:, dimX:, dimY), Inner, dimY, ntracer, nx, ny, l0)

      call linear(m(2)%columns_, m(2)%columns_, xd(1:, dimX), idx(1:, dimX:, dimX), Outer, dimX, ntracer, nx, ny, l0)
      call linear(m(2)%columns_, m(2)%columns_, xd(1:, dimY), idx(1:, dimX:, dimY), Outer, dimY, ntracer, nx, ny, l0)

      MeshColumn => m(2)%columns_%head_
      do i = 1, m(2)%columns_%n_
        MeshColumn%prop_%q_(:, k) = 0.5*(MeshColumn%prop_%qmid_(:, dimX, Outer) + MeshColumn%prop_%qmid_(:, dimY, Outer))
        MeshColumn%prop_%qmid_(1:, 1:, 1:) = 0._dp
        MeshColumn => MeshColumn%next_
      enddo
    enddo
  end subroutine predict

  subroutine linear(columns, OldColumns, xd, idx, InOrOut, dim, ntracer, nx, ny, l0)
    use class_AMRTypes,    only: ColumnList, ColumnStack, column
    use class_kind,        only: dp
    use class_parameter,   only: dimX, dimY, Inner, pi
    use class_column,      only: ColumnList_get_column
    use class_column,      only: column_get_leaves
    type(ColumnList),      intent(in)    :: columns, OldColumns
    integer,               intent(in)    :: ntracer, InOrOut, dim, idx(1:, dimX:), nx, ny, l0
    real(dp),              intent(inout) :: xd(1:)
    type(column), pointer                :: MeshColumn, DepartColumn, LeafColumn
    type(ColumnStack)                    :: ColumnLeaves
    real(dp)                             :: dx, xmin, q(ntracer + 1, 2), x0
    integer                              :: i, j, itmp(dimX:dimY), nleaves, Lscale, adjidx
    MeshColumn => columns%head_
    do i = 1, columns%n_
      itmp = idx(i, dimX:)
      dx = 0._dp
      do j = 1, 2
        DepartColumn => ColumnList_get_column(OldColumns, itmp(dimX), itmp(dimY), MeshColumn%l_, l0)
        ColumnLeaves = column_get_leaves(DepartColumn)
        xmin = 2*pi; q(:, j) = 0._dp; nLeaves = 0
        LeafColumn => null()
        do while (associated(ColumnLeaves%head_))
          LeafColumn => ColumnLeaves%head_
          ColumnLeaves%head_ => LeafColumn%nextLeaf_
          LeafColumn%nextLeaf_ => null()
          q(1:, j) = q(1:, j) + LeafColumn%prop_%qmid_(:, 3 - dim, Inner)
          xmin  = min(xmin, LeafColumn%prop_%x_(dim))
          nLeaves = nLeaves + 1
        enddo
        q(:, j) = q(:, j)/nLeaves
        if (MeshColumn%l_ > DepartColumn%l_) then
          Lscale = 2**(MeshColumn%l_ - DepartColumn%l_)
          dx   = dx + 0.5_dp*DepartColumn%prop_%dx_(dim)/Lscale
          xmin = DepartColumn%prop_%x_(dim) + (idx(i, dim) - DepartColumn%prop_%i_(dim)*Lscale)*DepartColumn%prop_%dx_(dim)/Lscale
        else
          Lscale = 2**(DepartColumn%l_ - MeshColumn%l_)
          dx = dx + 0.5_dp*LeafColumn%prop_%dx_(dim)*Lscale
        endif
        if (j == 1) then
          x0 = xmin + dx
          if ((xd(i) - xmin)/dx > 1._dp) then
            adjidx = 1
          else
            adjidx = -1
          endif
          itmp(1) = idx(i, dimX) + (2 - dim)*adjidx
          itmp(2) = idx(i, dimY) + (dim - 1)*adjidx
          itmp = boundary_check(itmp, MeshColumn%l_, nx, ny, l0)
        endif
      enddo
      if (dim == dimX) dx = dx/MeshColumn%prop_%cosp_
      MeshColumn%prop_%qmid_(:, dim, InOrOut) = (q(:, 2) - q(:, 1))*adjidx/dx*(xd(i) - x0) + q(:, 1)
      MeshColumn => MeshColumn%next_
    enddo
  end subroutine linear

  subroutine DepartPosition(columns, x, dt, xd, k)
    use class_AMRTypes,    only: ColumnList, column, face
    use class_kind,        only: dp
    use class_parameter,   only: dimX, dimY, ae, pi
    type(ColumnList),      intent(in)    :: columns
    integer,               intent(in)    :: k
    real(dp),              intent(in)    :: x(0:), dt
    real(dp),              intent(inout) :: xd(1:, dimX:)
    type(column), pointer                :: MeshColumn
    type(face),   pointer                :: MeshFace
    real(dp)                             :: uc(2), ucc
    integer                              :: nx, i, j, indx
    nx = size(x) - 1
    MeshColumn => columns%head_
    do i = 1, columns%n_
      uc(:) = 0._dp
      do indx = 1, 2
        MeshFace => MeshColumn%prop_%faces_(dimX, indx)%head_
        do j = 1, MeshColumn%prop_%faces_(dimX, indx)%n_
          uc(indx) = uc(indx) + MeshFace%prop_%u_(k)
          MeshFace => MeshFace%nextc_(3 - indx)%p_
        enddo
        uc(indx) = uc(indx)/MeshColumn%prop_%faces_(dimX, indx)%n_
      enddo
      ucc = 0.5_dp*sum(uc)

      ! get xd
      xd(i, dimX) = MeshColumn%prop_%x_(dimX) + 0.5_dp*MeshColumn%prop_%dx_(dimX) - ucc*dt/MEshColumn%prop_%cosp_/ae
      xd(i, dimX) = modulo(xd(i, dimX), 2._dp*pi)
      if (xd(i, dimX) > x(nx)) xd(i, dimX) = xd(i, dimX) - 2._dp*pi

      uc(:) = 0._dp
      do indx = 1, 2
        MeshFace => MeshColumn%prop_%faces_(dimY, indx)%head_
        do j = 1, MeshColumn%prop_%faces_(dimY, indx)%n_
          uc(indx) = uc(indx) + MeshFace%prop_%u_(k)
          MeshFace => MeshFace%nextc_(3 - indx)%p_
        enddo
        uc(indx) = uc(indx)/MeshColumn%prop_%faces_(dimY, indx)%n_
      enddo
      ucc = 0.5_dp*sum(uc)
      xd(i, dimY) = MeshColumn%prop_%x_(dimY) + 0.5_dp*MeshColumn%prop_%dx_(dimY) - ucc*dt/ae
      if (xd(i, dimY) > 0.5*pi) xd(i, dimY) = pi - xd(i, dimY)
      if (xd(i, dimY) < -0.5*pi) xd(i, dimY) = -pi - xd(i, dimY)
      MeshColumn => MeshColumn%next_
    enddo
  end subroutine DepartPosition

  subroutine DepartIndex(columns, xd, x, dxMax, idx, l0, dim)
    use class_AMRTypes,    only: ColumnList, column
    use class_kind,        only: dp
    use class_parameter,   only: dimX, dimY, pi
    type(ColumnList),      intent(in)    :: columns
    real(dp),              intent(in)    :: xd(1:), x(0:), dxMax
    integer,               intent(inout) :: idx(1:, dimX:)
    integer,               intent(in)    :: l0, dim
    type(column), pointer                :: MeshColumn
    real(dp)                             :: xd_tmp(columns%n_)
    real(dp)                             :: dx
    integer                              :: i, Lscale, n, i_b, im
    xd_tmp(:) = xd(:)
    MeshColumn => columns%head_
    do i = 1, columns%n_
      ! get index
      idx(i, dim) = INT((xd_tmp(i) - x(0))/dxMax)

      do while (xd_tmp(i) > x(idx(i, dim) + 1) )
        idx(i, dim) = idx(i, dim) + 1
      end do

      if (.not. ((xd_tmp(i) >= x(idx(i, dim))) .and. (xd_tmp(i) <= x(idx(i, dim) + 1)))) then
        print *,"error for index after correction", &
        dim, idx(i, dim), xd(i)*180./pi, x(idx(i, dim))*180./pi, x(idx(i, dim) + 1)*180./pi
      endif
      Lscale = 2**(MeshColumn%l_ - l0)
      dx = (x(idx(i, dim) + 1) - x(idx(i, dim)))/Lscale
      n = 1
      do while (xd_tmp(i) > x(idx(i, dim)) + n*dx)
        n = n + 1
      end do
      idx(i, dim) = idx(i, dim)*Lscale + n - 1
      idx(i, 3 - dim) = MeshColumn%prop_%i_(3 - dim)
      if (dim == dimY) then
        if ((xd(i) > 0.5*pi) .or. (xd(i) < -0.5*pi)) then
          i_b = (size(x) - 1)*Lscale - 1
          im = (i_b + 1)/2          
          idx(i, 3 - dim) = mod(idx(i, 3 - dim) + im, i_b + 1)
        endif
      endif
      MeshColumn => MeshColumn%next_
    enddo
  end subroutine DepartIndex

  function boundary_check(i, l, nx, ny, l0) result(val)
  ! shifting the index based on boundary condition
    use class_parameter,            only: dimX, dimY
    integer, intent(in)                :: i(:)
    integer, intent(in)                :: l, nx, ny
    integer, intent(in)                :: l0
    integer                            :: i_b, j_b, im
    integer                            :: val(2)
    i_b = ishft(nx, l - l0) - 1
    j_b = ishft(ny, l - l0) - 1
    im = ishft(i_b + 1, -1)
    val(:) = i(:)
    if (i(dimY) > j_b) then
      val(dimX) = mod(val(dimX) + im, i_b + 1)
      val(dimY) = j_b - (val(dimY) - j_b) + 1 
    endif

    if (i(dimY) < 0) then
      val(dimX) = mod(val(dimX) + im, i_b + 1)
      val(dimY) = -val(dimY) - 1
    endif

    if (val(dimX) > i_b) val(dimX) = (val(dimX) - i_b - 1)
    if (val(dimX) < 0) val(dimX) = i_b + (val(dimX) + 1)
  end function boundary_check
end module class_predict