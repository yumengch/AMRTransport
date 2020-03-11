module class_setupWind
implicit none
private
public :: setup_interp, setup_ideal, setup_GetInputUV
contains
  subroutine setup_interp(MeshFace, x, y, l0, nx, ny, nz, uInput, vInput, dim)
    ! interpolate velocity onto high resolution wind
    use class_AMRTypes,             only: face
    use class_kind,                 only: dp
    use class_parameter,            only: dimX, dimY, pi
    real(dp),            intent(in)    :: x(0:), y(0:)
    real(dp),            intent(in)    :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)
    integer,             intent(in)    :: l0, nx, ny, nz, dim
    type(face),          intent(inout) :: MeshFace
    real(dp)                           :: dx0, dx, dy0, dy
    real(dp)                           :: uL(nz), uR(nz), u0, u1
    integer                            :: i0(dimX:dimY), shift(dimX:dimY), itmp(dimX:dimY)
    integer                            :: k, Lscale
    dx0 = 2._dp*pi/nx
    if (dim == dimX) then
      ! get velocity
      Lscale = 2**(MeshFace%prop_%l_ - l0)
      i0(:) = MeshFace%prop_%i_(:)/Lscale
      ! determine the up/down value should be used for interpolation
      shift(:) = 1
      if (MeshFace%prop_%i_(dimY) - i0(dimY)*Lscale < Lscale/2) shift(dimY) = -1
      ! get the distance between two values
      dx = MeshFace%prop_%x_ - x(i0(dimX))
      dy0 = 0.5_dp*(y(i0(dimY) + 1) - y(i0(dimY))) + &  
        0.5_dp*abs(y( merge(i0(dimY) + shift(dimY), 1, i0(dimY) + shift(dimY) >= 0) ) - &
          y(merge(i0(dimY) + shift(dimY) + 1, ny - 1, i0(dimY) + shift(dimY) + 1 <= ny)))
      dy = 0.5*(MeshFace%prop_%y_(1) + MeshFace%prop_%y_(2)) - 0.5_dp*(y(i0(dimY)) + y(i0(dimY) + 1))
      ! get the velocity filed of intermediate interpolations
      itmp(:) = i0(:) + [0, shift(dimY)]
      itmp(:) = boundary_check(itmp, l0, nx, ny, l0)
      uL = uInput(itmp(dimX) + 1, itmp(dimY) + 1, :)
      if (itmp(dimX) /= i0(dimX)) uL(:) = -uL(:)

      itmp(dimX) = modulo(i0(dimX) + shift(dimX), nx)
      itmp(dimY) = i0(dimY) + shift(dimY)
      itmp(:) = boundary_check(itmp, l0, nx, ny, l0)
      uR = uInput(itmp(dimX) + 1, itmp(dimY) + 1, :)
      if (itmp(dimX) /= modulo(i0(dimX) + shift(dimX), nx)) uR(:) = -uR(:)
      ! interpolation
      do k = 1, nz
        u0 = linear_interpolation(uInput(i0(dimX) + 1, i0(dimY) + 1, k), &
          uL(k), dy0, dy, shift(dimY))
        u1 = linear_interpolation(uInput(modulo(i0(dimX) + shift(dimX), &
          nx) + 1, i0(dimY) + 1, k), uR(k), dy0, dy, shift(dimY))
        MeshFace%prop_%u_(k) = linear_interpolation(u0, u1, dx0, dx, shift(dimX))
      enddo
    else
    ! same for Y
      Lscale = 2**(MeshFace%prop_%l_ - l0)
      i0(:) = MeshFace%prop_%i_(:)/Lscale
      if ( (i0(dimY) == ny) .or. (MeshFace%prop_%i_(dimY) == 0) ) then
        MeshFace%prop_%u_(:) = 0._dp
        return
      endif
      ! determine the up/down value should be used for interpolation
      shift(:) = 1
      if (MeshFace%prop_%i_(dimX) - i0(dimX)*Lscale < Lscale/2) shift(dimX) = -1
      ! get the distance between two values
      dy0 = y(i0(dimY) + 1) - y(i0(dimY))
      dy  = MeshFace%prop_%x_ - y(i0(dimY))
      dx  = 0.5*(MeshFace%prop_%y_(1) + MeshFace%prop_%y_(2)) - x(i0(dimX)) - 0.5_dp*dx0
      do k = 1, nz
        u0 = linear_interpolation(vInput(i0(dimX) + 1, i0(dimY), k), &
          vInput(modulo(i0(dimX) + shift(dimX), nx) + 1, i0(dimY), k), dx0, dx, shift(dimX))
        u1 = linear_interpolation(vInput(i0(dimX) + 1, i0(dimY) + 1, k), &
          vInput(modulo(i0(dimX) + shift(dimX), nx) + 1, i0(dimY) + 1, k), dx0, dx, shift(dimX))
        MeshFace%prop_%u_(k) = linear_interpolation(u0, u1, dy0, dy, shift(dimY))
      enddo
      
    endif
  end subroutine setup_interp

  subroutine setup_GetInputUV(TestName, t, dt, beta, x, y, nx, ny, nz, uInput, vInput)
    ! get coarse resolution velocity
    use class_parameter,          only: dimX, dimY, pi
    use class_kind,               only: dp
    character(len =*), intent(in)    :: TestName
    real(dp),          intent(inout) :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)
    real(dp),          intent(in)    :: x(0:), y(0:), dt, beta
    integer,           intent(in)    :: nx, ny, nz, t
    integer                          :: i, j, k
    do i = 1, nx
      do j = 1, ny
        do k = 1, nz
          if ((TestName == "solid") .or. (TestName == "deform_gaussian") .or. &
                  (TestName == "deform_cosbell") .or. (TestName == "uniform")) then
            uInput(i, j, k) = -(StreamFunction(TestName, beta, modulo(x(i - 1), 2*pi), y(j), dt, t) - &
              StreamFunction(TestName, beta, modulo(x(i - 1), 2*pi), y(j - 1), dt, t))/(y(j) - y(j - 1))  
          else
            uInput(i, j, k) = TestCaseU(TestName, beta, modulo(x(i - 1), 2*pi), 0.5_dp*(y(j) + y(j - 1)), dt, t, dimX)
          endif
        enddo
      enddo
    enddo

    do i = 1, nx
      do j = 0, ny
        do k = 1, nz
          if ((TestName == "solid") .or. (TestName == "deform_gaussian") .or. &
                  (TestName == "deform_cosbell") .or. (TestName == "uniform")) then
            vInput(i, j, k) = (StreamFunction(TestName, beta, modulo(x(i), 2*pi), y(j), dt, t) - &
              StreamFunction(TestName, beta, modulo(x(i - 1), 2*pi), y(j), dt, t))/(x(i) - x(i - 1))
            if ((j /= 0) .or. (j /= ny)) vInput(i, j, 1) = vInput(i, j, 1)/cos(y(j))          
          else
            vInput(i, j, k) = TestCaseU(TestName, beta, modulo(0.5_dp*(x(i - 1) + x(i)), 2*pi), y(j), dt, t, dimY)
          endif
        enddo
      enddo
    enddo
  end subroutine setup_GetInputUV

  subroutine setup_ideal(MeshFace, TestName, nz, t, dt, beta, dim)
    ! get velocity for each face
    use class_AMRTypes,             only: face
    use class_kind,                 only: dp
    use class_parameter,            only: dimX, dimY, pi
    character(len =*),   intent(in)    :: TestName
    real(dp),            intent(in)    :: dt, beta
    integer,             intent(in)    :: t, dim, nz
    type(face),          intent(inout) :: MeshFace
    integer                            :: k
    if (TestName == "Old") return
    do k = 1, nz
      if (dim == dimX) then
        if ((TestName == "solid") .or. (TestName == "deform_gaussian") .or. &
                (TestName == "deform_cosbell") .or. (TestName == "uniform")) then
          MeshFace%prop_%u_(k) = -(StreamFunction(TestName, beta, &
            modulo(MeshFace%prop_%x_, 2*pi), MeshFace%prop_%y_(2), dt, t) - &
            StreamFunction(TestName, beta, modulo(MeshFace%prop_%x_, 2*pi), &
              MeshFace%prop_%y_(1), dt, t))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
        else
          MeshFace%prop_%u_(k) = TestCaseU(TestName, beta, &
            modulo(MeshFace%prop_%x_, 2*pi), 0.5_dp*(MeshFace%prop_%y_(1) + MeshFace%prop_%y_(2)), dt, t, dimX)
        endif
      else
          if ((TestName == "solid") .or. (TestName == "deform_gaussian") .or. &
                  (TestName == "deform_cosbell") .or. (TestName == "uniform")) then
            MeshFace%prop_%u_(k) = (StreamFunction(TestName, beta, modulo(MeshFace%prop_%y_(2), 2*pi), MeshFace%prop_%x_, dt, t) - &
              StreamFunction(TestName, beta, modulo(MeshFace%prop_%y_(1), 2*pi), &
                MeshFace%prop_%x_, dt, t))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
            if (.not. associated(MeshFace%pole_)) MeshFace%prop_%u_(k) = MeshFace%prop_%u_(k)/cos(MeshFace%prop_%x_)          
          else
            MeshFace%prop_%u_(k) = TestCaseU(TestName, beta, &
              modulo(0.5_dp*(MeshFace%prop_%y_(2) + MeshFace%prop_%y_(1)), 2*pi), MeshFace%prop_%x_, dt, t, dimY)
          endif
      endif
    enddo
  end subroutine setup_ideal

  function TestCaseU(TestName, beta, x, y, dt, t, dim)  result(q)
  ! give streamfunction to each cell
    use class_kind,              only: dp
    use class_parameter,         only: dimX, ae, pi
    real(dp), intent(in)            :: x, y, beta, dt
    integer, intent(in)             :: t, dim
    character(len = *), intent(in)  :: TestName
    real(dp)                        :: rho, x_r, y_r, w, x0, y0
    real(dp)                        :: u0, T0, sinxr
    real(dp)                        :: q

    T0 = 12._dp*3600._dp*24._dp; 
    u0 = 2*pi*ae/T0
    x0 = 1.5*pi; y0 = 0._dp
    q = 0._dp
    if (TestName == "moving") then
      x_r = atan2(cos(y0)*sin(x0 - pi), cos(y0)*sin(0.5*pi - beta)*cos(x0 - pi) - cos(0.5*pi - beta)*sin(y0))
      if (x_r < 0._dp) x_r = x_r + 2._dp*pi
      y_r = asin(sin(y0)*sin(0.5*pi - beta) + cos(y0)*cos(0.5*pi - beta)*cos(x0 - pi))

      x0 = atan2(cos(y_r)*sin(x_r + u0*dt*t/ae), &
        sin(y_r)*cos(0.5*pi - beta) + cos(y_r)*cos(x_r + u0*dt*t/ae)*sin(0.5*pi - beta)) + pi
      y0 = asin(sin(y_r)*sin(0.5*pi - beta) - cos(y_r)*cos(0.5*pi - beta)*cos(x_r + u0*dt*t/ae))


      sinxr = (sin(y)*sin(y0) + cos(y)*cos(y0)*cos(x - x0))**2
      if (sinxr > 1.0_dp) sinxr = 1.0_dp
      sinxr = 1._dp - sinxr
      if (sinxr < 0._dp) sinxr = 0._dp
      rho = 3._dp*sqrt(sinxr)

      if (abs(rho) > epsilon(1._dp)) then
        w = (u0*1.5_dp*sqrt(3.)*(1 - (tanh(rho))**2)*tanh(rho))/rho
      else
        w = 0.
      endif

      if (dim == dimX) then
        q = u0*(cos(y)*cos(beta) + sin(y)*cos(x)*sin(beta)) + &
          w*(sin(y0)*cos(y) - cos(y0)*cos(x - x0)*sin(y))
      else
        q = -u0*sin(x)*sin(beta) + &
         w*(cos(y0)*sin(x - x0))
      endif
    else if (TestName == "accelerating") then
      x_r = atan2(cos(y)*sin(x - pi), cos(y)*sin(0.5*pi - beta)*cos(x - pi) - cos(0.5*pi - beta)*sin(y))
      if (x_r < 0._dp) x_r = x_r + 2._dp*pi

      if ((x_r > 0.5*pi) .and. (x_r < 1.5*pi)) then
        if (dim == dimX) then
          q = u0*(cos(y)*cos(beta) + sin(y)*cos(x)*sin(beta))
        else
          q = -u0*sin(x)*sin(beta)
        endif
      else if ((x_r >= 1.5*pi) .and. (x_r < 2*pi)) then
        u0 = u0*(1.5*cos(4.*pi - 2.*x_r) + 2.5)
        if (dim == dimX) then

          q = u0*(cos(y)*cos(beta) + sin(y)*cos(x)*sin(beta))
        else
          q = -u0*sin(x)*sin(beta)
        endif
      else
        u0 = u0*(1.5*cos(2*x_r) + 2.5)
        if (dim == dimX) then
          q = u0*(cos(y)*cos(beta) + sin(y)*cos(x)*sin(beta))
        else
          q = -u0*sin(x)*sin(beta)
        endif
      endif
    else if (TestName == "divergent") then
      u0 = 5._dp*ae/T0
      if (dim == dimX) then
        q = -u0*((sin((x - 2*pi*t*dt/T0)*0.5)*cos(y))**2)*sin(2*y)*cos(pi*dt*t/T0) + 2*pi*ae*cos(y)/T0
      else
        q = 0.5*u0*sin(x - 2*pi*t*dt/T0)*((cos(y))**3)*cos(pi*dt*t/T0)
      endif
    endif
  end function TestCaseU

  function StreamFunction(TestName, beta, x, y, dt, t) result(q)
  ! give streamfunction to each cell
    use class_kind,              only: dp
    use class_parameter,         only: ae, pi
    real(dp), intent(in)            :: beta, x, y, dt
    integer, intent(in)             :: t
    character(len = *), intent(in)  :: TestName
    real(dp)                        :: u0, T0
    real(dp)                        :: q

    T0 = 12._dp*3600._dp*24._dp; 
    if (TestName == "solid") then
      u0 = 2.0*pi*ae/T0
      q =  -u0*(sin(y)*cos(beta) - cos(x)*cos(y)*sin(beta))
    else if ((TestName == "deform_gaussian") .or. (TestName == 'deform_cosbell')) then
      u0 = 10._dp*ae/T0
      q = u0*((sin(x - 2*pi*t*dt/T0)*cos(y))**2)*cos(pi*t*dt/T0) - (2*pi*ae/T0)*sin(y)
    else if (TestName == "uniform") then
      u0 = 2.0*pi*ae/T0
      q =  -u0*(sin(y)*cos(beta) - cos(x)*cos(y)*sin(beta))
    else
      print *, "Non-existence Initial Profile"
      stop
    endif
  end function StreamFunction

  function linear_interpolation(u0, u1, dx0, dx, ishift) result(u)
    use class_parameter,               only: dp
    real(dp), intent(in)          :: u0, u1
    real(dp), intent(in)          :: dx0, dx
    integer,  intent(in)          :: ishift
    real(dp)                      :: u
    u = (u1 - u0)/dx0*dx*ishift + u0
  end function linear_interpolation

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
end module class_setupWind