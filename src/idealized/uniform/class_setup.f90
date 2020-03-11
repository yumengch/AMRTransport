module class_Setup
implicit none
private
public :: setup_init, norms
contains
  subroutine setup_init(columns, nz, TestName)
    ! setup initial value
    use class_AMRTypes,       only: ColumnList, column
    use class_kind, only: dp
    character(len = *), intent(in)  :: TestName
    type(ColumnList), intent(in) :: columns
    integer,          intent(in) :: nz
    type(column), pointer        :: MeshColumn
    integer                      :: i, k, itracer
    ! concentrations
    MeshColumn => columns%head_
    do i = 1, columns%n_
      do k = 1, nz
        itracer = 1
        MeshColumn%prop_%q0_(itracer, k) = TestCaseQ(TestName, MeshColumn%prop_%x_(:) + 0.5_dp*MeshColumn%prop_%dx_(:))
        MeshColumn%prop_%q_(itracer, k) = MeshColumn%prop_%q0_(itracer, k)
      enddo
      MeshColumn => MeshColumn%next_
    end do
  end subroutine setup_init

  function TestCaseQ(TestName, x) result(q)
  ! initial tracer distribution
    use class_parameter, only: dimX, dimY, ae, pi
    use class_kind,      only: dp
    real(dp), intent(in)            :: x(1:)
    character(len = *), intent(in)  :: TestName
    real(dp)                        :: x0, y0, x1,  y1
    real(dp)                        :: q0, r, r0, coeff
    real(dp)                        :: q
    real(dp)                        :: a, b, c, a0, b0, c0, a1, b1, c1
    real(dp)                        :: y_r, rho, x_r
    x0 = 1.5*pi; y0 = 0._dp; x1 = 0.; y1 = 0.;
    if (TestName == "solid") then
      q0 = 1.
      r = haversine(y0, x0, x(dimY), x(dimX))
      r0 = ae/3.
      if (r < r0) then
        q = 0.5*q0*(1+cos(pi*r/(r0)))
      else
        q = 0.0_dp
      endif

    else if (TestName == "deform_gaussian") then

      x0 = 5*pi/6; y0 = 0; x1 = 7*pi/6; y1 = 0
      ! x0 = 3*pi/4; y0 = 0; x1 = 5*pi/4; y1 = 0
      q0 = 0.95_dp; coeff = 5._dp
      a = cos(x(dimX))*cos(x(dimY)); b = sin(x(dimX))*cos(x(dimY)); c = sin(x(dimY))
      a0 = cos(x0)*cos(y0); b0 = sin(x0)*cos(y0); c0 = sin(y0)
      a1 = cos(x1)*cos(y1); b1 = sin(x1)*cos(y1); c1 = sin(y1)
      q = q0*exp(- coeff*((a - a0)**2 + (b - b0)**2 + (c - c0)**2)) + &
        q0*exp(- coeff*((a - a1)**2 + (b - b1)**2 + (c - c1)**2))
    else if ((TestName == "deform_cosbell") .or. (TestName == "divergent")) then

      q = 0.1_dp
      x0 = 5*pi/6; y0 = 0; x1 = 7*pi/6; y1 = 0
      q0 = 0.9_dp; coeff = 5._dp
      r = haversine(y0, x0, x(dimY), x(dimX))
      if (r < 0.5*ae) then
        q = q + q0*0.5*(1+cos(pi*r/(0.5*ae)))
      endif
      r = haversine(y1, x1, x(dimY), x(dimX))
      if (r < 0.5*ae) then
        q = q + q0*0.5*(1+cos(pi*r/(0.5*ae)))
      endif

    else if (TestName == "moving") then
      x_r = atan2(cos(x(dimY))*sin(x(dimX) - x0), cos(x(dimY))*sin(y0)*cos(x(dimX) - x0) - cos(y0)*sin(x(dimY)))
      if (x_r < 0._dp) x_r = x_r + 2._dp*pi
      y_r = asin(sin(x(dimY))*sin(y0) + cos(x(dimY))*cos(y0)*cos(x(dimX) - x0))

      rho = 3._dp*cos(y_r)
      q = 1 - tanh((rho/5.)*sin(x_r))
      
    else if (TestName == "uniform") then

      q = 1.0_dp
    else if (TestName == "accelerating") then
      q0 = 1.
      r = haversine(y0, x0, x(dimY), x(dimX))
      if (r < ae/3.) then
        q = 0.5*q0*(1+cos(pi*r/(ae/3.)))
      else
        q = 0.0_dp
      endif
    else
      print *, "Non-existence Initial Profile"
      stop
    endif
  end function TestCaseQ

  function analy_sol(initial, x, y, t, beta, dt)  result(q)
  ! give streamfunction to each cell
    use class_kind,               only: dp
    use class_parameter,          only: ae, pi
    real(dp), intent(in)            :: x, y, beta, dt
    integer, intent(in)             :: t
    character(len = *), intent(in):: initial
    real(dp)                        :: x0, x1, y0, y1
    real(dp)                        :: x_c, y_c
    real(dp)                        :: x_d, y_d
    real(dp)                        :: rho, x_r, y_r, w
    real(dp)                        :: a, b, c, a0, b0, c0, a1, b1, c1
    real(dp)                        :: u0, T0, q0, r, coeff, r0
    real(dp)                        :: q
    q = 0._dp
    T0 = 12._dp*3600._dp*24._dp; 
    u0 = 2*pi*ae/T0
    x0 = 1.5*pi; y0 = 0._dp

    ! get the position of the center of vortex/cosbell
    x_r = atan2(cos(y0)*sin(x0 - pi), cos(y0)*sin(0.5*pi - beta)*cos(x0 - pi) - cos(0.5*pi - beta)*sin(y0))
    if (x_r < 0._dp) x_r = x_r + 2._dp*pi
    y_r = asin(sin(y0)*sin(0.5*pi - beta) + cos(y0)*cos(0.5*pi - beta)*cos(x0 - pi))

    x_c = atan2(cos(y_r)*sin(x_r + u0*dt*t/ae), &
      sin(y_r)*cos(0.5*pi - beta) + cos(y_r)*cos(x_r + u0*dt*t/ae)*sin(0.5*pi - beta)) + pi
    y_c = asin(sin(y_r)*sin(0.5*pi - beta) - cos(y_r)*cos(0.5*pi - beta)*cos(x_r + u0*dt*t/ae))

    if (initial == "solid") then
      q0 = 1.
      r = haversine(y_c, x_c, y, x)
      r0 = ae/3.
      if (r < r0) then
        q = 0.5*q0*(1+cos(pi*r/r0))
      else
        q = 0.0_dp
      endif
    else if (initial == "moving") then
      ! get angular velocity controlling the vortex
      rho = 3._dp*sqrt(1 - (sin(y)*sin(y_c) + cos(y)*cos(y_c)*cos(x - x_c))**2)
      if (abs(rho) >= epsilon(1._dp)) then
        w = (u0*1.5_dp*sqrt(3.)*(1 - (tanh(rho))**2)*tanh(rho))/rho
      else
        w = 0.
      endif 

     ! get the departure position of solid body rotation
      x_r = atan2(cos(y)*sin(x - pi), cos(y)*sin(0.5*pi - beta)*cos(x - pi) - cos(0.5*pi - beta)*sin(y))
      if (x_r < 0._dp) x_r = x_r + 2._dp*pi
      y_r = asin(sin(y)*sin(0.5*pi - beta) + cos(y)*cos(0.5*pi - beta)*cos(x - pi))


      x_d = atan2(cos(y_r)*sin(x_r - u0*dt*t/ae), &
        sin(y_r)*cos(0.5*pi - beta) + cos(y_r)*cos(x_r - u0*dt*t/ae)*sin(0.5*pi - beta)) + pi
      y_d = asin(sin(y_r)*sin(0.5*pi - beta) - cos(y_r)*cos(0.5*pi - beta)*cos(x_r - u0*dt*t/ae))

      ! get rho and x_r for tracer distribution without vortex
      rho = 3._dp*sqrt(1 - (sin(y_d)*sin(y0) + cos(y_d)*cos(y0)*cos(x_d - x0))**2)
      x_r = atan2(cos(y_d)*sin(x_d - x0), cos(y_d)*sin(y0)*cos(x_d - x0) - cos(y0)*sin(y_d))
      if (x_r < 0._dp) x_r = x_r + 2._dp*pi
      
      q = 1 - tanh((rho/5.) *sin(x_r - w*t*dt/ae))

    else if (initial == "deform_gaussian") then

      x0 = 5*pi/6; y0 = 0; x1 = 7*pi/6; y1 = 0
      ! x0 = 3*pi/4; y0 = 0; x1 = 5*pi/4; y1 = 0
      q0 = 0.95_dp; coeff = 5._dp
      a = cos(x)*cos(y); b = sin(x)*cos(y); c = sin(y)
      a0 = cos(x0)*cos(y0); b0 = sin(x0)*cos(y0); c0 = sin(y0)
      a1 = cos(x1)*cos(y1); b1 = sin(x1)*cos(y1); c1 = sin(y1)
      q = q0*exp(- coeff*((a - a0)**2 + (b - b0)**2 + (c - c0)**2)) + &
        q0*exp(- coeff*((a - a1)**2 + (b - b1)**2 + (c - c1)**2))
      ! q = 1._dp
    else if ((initial == "deform_cosbell") .or. (initial == "divergent")) then

      q = 0.1_dp
      x0 = 5*pi/6; y0 = 0; x1 = 7*pi/6; y1 = 0
      q0 = 0.9_dp; coeff = 5._dp
      r = haversine(y0, x0, y, x)
      if (r < 0.5*ae) then
        q = q + q0*0.5*(1+cos(pi*r/(0.5*ae)))
      endif
      r = haversine(y1, x1, y, x)
      if (r < 0.5*ae) then
        q = q + q0*0.5*(1+cos(pi*r/(0.5*ae)))
      endif

    else if ( initial == "accelerating") then
      q0 = 1.
      r = haversine(y0, x0, y, x)
      if (r < ae/3.) then
        q = 0.5*q0*(1+cos(pi*r/(ae/3.)))
      else
        q = 0.0_dp
      endif

    else if (initial == "uniform") then

      q = 1.0_dp

    endif
  end function analy_sol

  subroutine norms(columns, TestName, t, beta, dt)
    use class_AMRTypes,        only: ColumnList, column
    use class_column,          only: column_dA
    use class_kind,            only: dp
    use class_parameter,       only: dimX, dimY, pi
    type(ColumnList),   intent(in) :: columns
    integer,            intent(in) :: t
    real(dp),           intent(in) :: beta, dt
    character(len = *), intent(in) :: TestName                               ! current time steps
    real(dp)                       :: l1, l2, linf                           ! error norms
    real(dp)                       :: l1_div, l2_div, linf_div, error, exact ! divisions of error norms and exact solutions
    real(dp)                       :: dA, qmax, qmin, mass
    type(column), pointer          :: MeshColumn
    integer                        :: i

    l1 = 0.; l1_div = 0.; l2 = 0.; l2_div = 0.
    linf = 0.; linf_div = 0.; exact = 0._dp; dA = 0._dp
    qmax = -999999; qmin = 999999;
    ! concentrations
    mass = 0._dp
    MeshColumn => columns%head_
    do i = 1, columns%n_
        exact = analy_sol(TestName, modulo(MeshColumn%prop_%x_(dimX) + &
          0.5_dp*MeshColumn%prop_%dx_(dimX), 2*pi), MeshColumn%prop_%x_(dimY) &
          + 0.5_dp*MeshColumn%prop_%dx_(dimY), t, beta, dt)
        error = abs(MeshColumn%prop_%q_(1, 1) - exact)
        dA    = column_dA(MeshColumn)
        qmax  = max(MeshColumn%prop_%q_(1, 1), qmax)
        qmin  = min(MeshColumn%prop_%q_(1, 1), qmin)
        ! L1 norm
        l1 = l1 + error*dA
        l1_div = l1_div + abs(exact)*dA
        ! L2 norm
        l2 = l2 + error**2*dA
        l2_div = l2_div + exact**2*dA
        ! l_inf norm
        linf = max(linf, error)
        linf_div = max(linf_div, abs(exact))
        mass = mass + MeshColumn%prop_%q_(1, 1)*dA
        MeshColumn => MeshColumn%next_
    end do
    ! normalize error norms
    l1 = l1/l1_div
    l2 = sqrt(l2/l2_div)
    linf = linf/linf_div
    print *, "at time", t, "l1 : ", l1, "l2 : ", l2, "linf : ", linf, "mass", mass
    print *, 'max:', qmax, "min", qmin
  end subroutine norms

  function haversine(deglat1,deglon1,deglat2,deglon2) result (dist)
    ! great circle distance -- adapted from Matlab
    use class_parameter, only: ae
    use class_kind, only: dp
    REAL(dp),intent(in) :: deglat1,deglon1,deglat2,deglon2
    REAL(dp) :: a,dist,dlat,dlon,lat1,lat2

    dlat = deglat2-deglat1
    dlon = deglon2-deglon1
    lat1 = deglat1
    lat2 = deglat2
    ! a = (sin(dlat/2))**2 + cos(lat1)*cos(lat2)*(sin(dlon/2))**2
    ! c = 2._dp*asin(sqrt(a))
    a = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon)
    dist = ae*acos(a)
  end function haversine
end module class_Setup
