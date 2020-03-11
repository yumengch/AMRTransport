module class_coord
use class_parameter, only: dimX, dimY, dimZ
use class_kind,      only: dp
implicit none
private
public :: coord, coord_clean
type coord
  real(dp),        allocatable  :: x_(:), y_(:), z_(:)
  real(dp)                      :: dxMax_(dimX:dimZ)
  integer                       :: l0_
end type coord

interface coord
  procedure :: coordinate_constructor
end interface coord

contains
  function coordinate_constructor(nx, ny, nz) result(val)
    use class_parameter,          only: pi
    integer,          intent(in)     :: nx, ny, nz
    real(dp)                         :: dx
    real(dp)                         :: y(ny), zgmu(ny), zgw(ny)

    integer                          :: i, j

    type(coord)                      :: val

    ! initialize the base coordinate
    allocate(val%x_(0:nx))
    allocate(val%y_(0:ny))
    allocate(val%z_(1:nz + 1))
    ! initilize coordinate in x direction
    dx = 2._dp*pi/real(nx, dp)
    val%x_(0) = -0.5_dp*dx
    do i = 1, nx
      val%x_(i) = val%x_(i - 1) + dx
    enddo
    val%dxMax_(dimX) = dx

    ! initalize coordinate in y direction
    call gauaw(zgmu, zgw, ny)
    val%y_(0) = -0.5_dp*pi

    do j = 1, ny
      y(ny + 1 - j) = asin(zgmu(j))
    enddo

    val%y_(1:ny-1) = 0.5_dp*(y(1:ny - 1) + y(2:ny))
    val%y_(ny) = 0.5_dp*pi

    val%dxMax_(dimY) = 0._dp
    do j = 1, ny
      val%dxMax_(dimY) = max(val%y_(j) - val%y_(j - 1), val%dxMax_(dimY))
    enddo

    ! initialize coordinate in z direction
    do i = 1, nz
      val%z_(i) = i
    enddo
    val%l0_     = ceiling(log(real(max(nx, ny)))/log(2.))
  end function coordinate_constructor

  subroutine gauaw(pl,pw, kn)
    use class_parameter, only: pi
    real(dp), intent(out) :: pl(:) !< abscissas of gauss integration
    real(dp), intent(out) :: pw(:) !< weights of the gaussian integration
    integer,  intent(in)  :: kn    !< number of gauss abscissas

    real(dp) :: zfn(0:kn,0:kn), zfnlat(0:kn)
    integer :: iter(kn)

    real(dp) :: z, zfnn
    integer :: jgl, iodd, ik, jn, ins2, isym

    !
    ! 1.0 initialize fourier coefficients for ordinary legendre polynomials
    !
    ! belousov, swarztrauber, and echam use zfn(0,0) = sqrt(2)
    ! ifs normalisation chosen to be 0.5*integral(pnm**2) = 1 (zfn(0,0) = 2.0)
    !

    zfn(0,0) = sqrt(2.0_dp)
    do jn = 1, kn
      zfnn = zfn(0,0)
      do jgl = 1, jn
        zfnn = zfnn*sqrt(1.0_dp-0.25_dp/real(jgl**2,dp))
      enddo

      zfn(jn,jn) = zfnn

      iodd = mod(jn, 2)
      do jgl = 2, jn-iodd, 2
        zfn(jn,jn-jgl) = zfn(jn,jn-jgl+2) &
             *real((jgl-1)*(2*jn-jgl+2),dp)/real(jgl*(2*jn-jgl+1),dp)
      enddo
    enddo

    !
    ! 2.0 gaussian latitudes and weights
    !

    iodd = mod(kn, 2)
    ik = iodd
    do jgl = iodd,kn, 2
      zfnlat(ik) = zfn(kn,jgl)
      ik = ik+1
    enddo

    ins2 = kn/2+mod(kn,2)

    !
    ! 2.1 find first approximation of the roots of the
    !     legendre polynomial of degree kn.
    !

    do jgl = 1, ins2
      z = real(4*jgl-1,dp)*pi/real(4*kn+2,dp)
      pl(jgl) = z+1.0_dp/(tan(z)*real(8*kn**2,dp))
    enddo

    ! 2.2 computes roots and weights for transformed theta

    do jgl = ins2, 1, -1
      call gawl(zfnlat, pl(jgl), pw(jgl), kn, iter(jgl))
    enddo

    ! convert to physical latitude

    do jgl = 1, ins2
      pl(jgl)  = cos(pl(jgl))
    enddo

    do jgl = 1, kn/2
      isym = kn-jgl+1
      pl(isym) = -pl(jgl)
      pw(isym) = pw(jgl)
    enddo

    contains

      subroutine gawl(pfn, pl, pw, kn, kiter)

        integer,  intent(in)    :: kn
        real(dp), intent(in)    :: pfn(0:kn)

        real(dp), intent(out)   :: pw
        integer,  intent(out)   :: kiter

        real(dp), intent(inout) :: pl

        real(dp) :: pmod

        integer :: iflag, itemax, jter, iodd
        real(dp) :: zw
        real(dp) :: zdlx,zdlxn

        ! 1.0 initialization.
        zw = 0
        iflag  =  0
        itemax = 20

        iodd   = mod(kn, 2)

        zdlx   = pl

        ! 2.0 newton iteration

        do jter = 1, itemax+1
          kiter = jter
          call cpledn(kn, iodd, pfn, zdlx, iflag, zw, zdlxn, pmod)
          zdlx = zdlxn
          if (iflag == 1) exit
          if (abs(pmod) <= epsilon(zw)*1000.0_dp) iflag = 1
        enddo

        pl  = zdlxn
        pw  = zw

      end subroutine gawl

      subroutine cpledn(kn, kodd, pfn, pdx, kflag, pw, pdxn, pxmod)

        integer,  intent(in) :: kn
        integer,  intent(in) :: kodd
        real(dp), intent(in) :: pfn(0:kn)
        real(dp), intent(in) :: pdx
        integer,  intent(in) :: kflag

        real(dp), intent(out) :: pw
        real(dp), intent(out) :: pxmod

        real(dp), intent(inout) :: pdxn

        real(dp) :: zdlx, zdlk, zdlldn, zdlxn, zdlmod

        integer :: ik, jn

        ! 1.0 newton iteration step

        zdlx = pdx
        zdlk = 0.0_dp
        if (kodd == 0) zdlk = 0.5_dp*pfn(0)
        zdlxn  = 0.0_dp
        zdlldn = 0.0_dp

        ik = 1

        if (kflag == 0) then
          do jn = 2-kodd, kn, 2
            ! normalised ordinary legendre polynomial == \overbar{p_n}^0
            zdlk = zdlk + pfn(ik)*cos(real(jn,dp)*zdlx)
            ! normalised derivative == d/d\theta(\overbar{p_n}^0)
            zdlldn = zdlldn - pfn(ik)*real(jn,dp)*sin(real(jn,dp)*zdlx)
            ik = ik+1
          enddo
          ! newton method
          zdlmod = -zdlk/zdlldn
          zdlxn  = zdlx+zdlmod
          pdxn   = zdlxn
          pxmod  = zdlmod
        endif

        ! 2.0 computes weight

        if (kflag == 1) then
          do jn = 2-kodd, kn, 2
            ! normalised derivative
            zdlldn = zdlldn - pfn(ik)*real(jn,dp)*sin(real(jn,dp)*zdlx)
            ik = ik+1
          enddo
          pw = real(2*kn+1,dp)/zdlldn**2
        endif

      end subroutine cpledn

  end subroutine gauaw

  subroutine coord_clean(self)
    type(coord), intent(inout) :: self
    deallocate(self%x_, self%y_, self%z_)
  end subroutine coord_clean
end module class_coord