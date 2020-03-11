module class_parameter
use class_kind, only: dp
implicit none

real(dp), parameter  :: ae = 6.371229e6_dp, pi = 3.14159265358979323846264338327950288419717_dp
integer,  parameter  :: dimX = 1, dimY = 2, dimZ = 3 ! parameters for the dimensions
integer,  parameter  :: inner = 1, outer = 2
 
integer,  parameter  :: ntracer = 0
integer,  parameter  :: nz = 1
integer              :: nx, ny
integer              :: MaxRef, nLimit, nof, nt
real(dp)             :: theta_c, theta_r
character(len = 256) :: TestName
real(dp)             :: beta, dt
namelist /grid/ nx, ny, nt, MaxRef, nLimit, theta_c, theta_r
namelist /output/ nof
namelist /test/ TestName, beta
contains
  subroutine read_var()
    open(unit = 10, file = "test.nml")
    read(10, nml=grid)
    read(10, nml = output)
    read(10, nml = test)
    close(10)
    beta = beta*pi
    dt   = 12._dp*3600._dp*24._dp/nt
  end subroutine read_var
end module class_parameter