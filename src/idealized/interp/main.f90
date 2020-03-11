program main
use class_kind,      only: dp
use class_mesh,      only: mesh, mesh_initOther, mesh_sync
use class_setup,     only: setup_init, norms
use class_refine,    only: mesh_refineInit, mesh_refine
use class_FFSL,      only: FFSL, update
use class_parameter, only: read_var, nx, ny, nz, ntracer, MaxRef, nLimit, theta_r, theta_c, TestName, beta, dt, nt, nof
use class_netcdf,    only: netcdf_init, netcdf_writevar
use class_setupwind, only: setup_GetInputUV
implicit none
type(mesh) :: m(2)
integer    :: it
real(dp), allocatable   :: uInput(:, :, :), vInput(:, :, :)
call read_var()
allocate(uInput(nx, ny, nz), vInput(nx, 0:ny, nz))
uInput(:, :, :) = 0._dp
vInput(:, 0:, :) = 0._dp
m(1) = mesh(nx, ny, nz, uInput, vInput)
m(2) = mesh(nx, ny, nz, uInput, vInput)
call mesh_initOther(m(2), m(1))
call setup_init(m(1)%columns_, nz, trim(adjustl(TestName)))
call setup_init(m(2)%columns_, nz, trim(adjustl(TestName)))
call mesh_refineInit(m, nx, ny, nz, ntracer, MaxRef, theta_r, theta_c, dt, uInput, vInput, nLimit)
call mesh_sync(m(1), nx, ny, nz, uInput, vInput)
call update(m)
call norms(m(2)%columns_, trim(adjustl(TestName)), 0, beta, dt)
call netcdf_init(nx, ny, nz, ntracer, MaxRef)
call netcdf_writevar(m(2)%columns_, 0, dt, nz, ntracer)
print *, 'nColumns', 0, m(2)%columns_%n_
do it = 1, nt
  call setup_GetInputUV(trim(adjustl(TestName)), it, dt, beta, m(1)%coord_%x_, m(1)%coord_%y_, nx, ny, nz, uInput, vInput)
  call mesh_refine(m, nx, ny, nz, ntracer, MaxRef, theta_r, theta_c, dt, uInput, vInput, nLimit)
  call FFSL(m, nx, ny, nz, ntracer, MaxRef, dt)
  call mesh_sync(m(1), nx, ny, nz, uInput, vInput)
  call update(m)
  call norms(m(1)%columns_, trim(adjustl(TestName)), it, beta, dt)
  print *, 'nColumns', it, m(2)%columns_%n_
  if (mod(it, nt/nof) == 0) &
  call netcdf_writevar(m(2)%columns_, it, dt, nz, ntracer)
  ! if (it == 1) stop
enddo
end program main
