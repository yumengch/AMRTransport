program main
use class_mesh,      only: mesh, mesh_initOther, mesh_sync
use class_setup,     only: setup_init, norms
use class_refine,    only: mesh_refineInit, mesh_refine
use class_FFSL,      only: FFSL, update
use class_parameter, only: read_var, nx, ny, nz, ntracer, MaxRef, nLimit, theta_r, theta_c, TestName, beta, dt, nt, nof
use class_netcdf,    only: netcdf_init, netcdf_writevar

implicit none
type(mesh) :: m(2)
integer    :: it
real(dp)    :: start_scheme, end_scheme
real(dp)    :: start_refine, end_refine
real(dp)    :: scheme_time, refine_time
call read_var()
m(1) = mesh(nx, ny, nz, dt, beta, trim(adjustl(TestName)))
m(2) = mesh(nx, ny, nz, dt, beta, trim(adjustl(TestName)))
call mesh_initOther(m(2), m(1))

call setup_init(m(2)%columns_, nz, trim(adjustl(TestName)))
call mesh_refineInit(m, nx, ny, nz, 0, ntracer, MaxRef, theta_r, theta_c, dt, beta, trim(adjustl(TestName)), nLimit)
call mesh_sync(m(1), nz, 0, dt, beta)

call setup_init(m(2)%columns_, nz, trim(adjustl(TestName)))
call setup_init(m(1)%columns_, nz, trim(adjustl(TestName)))
call norms(m(2)%columns_, trim(adjustl(TestName)), 0, beta, dt)
call netcdf_init(nx, ny, nz, ntracer, MaxRef)
call netcdf_writevar(m(2)%columns_, 0, dt, nz, ntracer)
print *, 'nColumns', 0, m(2)%columns_%n_
do it = 1, nt
  refine_time = 0.
  scheme_time = 0.
  call cpu_time(start_refine)
  call mesh_refine(m, nx, ny, nz, it, ntracer, MaxRef, theta_r, theta_c, dt, beta, trim(adjustl(TestName)), nLimit)
  call cpu_time(end_refine)
  refine_time = refine_time + (start_refine - end_refine)

  call cpu_time(start_scheme)
  call FFSL(m, nx, ny, nz, ntracer, MaxRef, dt)
  call cpu_time(end_scheme)
  scheme_time = scheme_time + (start_scheme - end_scheme)

  call cpu_time(start_refine)
  call mesh_sync(m(1), nz, it, dt, beta)
  call cpu_time(end_refine)
  refine_time = refine_time + (start_refine - end_refine)

  call cpu_time(start_scheme)
  call update(m)
  call cpu_time(end_scheme)
  scheme_time = scheme_time + (start_scheme - end_scheme)

  print *, 'nColumns', it, m(2)%columns_%n_, refine_time, scheme_time
  call norms(m(1)%columns_, trim(adjustl(TestName)), it, beta, dt)
  if (mod(it, nt/nof) == 0) &
  call netcdf_writevar(m(2)%columns_, it, dt, nz, ntracer)
enddo
end program main
