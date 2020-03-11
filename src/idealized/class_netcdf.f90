module class_netcdf
  ! use class_kind,              only: dp
  implicit none
  private
  public :: netcdf_init, netcdf_writevar
  type strings
    character(len = 256) :: c
  end type strings
  integer,  parameter   :: dp = selected_real_kind (15, 307)!selected_real_kind(pd,rd)
  integer                       :: ncid, rec_varid
  integer,  allocatable         :: varid(:)
  real(dp), allocatable         :: time(:), timeOld(:)
  real(dp), allocatable         :: lons(:), lats(:), verts(:, :, :), q(:, :, :)
  integer,  allocatable         :: i_(:), j_(:), l_(:)
  integer                       :: nColumnsTotal
  integer                       :: it
  character(len = *), parameter :: filename = 'AMRDUST.nc' 
contains
  subroutine netcdf_init(nx, ny, nz, ntracer, MaxRef)
    use class_parameter,      only: dimX, dimY
    integer, intent(in)          :: nx, ny, nz, ntracer, MaxRef
    integer                      :: Lscale
    character(len = 256)         :: cyear, cmonth, cday, chour, cminute, csecond
    character(len = 256)         :: chourtmp, cminutetmp, csecondtmp
    character(len = 256)         :: Tunit
    integer                      :: year,month,day,hour,minute,second
    type(strings), allocatable   :: names(:), units(:)
    integer                      :: itracer

    if (allocated(names)) deallocate(names)
    allocate(names(8 + ntracer + 1))
    if (allocated(units)) deallocate(units)
    allocate(units(8 + ntracer + 1))
    if (allocated(varid)) deallocate(varid)
    allocate(varid(ntracer + 8 + 1))

    Lscale = 2**MaxRef
    nColumnsTotal = Lscale*nx*Lscale*ny

    it = 0
    year = 2006; month = 10; day =1; hour = 0; minute = 0; second = 0
    write(cyear  , *) year
    write(cmonth , *) month
    write(cday   , *) day
    write(chourtmp  , *) hour
    if (hour < 10) then
      write(chour  , *) '0', trim(adjustl(chourtmp))
    else
      write(chour  , *) trim(adjustl(chourtmp))
    endif
    write(cminutetmp, *) minute
    if (minute < 10) then
      write(cminute, *) '0', trim(adjustl(cminutetmp))
    else
      write(cminute, *) trim(adjustl(cminutetmp))
    endif
    write(csecondtmp, *) second
    if (second < 10) then
      write(csecond, *) '0', trim(adjustl(csecondtmp))
    else
      write(csecond, *) trim(adjustl(csecondtmp))
    endif
    write(Tunit  , *) 'days since ', trim(adjustl(cyear)), '-', &
                      trim(adjustl(cmonth)), '-', trim(adjustl(cday)), &
                      ' ', trim(adjustl(chour)), ':', trim(adjustl(cminute)), &
                      ':', trim(adjustl(csecond))

    print *, 'netcdf_init:', trim(adjustl(Tunit))
    names(1)%c   = 'nCols'
    names(2:4)%c = ['i', 'j', 'l']
    names(5:6)%c = ['lon', 'lat']
    names(7:8)%c = ['vertx', 'verty']
    units(1:8)%c = ['none', 'none', 'none', 'none', 'none','none', 'none', 'none']
    names(9)%c = 'q'
    do itracer = 9, ntracer + 1 + 8
      ! names(itracer)%c = 'q'
      units(itracer)%c = 'kg kg-1'
    enddo

    ! setup output netcdf
    call writenc_setup(nColumnsTotal, nz, 4, ntracer + 8 + 1, trim(adjustl(Tunit)), names, units)
    deallocate(names)
    deallocate(units)
    allocate(i_(nColumnsTotal), j_(nColumnsTotal), l_(nColumnsTotal))
    allocate(lons(nColumnsTotal))
    allocate(lats(nColumnsTotal))
    allocate(verts(4, nColumnsTotal, dimX:dimY))
    allocate(q(ntracer + 1, nz, nColumnsTotal))        
  end subroutine netcdf_init

  subroutine writenc_setup(nColumns, nz, nvert, nvars, Tunit, varname, varunits)
    use netcdf
    implicit none
    ! This is the name of the data file we will create.
    integer,             intent(in) :: nz, nColumns, nvars, nvert
    character (len = *), intent(in) :: Tunit
    type(strings)      , intent(in) :: varname(:), varunits(:)
    integer                         :: lvl_dimid, col_dimid, rec_dimid, vert_dimid
   ! integer                         :: lvl_varid, col_varid, vert_varid
    integer                         :: dimsid1D(1), dimsid2D(2), dimsid3D(3), dimsid3DV(3)
    ! Loop indices
    integer                         :: i
    ! Create the file. 
    print *, 'write_setup: AMR output filename:', filename
    call check( nf90_create(filename, nf90_clobber, ncid) )
    
    call check( nf90_def_dim(ncid, "time",  nf90_unlimited, rec_dimid)  )
    call check( nf90_def_dim(ncid, "lev",   nz,             lvl_dimid)  )
    call check( nf90_def_dim(ncid, "col",   nColumns,       col_dimid)  )
    call check( nf90_def_dim(ncid, "vert",  nvert,          vert_dimid) )

    ! call check( nf90_def_var(ncid, "lev",  nf90_double, lvl_dimid,  lvl_varid) )
    ! call check( nf90_def_var(ncid, "col",  NF90_INT,    col_dimid,  col_varid) )
    ! call check( nf90_def_var(ncid, "vert", NF90_INT,    vert_dimid, vert_varid) )
    call check( nf90_def_var(ncid, "time",  nf90_double, rec_dimid,  rec_varid) )

    ! call check( nf90_put_att(ncid, lvl_varid, "units",    "level") )
    ! call check( nf90_put_att(ncid, lvl_varid, "positive", "down") )
    ! call check( nf90_put_att(ncid, col_varid, "units",    "none") )
    call check( nf90_put_att(ncid, rec_varid, "units", Tunit) )
    call check( nf90_put_att(ncid, rec_varid, 'calendar', "proleptic_gregorian") )
    call check( nf90_put_att(ncid, rec_varid, 'axis', "T") )

    dimsid1D  = [rec_dimid]
    dimsid3DV = [vert_dimid, col_dimid, rec_dimid]
    dimsid2D  = [col_dimid, rec_dimid]
    dimsid3D  = [lvl_dimid, col_dimid, rec_dimid]
    call check( nf90_def_var(ncid, trim(adjustl(varname(1)%c)), NF90_INT,    dimsid1D,  varid(1)) )
    call check( nf90_def_var(ncid, trim(adjustl(varname(2)%c)), NF90_INT,    dimsid2D,  varid(2)) )
    call check( nf90_def_var(ncid, trim(adjustl(varname(3)%c)), NF90_INT,    dimsid2D,  varid(3)) )
    call check( nf90_def_var(ncid, trim(adjustl(varname(4)%c)), NF90_INT,    dimsid2D,  varid(4)) )    
    call check( nf90_def_var(ncid, trim(adjustl(varname(5)%c)), NF90_double, dimsid2D,  varid(5)) )
    call check( nf90_def_var(ncid, trim(adjustl(varname(6)%c)), NF90_double, dimsid2D,  varid(6)) )    
    call check( nf90_def_var(ncid, trim(adjustl(varname(7)%c)), NF90_double, dimsid3DV, varid(7)) )
    call check( nf90_def_var(ncid, trim(adjustl(varname(8)%c)), NF90_double, dimsid3DV, varid(8)) )    
    do i = 9, nvars
      call check( nf90_def_var(ncid, trim(adjustl(varname(i)%c)), NF90_double, dimsid3D, varid(i)) )
    enddo
    
    do i = 1, nvars
      call check( nf90_put_att(ncid, varid(i), "units", trim(adjustl(varunits(i)%c))) )
    enddo
    ! End define mode.
    call check( nf90_enddef(ncid) )
    ! call check( nf90_put_var(ncid, lvl_varid, [(i, i = 1, nz)]) )
    ! call check( nf90_put_var(ncid, col_varid, [(i, i = 1, nColumns)]) )
    ! call check( nf90_put_var(ncid, vert_varid, [(i, i = 1, nvert)]) )
    call check( nf90_close(ncid) )
  end subroutine writenc_setup

  subroutine netcdf_writevar(columns, t, dt, nz, ntracer)
    use class_parameter,        only: dimX, dimY
    use class_AMRTypes,         only: ColumnList, column
    type(ColumnList), intent(in) :: columns
    integer,          intent(in) :: t, ntracer, nz
    real(dp),         intent(in) :: dt
    type(Column),       pointer  :: MeshColumn
    integer                      :: i
 
    MeshColumn => columns%head_
    do i = 1, columns%n_
      i_(i)    = MeshColumn%prop_%i_(dimX)
      j_(i)    = MeshColumn%prop_%i_(dimY)
      l_(i)    = MeshColumn%l_
      lons(i) = MeshColumn%prop_%x_(dimX) + 0.5*MeshColumn%prop_%dx_(dimX)
      lats(i) = MeshColumn%prop_%x_(dimY) + 0.5*MeshColumn%prop_%dx_(dimY)
      verts(1, i, dimX) = MeshColumn%prop_%x_(dimX)
      verts(2, i, dimX) = MeshColumn%prop_%x_(dimX) + MeshColumn%prop_%dx_(dimX)
      verts(3, i, dimX) = MeshColumn%prop_%x_(dimX) + MeshColumn%prop_%dx_(dimX)
      verts(4, i, dimX) = MeshColumn%prop_%x_(dimX)
      verts(1, i, dimY) = MeshColumn%prop_%x_(dimY)
      verts(2, i, dimY) = MeshColumn%prop_%x_(dimY)
      verts(3, i, dimY) = MeshColumn%prop_%x_(dimY) + MeshColumn%prop_%dx_(dimY)
      verts(4, i, dimY) = MeshColumn%prop_%x_(dimY) + MeshColumn%prop_%dx_(dimY)
      q(1, :, i)        = MeshColumn%prop_%q_(1, :)
      ! q(2, :, i)        = MeshColumn%prop_%theta_
      MeshColumn => MeshColumn%next_
    enddo

    if (allocated(time)) then
      if (allocated(timeOld)) deallocate(timeOld)
      allocate(timeOld(it))
      timeOld(1:) = time(1:)
    endif
    if (allocated(time)) deallocate(time)
    allocate(time(it + 1))
    if (allocated(timeOld)) time(1:it) = timeOld(1:)

    ! call netcdf_time(timediff)
    time(it+1) = (t*dt)/(3600.*24.)
    it = it + 1
    call writenc_vars(ntracer + 8 + 1, columns%n_, nz)
  end subroutine netcdf_writevar

  subroutine writenc_vars(nvars, nColumns, nz)
    use class_parameter,             only: dimX, dimY
    use netcdf
    implicit none
    ! This is the name of the data file we will create.
    integer,             intent(in)    :: nvars, nz
    integer,             intent(in)    :: nColumns
    integer                            :: count1D(1), start1D(1)    
    integer                            :: count2D(2), start2D(2)
    integer                            :: count3D(3), start3D(3)
    integer                            :: count3DV(3), start3DV(3)
    ! Loop indices
    integer                            :: i
    call check( nf90_open(filename, nf90_write, ncid) )
    count1D = [1]; start1D = [1]

    count2D = [nColumnsTotal, 1]; start2D = [1, 1]

    count3D = [nz, nColumnsTotal, 1]; start3D = [1, 1, 1]

    count3DV = [4, nColumnsTotal, 1]; start3DV = [1, 1, 1]

    start3DV(3) = it; start3D(3)  = it
    start2D(2)  = it; start1D(1)  = it    
    ! write indexes
    call check( nf90_put_var(ncid, rec_varid, time) )
    call check( nf90_put_var(ncid, varid(1), [nColumns],          start = start1D,  count = count1D) )
    call check( nf90_put_var(ncid, varid(2), i_(1:),               start = start2D,  count = count2D) )
    call check( nf90_put_var(ncid, varid(3), j_(1:),               start = start2D,  count = count2D) )
    call check( nf90_put_var(ncid, varid(4), l_(1:),               start = start2D,  count = count2D) )
    call check( nf90_put_var(ncid, varid(5), lons(1:),            start = start2D,  count = count2D) )
    call check( nf90_put_var(ncid, varid(6), lats(1:),            start = start2D,  count = count2D) )
    call check( nf90_put_var(ncid, varid(7), verts(1:, 1:, dimX), start = start3DV, count = count3DV) )
    call check( nf90_put_var(ncid, varid(8), verts(1:, 1:, dimY), start = start3DV, count = count3DV) )
    do i = 9, nvars
      call check( nf90_put_var(ncid, varid(i), q(i - 8, 1:, 1:),  start = start3D,  count = count3D) )
    enddo
     call check( nf90_close(ncid) )
  end subroutine writenc_vars

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

  ! subroutine netcdf_time(time)
  !   USE mo_time_control,   ONLY: next_date, start_date, get_date_components
  !   real(dp), intent(out)     :: time
  !   integer                   :: start_year, start_month, start_day, start_hour, start_minute, start_second
  !   integer                   :: now_year, now_month, now_day, now_hour, now_minute, now_second
  !   integer                   :: year
  !   real(dp)                  :: start_day_time, now_day_time

  !   CALL get_date_components(start_date,start_year,start_month,start_day,start_hour,start_minute,start_second)
  !   CALL get_date_components(next_date,now_year,now_month,now_day,now_hour,now_minute,now_second)

  !   time = 0._dp
  !   do year = start_year, now_year - 1
  !     if (isLeapYear(year)) then
  !       time = time + 366._dp
  !     else
  !       time = time + 365._dp
  !     endif
  !   enddo

  !   if (start_year == now_year) then
  !     time = time + ( yearday(now_day, now_month, now_year) - yearday(start_day, start_month, start_year) )
  !   else
  !     time = time + yearday(now_day, now_month, now_year)
  !   endif

  !   start_day_time = 0._dp
  !   if ((start_year == now_year) .and. (start_month == now_month) .and. (start_day == now_day)) then
  !     start_day_time = start_hour/24._dp + start_minute/60._dp/24._dp + start_second/3600._dp/24._dp
  !   endif
  !   now_day_time = now_hour/24._dp + now_minute/60._dp/24._dp + now_second/3600._dp/24._dp
  !   time = time + (now_day_time - start_day_time)

  !   contains
  !     pure elemental logical function isLeapYear(year)
  !       integer,intent(in) :: year !! year

  !       isLeapYear = (mod(year,4) == 0 .and. .not. mod(year,100) == 0)&
  !               .or. (mod(year,400) == 0)
  !     end function isLeapYear

  !     pure elemental integer function yearday(day, month, year)
  !       integer, intent(in)  :: day, month, year

  !       integer :: imonth

  !       yearday = 0
  !       do imonth = 1, month-1
  !         yearday = yearday + daysInMonth(imonth, year)
  !       enddo
  !       yearday = yearday + day
  !     end function yearday

  !     pure elemental integer function daysInMonth(month,year)
  !       integer, intent(in) :: month !! month
  !       integer, intent(in) :: year  !! year
   
  !       integer, parameter,dimension(12) :: &
  !               days = [31,28,31,30,31,30,31,31,30,31,30,31]

  !       if(month < 1 .or. month > 12)then
  !         ! Should raise an error and abort here, however we want to keep
  !         ! the pure and elemental attributes. Make sure this function is 
  !         ! called with the month argument in range. 
  !         daysInMonth = 0
  !         return
  !       endif

  !       if(month == 2 .and. isLeapYear(year))then
  !         daysInMonth = 29
  !       else
  !         daysInMonth = days(month)
  !       endif
  !     end function daysInMonth
  ! end subroutine netcdf_time
end module class_netcdf