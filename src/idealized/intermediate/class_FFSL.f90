module class_FFSL
implicit none
private
public :: FFSL, update
contains
  subroutine FFSL(m, nx, ny, nz, ntracer, MaxRef, dt)
    ! flux form semi-Lagrangian scheme
    use class_mesh,           only: mesh
    use class_AMRTypes,       only: column
    use class_kind,           only: dp
    use class_column,         only: column_dA
    use class_parameter,      only: dimX, dimY, Inner, Outer
    use class_depart,         only: depart
    use class_CISL,           only: CISL
    ! use class_div, only: div
    type(mesh),            intent(inout) :: m(1:)
    real(dp),              intent(in)    :: dt
    integer,               intent(in)    :: nx, ny, nz, ntracer, MaxRef
    type(Column), pointer                :: MeshColumn
    real(dp)                             :: dA
    integer                              :: k, i

    do k = 1, nz
      ! initialize values
      MeshColumn => m(1)%columns_%head_
      do i = 1, m(1)%columns_%n_
        MeshColumn%prop_%qmid_(:, dimX, 0) = MeshColumn%prop_%q0_(:, k)
        MeshColumn%prop_%qmid_(:, dimY, 0) = MeshColumn%prop_%q0_(:, k)
        MeshColumn%prop_%area = 1._dp
        MeshColumn => MeshColumn%next_
      enddo

      ! one-dimensional CISL for inner operator
      call Depart(m, m(2)%faces_, dt, m(2)%coord_%dxMax_, nx, m(2)%coord_%x_, m(2)%coord_%y_, k, m(2)%coord_%l0_)
      call CISL(m(2)%columns_, m(1)%columns_, k, Inner, m(1)%coord_%l0_, MaxRef, nx, ny, m(1)%coord_%x_, m(1)%coord_%y_, ntracer)
      ! call div(m(2)%columns_, m(2)%columns_, k, Outer, m(1)%coord_%l0_, MaxRef, nx, ny, m(1)%coord_%x_, m(1)%coord_%y_, ntracer)
      ! inner operator splitting
      MeshColumn => m(2)%columns_%head_
      do i = 1, m(2)%columns_%n_
        dA = column_dA(MeshColumn)
        MeshColumn%prop_%qinter_(1, 1, dimX) = MeshColumn%prop_%qmid_(1, dimX, inner)/dA
        MeshColumn%prop_%qinter_(1, 1, dimY) = MeshColumn%prop_%qmid_(1, dimY, inner)/dA
        MeshColumn%prop_%qmid_(1, dimX, inner) = MeshColumn%prop_%q0_(1, k) + &
          0.5_dp*(MeshColumn%prop_%qmid_(1, dimX, inner) - MeshColumn%prop_%q0_(1, k)*MeshColumn%prop_%div_(dimX))/dA
        MeshColumn%prop_%qmid_(1, dimY, inner) = MeshColumn%prop_%q0_(1, k) + &
          0.5_dp*(MeshColumn%prop_%qmid_(1, dimY, inner) - MeshColumn%prop_%q0_(1, k)*MeshColumn%prop_%div_(dimY))/dA
        
        ! MeshColumn%prop_%qmid_(2:ntracer + 1, dimX, inner) = &
        !   MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%q0_(1, k) + &
        !   0.5_dp*(MeshColumn%prop_%qmid_(2:ntracer + 1, dimX, inner) - &
        !     MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%div_(dimX))/dA
        ! MeshColumn%prop_%qmid_(2:ntracer + 1, dimY, inner) = &
        !   MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%q0_(1, k) + &
        !   0.5_dp*(MeshColumn%prop_%qmid_(2:ntracer + 1, dimY, inner) - &
        !     MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%div_(dimY))/dA
        if (is_zero(MeshColumn, dimX, k)) then
          MeshColumn%prop_%qmid_(1, dimX, inner) = MeshColumn%prop_%q0_(1, k)
          ! MeshColumn%prop_%qmid_(2:ntracer + 1, dimX, inner) = &
          !   MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%q0_(1, k)
        endif
        if (is_zero(MeshColumn, dimY, k)) then
          MeshColumn%prop_%qmid_(1, dimY, inner) = MeshColumn%prop_%q0_(1, k)
          ! MeshColumn%prop_%qmid_(2:ntracer + 1, dimY, inner) = &
          !   MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%q0_(1, k)
        endif

        MeshColumn => MeshColumn%next_
      enddo
      ! pole average
      ! call PoleMeanMid(m(2)%columns_, ny, m(2)%coord_%l0_, ntracer)
      ! print *, "----------start outer ------------------"
      ! 1D CISL for outer operator
      call CISL(m(2)%columns_, m(2)%columns_, k, outer, m(2)%coord_%l0_, MaxRef, nx, ny, m(2)%coord_%x_, m(2)%coord_%y_, ntracer)

      ! final update for outer operator
      MeshColumn => m(2)%columns_%head_
      do i = 1, m(2)%columns_%n_
        dA = column_dA(MeshColumn)
        MeshColumn%prop_%qmid_(:, dimX, outer) = &
          MeshColumn%prop_%qmid_(:, dimX, outer)/dA - MeshColumn%prop_%qmid_(:, dimY, inner)
        MeshColumn%prop_%qmid_(:, dimY, outer) = &
          MeshColumn%prop_%qmid_(:, dimY, outer)/dA - MeshColumn%prop_%qmid_(:, dimX, inner)

        MeshColumn%prop_%q_(1, k) = &
          MeshColumn%prop_%q0_(1, k)+ MeshColumn%prop_%qmid_(1, dimX, outer) + MeshColumn%prop_%qmid_(1, dimY, outer)
        MeshColumn%prop_%qmid_(:, :, :) = 0._dp
        ! MeshColumn%prop_%q_(2:ntracer + 1, k) = MeshColumn%prop_%q0_(2:ntracer + 1, k)*MeshColumn%prop_%q0_(1, k) + &
        !   MeshColumn%prop_%qmid_(2:ntracer + 1, dimX, outer) + MeshColumn%prop_%qmid_(2:ntracer + 1, dimY, outer)
        ! mass = mass + MeshColumn%prop_%q_(1, k)*dA
        MeshColumn => MeshColumn%next_
      enddo
      ! pole average
      ! call PoleMeanFinal(m(2)%columns_, ny, ntracer, m(2)%coord_%l0_, k)
    enddo
    ! fct_x
    ! call fct_x(m(2)%columns_, ntracer, nz)
    ! vertical transport
  end subroutine FFSL

  subroutine update(m)
    ! time update
    use class_mesh,     only: mesh
    use class_AMRTypes, only: column
    type(mesh), intent(in):: m(1:)
    type(Column), pointer :: MeshColumn
    integer               :: i
    MeshColumn => m(2)%columns_%head_
    do i = 1, m(2)%columns_%n_
      MeshColumn%OtherMesh_%prop_%q0_(:, :) = MeshColumn%prop_%q_(:, :)
      MeshColumn%OtherMesh_%prop_%q_(:, :) = MeshColumn%prop_%q_(:, :)
      MeshColumn => MeshColumn%next_
    enddo
  end subroutine update

  function is_zero(MeshColumn, dim, k) result(val)
    ! check if the velocity is in opposite direction
    use class_AMRTypes,       only: face, column
    type(column), pointer, intent(in) :: MeshColumn
    integer,               intent(in) :: dim, k
    integer                           :: i, j
    type(face),   pointer             :: MeshFaceL, MeshFaceR
    logical                           :: val
    val = .false.
    MeshFaceL => MeshColumn%prop_%faces_(dim, 1)%head_
    do i = 1, MeshColumn%prop_%faces_(dim, 1)%n_
      MeshFaceR => MeshColumn%prop_%faces_(dim, 2)%head_
      do j = 1, MeshColumn%prop_%faces_(dim, 2)%n_
        if (MeshFaceR%prop_%u_(k)*MeshFaceL%prop_%u_(k) < 0) then
          val = .true.
          return
        endif
        MeshFaceR => MeshFaceR%nextc_(1)%p_
      enddo
      MeshFaceL => MeshFaceL%nextc_(2)%p_
    enddo
  end function is_zero

  ! subroutine PoleMeanMid(columns, ny, l0, ntracer)
  !   ! average pole values for intermediate values
  !   use class_AMRTypes,       only: ColumnList, column
  !   use class_kind,           only: dp
  !   use class_column,         only: column_dA
  !   use class_parameter,      only: dimX, dimY
  !   type(ColumnList),      intent(in) :: Columns
  !   integer,               intent(in) :: ny, l0, ntracer
  !   type(column), pointer             :: MeshColumn
  !   integer                           :: i, Lscale
  !   real(dp)                          :: dA
  !   real(dp)                          :: dA_n, dA_s, q_north(ntracer + 1, dimX:dimY, 2), q_south(ntracer + 1, dimX:dimY, 2)
  !   dA_n = 0._dp; dA_S = 0._dp
  !   q_north(:, :, :) = 0._dp; q_south(:, :, :) = 0._dp
  !   MeshColumn => Columns%head_
  !   do i = 1, Columns%n_
  !     dA = column_dA(MeshColumn)
  !     Lscale = 2**(MeshColumn%l_ - l0)
  !     if (MeshColumn%prop_%i_(dimY)/Lscale == ny - 1) then
  !       dA_N = dA_N + dA
  !       q_north(:, :, :) = q_north(:, :, :) + MeshColumn%prop_%qmid_(:, :, 1:)*dA
  !     else if (MeshColumn%prop_%i_(dimY)/Lscale == 0) then
  !       dA_S = dA_S + dA
  !       q_south(:, :, :) = q_south(:, :, :) + MeshColumn%prop_%qmid_(:, :, 1:)*dA
  !     endif
  !     MeshColumn => MeshColumn%next_
  !   enddo
  !   q_north(:, :, :) = q_north(:, :, :)/dA_N
  !   q_south(:, :, :) = q_south(:, :, :)/dA_S

  !   MeshColumn => Columns%head_
  !   do i = 1, Columns%n_
  !     Lscale = 2**(MeshColumn%l_ - l0)
  !     if (MeshColumn%prop_%i_(dimY)/Lscale == ny - 1) then
  !       MeshColumn%prop_%qmid_(:, :, 1:) = q_north(:, :, :)
  !     else if (MeshColumn%prop_%i_(dimY)/Lscale == 0) then
  !       MeshColumn%prop_%qmid_(:, :, 1:) = q_south(:, :, :)
  !     endif
  !     MeshColumn => MeshColumn%next_
  !   enddo
  ! end subroutine PoleMeanMid

  ! subroutine PoleMeanFinal(columns, ny, ntracer, l0, k)
  !   ! average pole values in the end
  !   use class_AMRTypes,       only: ColumnList, column
  !   use class_kind,           only: dp
  !   use class_column,         only: column_dA
  !   use class_parameter,      only: dimY
  !   type(ColumnList),      intent(in) :: Columns
  !   integer,               intent(in) :: ny, ntracer, l0, k
  !   type(column), pointer             :: MeshColumn
  !   integer                           :: i, Lscale
  !   real(dp)                          :: dA
  !   real(dp)                          :: dA_n, dA_s, q_north(ntracer + 1), q_south(ntracer + 1)

  !   dA_n = 0._dp; dA_S = 0._dp
  !   q_north(:) = 0._dp; q_south(:) = 0._dp
  !   MeshColumn => Columns%head_
  !   do i = 1, Columns%n_
  !     dA = column_dA(MeshColumn)
  !     Lscale = 2**(MeshColumn%l_ - l0)
  !     if (MeshColumn%prop_%i_(dimY)/Lscale == ny - 1) then
  !       dA_N = dA_N + dA
  !       q_north(:) = q_north(:) + MeshColumn%prop_%q_(:, k)*dA
  !     else if (MeshColumn%prop_%i_(dimY)/Lscale == 0) then
  !       dA_S = dA_S + dA
  !       q_south(:) = q_south(:) + MeshColumn%prop_%q_(:, k)*dA
  !     endif
  !     MeshColumn => MeshColumn%next_
  !   enddo
  !   q_north(:) = q_north(:)/dA_N
  !   q_south(:) = q_south(:)/dA_S

  !   MeshColumn => Columns%head_
  !   do i = 1, Columns%n_
  !     Lscale = 2**(MeshColumn%l_ - l0)
  !     if (MeshColumn%prop_%i_(dimY)/Lscale == ny - 1) then
  !       MeshColumn%prop_%q_(:, k) = q_north(:)
  !     else if (MeshColumn%prop_%i_(dimY)/Lscale == 0) then
  !       MeshColumn%prop_%q_(:, k) = q_south(:)
  !     endif
  !     MeshColumn => MeshColumn%next_
  !   enddo
  ! end subroutine PoleMeanFinal

  ! subroutine fct_x(columns, ntracer, nz)
  !   ! move negative values to adjacent cells
  !   use class_AMRTypes,       only: ColumnList, ColumnArray, face, column
  !   use class_kind,           only: dp
  !   use class_column,         only: column_dA
  !   use class_parameter,      only: dimX
  !   type(ColumnList), intent(in) :: columns
  !   integer,          intent(in) :: nz, ntracer
  !   type(ColumnArray)            :: AdjColumns(2)
  !   type(face),   pointer        :: MeshFace
  !   type(column), pointer        :: MeshColumn
  !   real(dp)                     :: d0, d1, dA, qsum, dAtmp
  !   integer                      :: ipx, i, j, itracer, k, ip, indx
  !   REAL(dp), PARAMETER          :: tiny = 1.0e-40_dp
  !   ! left
  !   ipx = 0
  !   MeshColumn => columns%head_
  !   do i = 1, columns%n_
  !     do k = 1, nz
  !       do itracer = 2, ntracer + 1
  !         if (MeshColumn%prop_%q_(itracer, k) < 0._dp) then
  !           do indx = 1, 2
  !             ipx = 1; qsum = 0._dp; dA = 0._dp
  !             MeshFace => MeshColumn%prop_%faces_(dimX, indx)%head_
  !             do j = 1, MeshColumn%prop_%faces_(dimX, indx)%n_
  !               AdjColumns(j)%p_ => MeshFace%prop_%columns_(1)%p_
  !               dAtmp = column_dA(AdjColumns(j)%p_)
  !               qsum = qsum + AdjColumns(j)%p_%prop_%q_(itracer, k)*dAtmp
  !               dA = dA + dAtmp
  !               MeshFace => MeshFace%nextc_(3 - indx)%p_
  !             enddo
  !             d0 = max(0._dp, qsum/dA)
  !             d1 = min(-MeshColumn%prop_%q_(itracer, k), d0)
  !             if (abs(qsum) >  epsilon(1._dp)) then
  !               dA = column_dA(AdjColumns(1)%p_)
  !               AdjColumns(1)%p_%prop_%q_(itracer, k) = &
  !                 AdjColumns(1)%p_%prop_%q_(itracer, k) - AdjColumns(1)%p_%prop_%q_(itracer, k)*d1*dA/qsum
  !               if (MeshColumn%prop_%faces_(dimX, 1)%n_ == 2) then
  !                 dA = column_dA(AdjColumns(2)%p_)
  !                 AdjColumns(2)%p_%prop_%q_(itracer, k) = &
  !                   AdjColumns(2)%p_%prop_%q_(itracer, k) - AdjColumns(2)%p_%prop_%q_(itracer, k)*d1*dA/qsum
  !               endif
  !             endif
  !             MeshColumn%prop_%q_(itracer, k) = MeshColumn%prop_%q_(itracer, k) + d1
  !           enddo
  !           MeshColumn%prop_%q_(itracer, k) = MeshColumn%prop_%q_(itracer, k) + d1 + tiny
  !         endif
  !       enddo
  !     enddo
  !     MeshColumn => MeshColumn%next_
  !   enddo

  !   if (ipx == 0) return
  !   do ip = 1, 2
  !     MeshColumn => columns%head_
  !     do i = 1, columns%n_
  !       do k = 1, nz
  !         do itracer = 2, ntracer + 1
  !           if (MeshColumn%prop_%q_(itracer, k) < 0._dp) then
  !             qsum = 0._dp; dA = 0._dp
  !             MeshFace => MeshColumn%prop_%faces_(dimX, 1)%head_
  !             do j = 1, MeshColumn%prop_%faces_(dimX, 1)%n_
  !               AdjColumns(j)%p_ => MeshFace%prop_%columns_(1)%p_
  !               dA= column_dA(AdjColumns(j)%p_)
  !               qsum = qsum + AdjColumns(j)%p_%prop_%q_(itracer, k)*dA
  !               MeshFace => MeshFace%nextc_(2)%p_
  !             enddo
  !             if (abs(qsum) >  epsilon(1._dp)) then
  !               dA = column_dA(AdjColumns(1)%p_)
  !               AdjColumns(1)%p_%prop_%q_(itracer, k) = &
  !                 AdjColumns(1)%p_%prop_%q_(itracer, k) - &
  !                   AdjColumns(1)%p_%prop_%q_(itracer, k)*MeshColumn%prop_%q_(itracer, k)*dA/qsum
  !               if (MeshColumn%prop_%faces_(dimX, 1)%n_ == 2) then
  !                 dA = column_dA(AdjColumns(2)%p_)
  !                 AdjColumns(2)%p_%prop_%q_(itracer, k) = &
  !                   AdjColumns(2)%p_%prop_%q_(itracer, k) - &
  !                     AdjColumns(2)%p_%prop_%q_(itracer, k)*MeshColumn%prop_%q_(itracer, k)*dA/qsum
  !               endif
  !             endif
  !             MeshColumn%prop_%q_(itracer, k) = 0._dp
  !           endif
  !         enddo
  !       enddo
  !       MeshColumn => MeshColumn%next_
  !     enddo
  !   enddo
  ! end subroutine fct_x
end module class_FFSL