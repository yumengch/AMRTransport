module class_CISL
implicit none
private
public :: CISL
contains
  subroutine CISL(columns, DepartColumns, klev, InOrOut, l0, MaxRef, nx, ny, x, y, ntracer)
    use class_kind,      only: dp
    use class_AMRTypes,       only: ColumnList, ColumnStack, column, face, FaceArray
    use class_parameter,      only: dimX, dimY, Inner, Outer, pi
    use class_column,    only: column_dA
    type(ColumnList), intent(in) :: columns, DepartColumns
    integer,          intent(in) :: klev, InOrOut, l0, MaxRef
    integer,          intent(in) :: nx, ny, ntracer
    real(dp),         intent(in) :: x(0:), y(0:)
    type(ColumnStack)            :: UnCheckedColumns
    type(column), pointer        :: ArrivalColumn, MeshColumn
    type(face),   pointer        :: MeshFace, DepartFace, AdjMeshFace
    type(FaceArray)              :: DepartFaces(2)
    real(dp)                     :: xd(2, 2), u
    real(dp)                     :: dx, cosp
    integer                      :: i, dim, indx, iFace, j
    integer                      :: indxLimit(2), jLimitFurther, lLimitFurther
    integer                      :: lLimit, jLimit, Lscale
    logical                      :: lstartColumn, lpass_left, lpass_right, lFound, loverlap
    logical                      :: ladd, lwind, InThisColumn(2)

    do dim = dimX, dimY
      ArrivalColumn => columns%head_
      do i = 1, columns%n_
        xd(:, :) = 4*pi
        ! get the departure positions and select the closer faces of the lagrangian volume
        ArrivalColumn%prop_%qmid_(:, dim, InOrOut) = 0._dp
        ArrivalColumn%prop_%div_(dim)              = 0._dp
        do indx = 1, 2
          MeshFace => ArrivalColumn%prop_%faces_(dim, indx)%head_
          indxLimit(indx) = 1
          cosp = 1._dp
          if (dim == dimX) cosp = &
            (sin(MeshFace%prop_%y_(2)) - sin(MeshFace%prop_%y_(1)))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
          u            = MeshFace%prop_%u_(klev)/cosp
          do iFace = 1, ArrivalColumn%prop_%faces_(dim, indx)%n_
            xd(indx, iFace) = MeshFace%prop_%xd_
            if (dim == dimX) cosp = &
              (sin(MeshFace%prop_%y_(2)) - sin(MeshFace%prop_%y_(1)))/(MeshFace%prop_%y_(2) - MeshFace%prop_%y_(1))
            lwind = MeshFace%prop_%u_(klev)/cosp < u
            if (indx == 2) lwind = MeshFace%prop_%u_(klev)/cosp > u
            if (indx == 2) DepartFaces(iFace)%p_ => MeshFace
            if (lwind) then
              indxLimit(indx) = iFace
              u = MeshFace%prop_%u_(klev)/cosp
            endif      
            MeshFace => MeshFace%nextc_(3 - indx)%p_
          enddo
        enddo

        ! additionally get jLimitFurther for the farther faces of the jLimit 
        MeshFace => ArrivalColumn%prop_%faces_(dim, 2)%head_
        if (indxLimit(2) == 1) MeshFace => ArrivalColumn%prop_%faces_(dim, 2)%tail_
        jLimitFurther = MeshFace%prop_%i_( 3 - dim)
        lLimitFurther = MeshFace%prop_%l_

        ! get unchecked columns and mark the Lagrangian Volumes
        MeshFace => ArrivalColumn%prop_%faces_(dim, 1)%head_
        do iFace = 1, ArrivalColumn%prop_%faces_(dim, 1)%n_
          ! Search of Lagrangian Volume always starts from the faces with smaller index
          call overlap_GetUnchecked(MeshFace, ArrivalColumn, InOrOut, dim, l0, MaxRef, UnCheckedColumns)
          jLimit = MeshFace%prop_%i_(3 - dim)
          lLimit = MeshFace%prop_%l_
          do while (associated(UnCheckedColumns%head_))
            MeshColumn => UnCheckedColumns%head_
            UnCheckedColumns%head_ => MeshColumn%nextStack_
            MeshColumn%nextStack_  => null()
            !---------------------------------------------
            ! check status of faces in current column
            ! --------------------------------------------
            ! check if the current colum is the starting column with starting face
            lstartColumn = xd(1, iFace) < MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim) .and. &
                                          xd(1, iFace) > MeshColumn%prop_%x_(dim)
            ! if the departure point is at the edge of the column
            lstartColumn = lstartColumn .or. &
              abs(MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim) - xd(1, iFace)) <= 4._dp*epsilon(1._dp)
            if (lstartColumn) then
              jLimit = MeshFace%prop_%i_(3 - dim)
              lLimit = MeshFace%prop_%l_
            endif
            ! check if we passed the first closer faces
            ! if we passed, change the limit of the Lagrangian volume to the same limit as the arrival column
            dx = MeshColumn%prop_%x_(dim) - xd(1, indxLimit(1))
            lpass_left = ( dx > 0 .and. abs(dx) < 2*pi - abs(dx) ) .or. &
                         (dx < 0  .and. abs(dx) > 2*pi - abs(dx)) .or. abs(dx) <= epsilon(1._dp)
            if (lpass_left) then
              jLimit = ArrivalColumn%prop_%i_(3 - dim)
              lLimit = ArrivalColumn%l_
            endif
            dx = MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim) - xd(2, indxLimit(2)) 
            lpass_right = ( dx > 0 .and. abs(dx) < 2*pi - abs(dx) ) .or. (dx < 0  .and. abs(dx) > 2*pi - abs(dx))
            lpass_right = lpass_right .or. abs(dx) <= epsilon(1._dp)

            ! check if the first LV is filled
            lFound = .false.
            indx = 1
            if (associated(MeshColumn%prop_%LV_(dim, 1)%ArrivalColumn, ArrivalColumn)) indx = 2

            ! check if the closer face is in the column
            InThisColumn(:) = .false.
            do j = 1, 2
              if (.not. associated(DepartFaces(j)%p_)) cycle            
              InThisColumn(j) = xd(2, j) < MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim) .and. &
                                            xd(2, j) > MeshColumn%prop_%x_(dim)
              InThisColumn(j) = InThisColumn(j) .or. &
              abs(MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim) - xd(2, j)) <= epsilon(1._dp)
              InThisColumn(j) = InThisColumn(j) .or. &
              abs(MeshColumn%prop_%x_(dim) - xd(2, j)) <= epsilon(1._dp)
              DepartFace => DepartFaces(j)%p_
              Lscale = 2**abs(DepartFace%prop_%l_ - MeshColumn%l_)
              loverlap = DepartFace%prop_%i_(3 - dim)/Lscale == MeshColumn%prop_%i_(3 - dim)
              if (DepartFace%prop_%l_ < MeshColumn%l_) &
                loverlap = DepartFace%prop_%i_(3 - dim) == MeshColumn%prop_%i_(3 - dim)/Lscale
              InThisColumn(j) = InThisColumn(j) .and. loverlap
            enddo
            !-------------------------------------------
            ! form Lagrangian volumes
            ! ------------------------------------------
            ! mark the closer larger index face
            if (InThisColumn(indxLimit(2))) then
              DepartFace => DepartFaces(indxLimit(2))%p_
              Lscale = 2**abs(DepartFace%prop_%l_ - lLimit)
              loverlap = DepartFace%prop_%i_(3 - dim)/Lscale == jLimit
              if (DepartFace%prop_%l_ < lLimit) loverlap = DepartFace%prop_%i_(3 - dim)== jLimit/Lscale
              if (loverlap) then
                MeshColumn%prop_%LV_(dim, indx)%jLimit = max(jLimit, MeshColumn%prop_%i_(3 - dim), DepartFace%prop_%i_(3 - dim))
                MeshColumn%prop_%LV_(dim, indx)%lLimit = max(lLimit, MeshColumn%l_, DepartFace%prop_%l_)
                MeshColumn%prop_%LV_(dim, indx)%x_(1)  = MeshColumn%prop_%x_(dim)
                if (lstartColumn) MeshColumn%prop_%LV_(dim, indx)%x_(1)  = xd(1, iFace)
                MeshColumn%prop_%LV_(dim, indx)%x_(2)  = xd(2, indxLimit(2))
                MeshColumn%prop_%LV_(dim, indx)%ArrivalColumn => ArrivalColumn
                lfound = .true.
                indx = indx + 1
              endif
            endif

            ! check if we have passed the closer larger index
            if ( (lpass_right) .and. (lLimit == ArrivalColumn%l_) .and. (ArrivalColumn%prop_%faces_(dim, 2)%n_ > 1)  ) then
              jlimit = jLimitFurther
              lLimit = lLimitFurther
              ! set false
              if ((MeshColumn%l_ <= ArrivalColumn%l_) .and. (.not. InThisColumn(3 - indxLimit(2)) )) &
                lFound = .false.
            endif

            ! check the farthest face
            if (associated(MeshColumn%prop_%LV_(dim, 1)%ArrivalColumn, ArrivalColumn)) indx = 2
            ! marker if the farther face is in the column
            if (InThisColumn(3 - indxLimit(2))) then
              DepartFace => DepartFaces(3 - indxLimit(2))%p_
              Lscale = 2**abs(DepartFace%prop_%l_ - lLimit)
              loverlap = DepartFace%prop_%i_(3 - dim)/Lscale == jLimit
              if (DepartFace%prop_%l_ < lLimit) loverlap = DepartFace%prop_%i_(3 - dim)== jLimit/Lscale
              if (loverlap) then              
                MeshColumn%prop_%LV_(dim, indx)%jLimit = max(jLimit, MeshColumn%prop_%i_(3 - dim), DepartFace%prop_%i_(3 - dim))
                MeshColumn%prop_%LV_(dim, indx)%lLimit = max(lLimit, MeshColumn%l_, DepartFace%prop_%l_)
                MeshColumn%prop_%LV_(dim, indx)%x_(1)  = MeshColumn%prop_%x_(dim)
                if (lstartColumn) MeshColumn%prop_%LV_(dim, indx)%x_(1)  = xd(1, iFace)
                MeshColumn%prop_%LV_(dim, indx)%x_(2)  = xd(2, 3 - indxLimit(2))
                MeshColumn%prop_%LV_(dim, indx)%ArrivalColumn => ArrivalColumn
                lfound = .true.
              endif
            endif

            ! if there is no face in the column
            if (.not. lfound) then
              MeshColumn%prop_%LV_(dim, indx)%jLimit = max(jLimit, MeshColumn%prop_%i_(3 - dim))
              MeshColumn%prop_%LV_(dim, indx)%lLimit = max(lLimit, MeshColumn%l_)
              MeshColumn%prop_%LV_(dim, indx)%x_(1)  = MeshColumn%prop_%x_(dim)
              if (lstartColumn) MeshColumn%prop_%LV_(dim, indx)%x_(1)  = xd(1, iFace)
              MeshColumn%prop_%LV_(dim, indx)%x_(2)  = MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim)
              MeshColumn%prop_%LV_(dim, indx)%ArrivalColumn => ArrivalColumn

              ! Find the adjacent column for the arrival column
              AdjMeshFace => MeshColumn%prop_%faces_(dim, 2)%head_
              do j = 1, MeshColumn%prop_%faces_(dim, 2)%n_
                ! check if the face is overlap with the lagrangian volume
                Lscale = 2**abs(AdjMeshFace%prop_%l_ - lLimit)
                loverlap = AdjMeshFace%prop_%i_( 3 - dim) == jLimit/Lscale
                if (AdjMeshFace%prop_%l_ > lLimit) loverlap = AdjMeshFace%prop_%i_( 3 - dim)/Lscale == jLimit
                ! if LV1 is in this column and LV2 is not present
                ladd = .true.
                if (associated(AdjMeshFace%prop_%columns_(2)%p_%prop_%LV_(dim, 1)%ArrivalColumn, ArrivalColumn)) then
                  ladd = AdjMeshFace%prop_%columns_(2)%p_%prop_%LV_(dim, 1)%lLimit >= lLimit
                  Lscale = 2**(AdjMeshFace%prop_%columns_(2)%p_%prop_%LV_(dim, 1)%lLimit - lLimit)
                  ladd = ladd .and. AdjMeshFace%prop_%columns_(2)%p_%prop_%LV_(dim, 1)%jLimit/Lscale /= jLimit .and. &
                    .not. associated(AdjMeshFace%prop_%columns_(2)%p_%prop_%LV_(dim, 2)%ArrivalColumn, ArrivalColumn)
                endif
                if ( (loverlap) .and. (ladd) )then
                  if (associated(MeshFace%prop_%columns_(2)%p_%nextStack_)) then
                    print *, "problematic implementation of LV for the next"
                  endif
                  AdjMeshFace%prop_%columns_(2)%p_%nextStack_ => UnCheckedColumns%head_     
                  UnCheckedColumns%head_ => AdjMeshFace%prop_%columns_(2)%p_
                endif
                AdjMeshFace => AdjMeshFace%nextc_(1)%p_
              enddo
            endif
        ! if ((ArrivalColumn%prop_%i_(1) == 22).and.(ArrivalColumn%prop_%i_(2) == 2).and. (ArrivalColumn%l_ == 8)) then
        !   print *, "----------------------------------------"
        !   print *, dim
        !   print *, MeshColumn%prop_%i_, MeshColumn%l_
        !   print *, ArrivalColumn%prop_%faces_(dim, 2)%n_
        !   print *, MeshFace%prop_%i_, MeshFace%prop_%l_
        !   print *, xd(2, :)
        !   print *, MeshColumn%prop_%x_(dim), MeshColumn%prop_%x_(dim) + MeshColumn%prop_%dx_(dim)
        !   print *, indxLimit(2)
        !   print *, jLimitFurther, lLimitFurther
        !   print *, jLimit, lLimit
        !   print *, indx, lfound, InThisColumn
        ! endif
            call CISL_stencil(MeshColumn, ArrivalColumn, DepartColumns, lFound,&
                               lstartColumn, x, y, nx, ny, ntracer, dim, l0, InOrOut)
          enddo
          MeshFace => MeshFace%nextc_(2)%p_
        enddo
        DepartFaces(1)%p_ => null()
        DepartFaces(2)%p_ => null()
        ArrivalColumn => ArrivalColumn%next_
      enddo
      if (InOrOut == Outer) then
        MeshColumn => Departcolumns%head_
        do i = 1, Departcolumns%n_
          MeshColumn%prop_%has_dim(:) = .false.
          if (allocated(MeshColumn%prop_%qghost_)) deallocate(MeshColumn%prop_%qghost_)
          MeshColumn%prop_%nGhost_ = 0
          MeshColumn%prop_%LV_(dimX, 1)%ArrivalColumn => null()
          MeshColumn%prop_%LV_(dimX, 2)%ArrivalColumn => null()
          MeshColumn%prop_%LV_(dimY, 1)%ArrivalColumn => null()
          MeshColumn%prop_%LV_(dimY, 2)%ArrivalColumn => null()
          MeshColumn => MeshColumn%next_
        enddo
      endif
    enddo
    if (InOrOut == Inner) then
      MeshColumn => Departcolumns%head_
      do i = 1, Departcolumns%n_
        MeshColumn%prop_%has_dim(:) = .false.
        if (allocated(MeshColumn%prop_%qghost_)) deallocate(MeshColumn%prop_%qghost_)
        MeshColumn%prop_%nGhost_ = 0
        MeshColumn%prop_%LV_(dimX, 1)%ArrivalColumn => null()
        MeshColumn%prop_%LV_(dimX, 2)%ArrivalColumn => null()
        MeshColumn%prop_%LV_(dimY, 1)%ArrivalColumn => null()
        MeshColumn%prop_%LV_(dimY, 2)%ArrivalColumn => null()        
        MeshColumn => MeshColumn%next_
      enddo 

      ! MeshColumn => Departcolumns%head_
      ! do i = 1, Departcolumns%n_
      !   ! if (abs(MeshColumn%prop_%area(dimX))> 1e-13) then
      !   !   print *, MeshColumn%prop_%i_, MeshColumn%l_, MeshColumn%prop_%area(dimX)
      !   ! endif
      !   if (abs(MeshColumn%prop_%area(dimY))> 1e-13) then
      !     print *, MeshColumn%prop_%i_, MeshColumn%l_, MeshColumn%prop_%area(dimX)
      !   endif        
      !   MeshColumn => MeshColumn%next_
      ! enddo
    endif

  end subroutine CISL

  subroutine overlap_GetUnchecked(MeshFace, ArrivalColumn, InOrOut, dim, l0, MaxRef, UnCheckedColumns)
    use class_kind,      only: dp
    use class_AMRTypes,       only: ColumnStack, column, face
    use class_parameter,      only: dimX, dimY        
    type(face),        pointer, intent(in)    :: MeshFace
    type(column),      pointer, intent(in)    :: ArrivalColumn
    integer,                    intent(in)    :: InOrOut, l0, MaxRef, dim
    type(ColumnStack),          intent(inout) :: UnCheckedColumns
    type(ColumnStack)                         :: StackColumns(MeshFace%prop_%l_:l0 + MaxRef)
    type(column),      pointer                :: StackColumn
    integer                                   :: l, i, ii(dimX:dimY), Lscale
    real(dp)                                  :: dx, x0
    logical                                   :: ladd
    ! put departure face into departure column
    StackColumns(MeshFace%prop_%l_)%head_ => MeshFace%prop_%DepartColumns_(InOrOut)%p_
    ! if the departure column is not active, get the possible scale factors
    dx = MeshFace%prop_%Departdx_(dim)
    x0 = MeshFace%prop_%Departx0_(dim)
    do l = MeshFace%prop_%l_, l0 + MaxRef
      ii(dim) = 0
      ! decide if the children of the departure cell
      if ((MeshFace%prop_%xd_ - x0)/dx > 0.5_dp) then
        ii(dim) = 1
        x0 = x0 + dx/2._dp
      endif
      dx = dx/2._dp

      do while (associated(StackColumns(l)%head_))
        StackColumn => StackColumns(l)%head_
        StackColumns(l)%head_ => StackColumn%nextStack_
        StackColumn%nextStack_ => null()
        if (associated(StackColumn%Child_)) then
          do i = 0, 1
            ii( 3 - dim) = i
            ! if there is still child, put the column into next level stack
            StackColumn%Child_(ii(dimX), ii(dimY))%nextStack_ => StackColumns(l + 1)%head_
            StackColumns(l + 1)%head_ => StackColumn%Child_(ii(dimX), ii(dimY))
          enddo
        else
          ladd = .true.
          ! put the face into column
          if (associated(StackColumn%prop_%LV_(dim, 1)%ArrivalColumn, ArrivalColumn)) then
            ladd = StackColumn%prop_%LV_(dim, 1)%lLimit >= MeshFace%prop_%l_
            Lscale = 2**max(StackColumn%prop_%LV_(dim, 1)%lLimit - MeshFace%prop_%l_, 0)         
            ladd = ladd .and. StackColumn%prop_%LV_(dim, 1)%jLimit/Lscale /= MeshFace%prop_%i_(3 - dim)
            ladd = ladd .and. .not. associated(ArrivalColumn, StackColumn%prop_%LV_(dim, 2)%ArrivalColumn)
          endif
          if (ladd) then
            StackColumn%nextStack_ => UnCheckedColumns%head_
            UnCheckedColumns%head_ => StackColumn
          endif
        endif
      enddo
    enddo
  end subroutine overlap_GetUnchecked

  subroutine CISL_stencil(MeshColumn, ArrivalColumn, columns, lFound, lstartColumn, x, y, nx, ny, ntracer, dim, l0, InOrOut)
    use class_kind,      only: dp
    use class_AMRTypes,       only: column, ColumnList
    use class_parameter,      only: dimX, dimY, Inner
    use class_column,    only: column_dA
    type(ColumnList),      intent(in)    :: columns
    type(column), pointer, intent(inout) :: MeshColumn, ArrivalColumn
    logical,               intent(in)    :: lFound, lstartColumn
    real(dp),              intent(in)    :: x(0:), y(0:)
    integer,               intent(in)    :: nx, ny, ntracer, dim, l0, InOrOut
    integer                              :: i
    real(dp)                             :: dA
    real(dp)                             :: ax(ntracer + 1), ay(ntracer + 1)
    real(dp)                             :: bx(ntracer + 1), by(ntracer + 1), xLV(2), yLV(2)
    logical                              :: skip(2)
    ! get area of the mesh
    dA = column_dA(MeshColumn)
    skip(:) = .false.
    do i = 1, 2
      if (.not. associated(ArrivalColumn, MeshColumn%prop_%LV_(dim, i)%ArrivalColumn)) skip(i) = .true.
      if (abs(MeshColumn%prop_%LV_(dim, i)%x_(2) - MeshColumn%prop_%LV_(dim, i)%x_(1)) <epsilon(1._dp)) &
        skip(i) = .true.
    enddo
    ! compute divergence
    if (InOrOut == Inner) then
      do i = 1, 2
        if (skip(i)) cycle
        if (dim == dimY) then
          yLV(:) = MeshColumn%prop_%LV_(dim, i)%x_(:)
          xLV(:) = get_range(MeshColumn%prop_%LV_(dim, i)%jLimit, MeshColumn%prop_%LV_(dim, i)%lLimit, x, l0)
        else
          xLV(:) = MeshColumn%prop_%LV_(dim, i)%x_(:)
          yLV(:) = get_range(MeshColumn%prop_%LV_(dim, i)%jLimit, MeshColumn%prop_%LV_(dim, i)%lLimit, y, l0)
        endif
        ! scale the lagrangian volume into [-0.5, 0.5]*[-0.5, 0.5] as reference column
        xLV(:) = (xLV(:) - MeshColumn%prop_%x_(dimX))/MeshColumn%prop_%dx_(dimX) - 0.5_dp
        yLV(:) = ( sin(yLV) - sin(MeshColumn%prop_%x_(dimY)) )/(sin(MeshColumn%prop_%x_(dimY) &
          + MeshColumn%prop_%dx_(dimY)) - sin(MeshColumn%prop_%x_(dimY))) - 0.5_dp      
        ArrivalColumn%prop_%div_( dim ) = ArrivalColumn%prop_%div_( dim ) + (xLV(2) - xLV(1) ) * (yLV(2) - yLV(1))*dA
        MeshColumn%prop_%area(dim) = MeshColumn%prop_%area(dim) -  (xLV(2) - xLV(1) ) * (yLV(2) - yLV(1))
        ! if ((ArrivalColumn%prop_%i_(1) == 22).and.(ArrivalColumn%prop_%i_(2) == 2).and. (ArrivalColumn%l_ == 8)) then
        !   print *, dim
        !   print *, MeshColumn%prop_%i_, MeshColumn%l_
        !   print *, XLV
        !   print *, yLV
        ! endif
        ! if ((MeshColumn%prop_%i_(1) == 22).and.(MeshColumn%prop_%i_(2) == 2).and. (MeshColumn%l_ == 8)) then
        !   print *, dim
        !   print *, ArrivalColumn%prop_%i_, ArrivalColumn%l_
        !   print *, XLV
        !   print *, yLV
        ! endif
      enddo
    endif
    ! make integrals
    do i = 1, 2
      if (skip(i)) cycle
      ! 2-D reconstruction function
      if (MeshColumn%prop_%LV_(dim, i)%lLimit > MeshColumn%l_) then
        MeshColumn%prop_%dimStart_ = dimX
        MeshColumn%prop_%dimEnd_   = dimX
        MeshColumn%prop_%StencilStart_(dimX:) = -2
        MeshColumn%prop_%StencilEnd_(dimX:)   = 2      
        if ( (.not. MeshColumn%prop_%has_dim(dimX)) .and. (.not. MeshColumn%prop_%has_dim(dimY)) ) then
          MeshColumn%prop_%dimEnd_   = dimY
        else if ( (MeshColumn%prop_%has_dim(dimX)) .and. (.not. MeshColumn%prop_%has_dim(dimY)) ) then
          MeshColumn%prop_%dimStart_ = dimY
          MeshColumn%prop_%dimEnd_   = dimY
        else if ( (.not. MeshColumn%prop_%has_dim(dimX)) .and. (MeshColumn%prop_%has_dim(dimY)) ) then
          MeshColumn%prop_%dimEnd_   = dimY
        else
          MeshColumn%prop_%dimEnd_   = dimX
          MeshColumn%prop_%StencilStart_(dimX) = 0
          MeshColumn%prop_%StencilEnd_(dimX)   = 0
        endif
        call CISL_integral(MeshColumn, columns, nx, ny, ntracer, dim, l0, InOrOut)
        ! integrate inside the current stack column with picewise parabolic method
        ax = MeshColumn%prop_%a_(dimX, :)
        bx = MeshColumn%prop_%b_(dimX, :)
        ay = MeshColumn%prop_%a_(dimY, :)
        by = MeshColumn%prop_%b_(dimY, :)
        ! compute the range of the current lagrangian volume
        if (dim == dimY) then
          yLV(:) = MeshColumn%prop_%LV_(dim, i)%x_(:)
          xLV(:) = get_range(MeshColumn%prop_%LV_(dim, i)%jLimit, MeshColumn%prop_%LV_(dim, i)%lLimit, x, l0)
        else
          xLV(:) = MeshColumn%prop_%LV_(dim, i)%x_(:)
          yLV(:) = get_range(MeshColumn%prop_%LV_(dim, i)%jLimit, MeshColumn%prop_%LV_(dim, i)%lLimit, y, l0)
        endif
        ! scale the lagrangian volume into [-0.5, 0.5]*[-0.5, 0.5] as reference column
        xLV(:) = (xLV(:) - MeshColumn%prop_%x_(dimX))/MeshColumn%prop_%dx_(dimX) - 0.5_dp
        yLV(:) = ( sin(yLV) - sin(MeshColumn%prop_%x_(dimY)) )/(sin(MeshColumn%prop_%x_(dimY) &
          + MeshColumn%prop_%dx_(dimY)) - sin(MeshColumn%prop_%x_(dimY))) - 0.5_dp
        ! update thMeshColumne value at arrival column
        ArrivalColumn%prop_%qmid_(:, dim, InOrOut) = ArrivalColumn%prop_%qmid_(:, dim, InOrOut) + &
          (PPM2D_remap(xLV(2), yLV, ax, ay, bx, by, &
            MeshColumn%prop_%qmid_(:, 3 - dim, InOrOut - 1), ntracer) &
            - PPM2D_remap(xLV(1), yLV, ax, ay, bx, by, &
              MeshColumn%prop_%qmid_(:, 3 - dim, InOrOut - 1), ntracer))*dA
      else if ( (.not. lFound) .and. (.not. lstartColumn) ) then
      ! No need to use PPM
        ArrivalColumn%prop_%qmid_(:, dim, InOrOut) = ArrivalColumn%prop_%qmid_(:, dim, InOrOut) + &
                                                     MeshColumn%prop_%qmid_(:, 3 - dim, InOrOut - 1)*dA
      else
      ! 1-D reconstruction function
        MeshColumn%prop_%dimStart_ = dim
        MeshColumn%prop_%dimEnd_   = dim      
        MeshColumn%prop_%StencilStart_ = 0
        MeshColumn%prop_%StencilEnd_   = 0   
        if (.not. MeshColumn%prop_%has_dim(dim)) then
          MeshColumn%prop_%StencilStart_(dim) = -2 
          MeshColumn%prop_%StencilEnd_(dim)   = 2
        endif
        call CISL_integral(MeshColumn, columns, nx, ny, ntracer, dim, l0, InOrOut)
        ! integrate inside the current stack column with picewise parabolic method
        ax = MeshColumn%prop_%a_(dim, :)
        bx = MeshColumn%prop_%b_(dim, :)
        ! compute the range of the current lagrangian volume
        xLV(:) = MeshColumn%prop_%LV_(dim, i)%x_(:)
        ! scale the lagrangian volume into [-0.5, 0.5]*[-0.5, 0.5] as reference column
        if (dim == dimX) then
          xLV(:) = (xLV(:) - MeshColumn%prop_%x_(dimX))/MeshColumn%prop_%dx_(dimX) - 0.5_dp
        else
          xLV(:) = ( sin(xLV) - sin(MeshColumn%prop_%x_(dimY)) )/(sin(MeshColumn%prop_%x_(dimY) &
            + MeshColumn%prop_%dx_(dimY)) - sin(MeshColumn%prop_%x_(dimY))) - 0.5_dp
        endif
        ! update thMeshColumne value at arrival column
        ArrivalColumn%prop_%qmid_(:, dim, InOrOut) = ArrivalColumn%prop_%qmid_(:, dim, InOrOut) + &
                (ppm1D_remap(xLV(2), ax, bx, MeshColumn%prop_%qmid_(:, 3 - dim, InOrOut - 1), ntracer) - &
                  ppm1D_remap(xLV(1), ax, bx, MeshColumn%prop_%qmid_(:, 3 - dim, InOrOut - 1), ntracer))*dA
      endif
      MeshColumn%prop_%LV_(dim, i)%x_(:) = 0._dp
    enddo
  end subroutine CISL_stencil

  subroutine CISL_integral(MeshColumn, columns, nx, ny, ntracer, dim, l0, InOrOut)
    ! compute the modified cell integrated semi-Lagrangian scheme
    use class_parameter, only: dimX, dimY
    use class_AMRTypes,  only: ColumnList, ColumnStack, column, LV
    use class_column,    only: column_dA
    type(ColumnList),           intent(in)    :: columns
    type(column), pointer,      intent(inout) :: MeshColumn
    integer,                    intent(in)    :: nx, ny, ntracer, l0, InOrOut, dim
    type(ColumnStack)                         :: StackColumns

    type(column),      pointer                :: StackColumn

    logical                                   :: completed

    ! put the current column into a stack, the stack can ensure the necessary ghost cells are computed
    StackColumns%head_ => MeshColumn
    do while (associated(StackColumns%head_))
      StackColumn => StackColumns%head_
      ! get the coefficients of the current column in the stack
      call GetCoefficients(StackColumns, StackColumn, columns, InOrOut, dim, nx, ny, ntracer, l0, completed)
      ! if the stack column needs ghost cells, retrieve from the stack again
      if (.not. completed ) cycle
      ! if necessary, compute the necessary ghost cells in the column
      if ( (StackColumn%prop_%nGhost_ /= 0) .and. (.not. allocated(StackColumn%prop_%qghost_)) ) &
        call get_Ghost(StackColumn, ntracer, InOrOut, dim)
      ! in order to make sure the next time step; reinit the values of lagrangian volumes
      StackColumn%prop_%StencilStart_(dimX:dimY) = 0; StackColumn%prop_%StencilEnd_(dimX:dimY) = 0
      StackColumn%prop_%dimStart_ = 0;                StackColumn%prop_%dimEnd_ = 0

      StackColumns%head_ => StackColumn%nextLV_
      StackColumn%nextLV_ => null()
    enddo

  end subroutine CISL_integral

  subroutine GetCoefficients(StackColumns, MeshColumn, columns, InOrOut, dimLV, nx, ny, ntracer, l0, completed)
    ! get coefficients of the current column
    use class_kind,      only: dp
    use class_parameter, only: dimX, dimY
    use class_AMRTypes,  only: column, ColumnStack, ColumnList
    use class_column,    only: column_dA, column_get_leaves, ColumnList_get_column
    type(ColumnList),           intent(in)    :: columns
    type(ColumnStack),          intent(inout) :: StackColumns
    type(column),      pointer, intent(inout) :: MeshColumn
    integer,                    intent(in)    :: ntracer, dimLV, l0, nx, ny, InOrOut
    logical,                    intent(inout) :: completed
    type(column),      pointer                :: StencilColumn, LeafColumn
    type(ColumnStack)                         :: ColumnLeaves
    real(dp)                                  :: dx, dy, dA
    integer                                   :: idx(dimX:dimY), iCol(dimX:dimY)
    integer                                   :: dim, Lscale, ii, jj, iStencil
    integer                                   :: dimStart, dimEnd
    integer                                   :: StencilStart, StencilEnd

    completed = .false.
    dimStart = MeshColumn%prop_%dimStart_
    dimEnd = MeshColumn%prop_%dimEnd_

    do dim = dimStart, dimEnd
      ! print *, MeshColumn%prop_%i_, MeshColumn%l_, dim
      StencilStart = MeshColumn%prop_%StencilStart_(dim)
      StencilEnd   = MeshColumn%prop_%StencilEnd_(dim)
      do iStencil = StencilStart, StencilEnd
        if (iStencil == 0) then
          MeshColumn%prop_%qx_(dim, :, iStencil) = MeshColumn%prop_%qmid_(:, 3 - dimLV, InOrOut - 1)
          MeshColumn%prop_%qdx_(dim, iStencil) = MeshColumn%prop_%dx_(dim)
          cycle
        endif 
        iCol(3 - dim) = MeshColumn%prop_%i_(3 - dim)
        iCol(dim) = MeshColumn%prop_%i_(dim) + iStencil
        iCol = boundary_check(iCol, MeshColumn%l_, nx, ny, l0)
        ! get the stencil
        StencilColumn => ColumnList_get_column(columns, iCol(dimX), iCol(dimY), MeshColumn%l_, l0)
        if (StencilColumn%l_ < MeshColumn%l_) then
          ! if the Stencil Column has larger size than the original column
          ! we need ghost cells
          if ((StencilColumn%prop_%has_dim(dimX)) .and. (StencilColumn%prop_%has_dim(dimY))) then
            ! get the index of the ghost cells for q in MeshColumn
            if ((StencilColumn%prop_%nGhost_ <  MeshColumn%l_ - StencilColumn%l_) &
              .or. (.not. allocated(StencilColumn%prop_%qghost_)) )then
              StencilColumn%prop_%nGhost_  = MeshColumn%l_ - StencilColumn%l_
              if (allocated(StencilColumn%prop_%qghost_)) deallocate(StencilColumn%prop_%qghost_)
              call get_Ghost(StencilColumn, ntracer, InOrOut, dimLV)
            endif
            idx(dim)     = MeshColumn%prop_%i_(dim) + iStencil
            idx(3 - dim) = MeshColumn%prop_%i_(3 - dim)
            idx(:)       = boundary_check(idx, MeshColumn%l_, nx, ny, l0)
            ! scale the index to the same level of the ghost column and get required the index in the ghost column
            Lscale       = 2**(StencilColumn%l_ + StencilColumn%prop_%nGhost_ - MeshColumn%l_)
            idx(:)       = idx(:)*Lscale - StencilColumn%prop_%i_(:)*(2**StencilColumn%prop_%nGhost_)
            ! get the spatial interval of the ghost cell
            dy           = StencilColumn%prop_%dx_(dimY)/(2**StencilColumn%prop_%nGhost_)
            dx           = StencilColumn%prop_%dx_(dimX)/(2**StencilColumn%prop_%nGhost_)
            ! get the required value of the specific column
            MeshColumn%prop_%qx_(dim, :, iStencil) = 0._dp
            do ii = idx(dimX), idx(dimX) + Lscale - 1
              do jj = idx(dimY), idx(dimY) + Lscale - 1
                dA = dx*(sin(StencilColumn%prop_%x_(dimY) + (jj + 1)*dy) - sin(StencilColumn%prop_%x_(dimY) + jj*dy) )
                MeshColumn%prop_%qx_(dim, :, iStencil) = MeshColumn%prop_%qx_(dim, :, iStencil) &
                  + StencilColumn%prop_%qghost_(ii, jj, :)*dA
              enddo
            enddo
            MeshColumn%prop_%qdx_(dim, iStencil)    = StencilColumn%prop_%dx_(dim)/(2**(MeshColumn%l_ - StencilColumn%l_))
            dA = (StencilColumn%prop_%dx_(dimX)/(2**(MeshColumn%l_ - StencilColumn%l_)))*&
              (sin(StencilColumn%prop_%x_(dimY) + (idx(dimY) + Lscale)*dy) - &
              sin(StencilColumn%prop_%x_(dimY) + idx(dimY)*dy) )
            MeshColumn%prop_%qx_(dim, :, iStencil) = MeshColumn%prop_%qx_(dim, :, iStencil)/dA
          else
            ! if there is no ghost cells yet, put the stencil column to the stack
            ! and get the ghost column
            MeshColumn%prop_%StencilStart_(dim) = iStencil
            MeshColumn%prop_%dimStart_     = dim
            StencilColumn%prop_%nGhost_  = MeshColumn%l_ - StencilColumn%l_
            StencilColumn%prop_%dimStart_ = dimX
            StencilColumn%prop_%dimEnd_   = dimY
            StencilColumn%prop_%StencilStart_(dimX:) = -2
            StencilColumn%prop_%StencilEnd_(dimX:)   = 2
            ! test the 
            StencilColumn%nextLV_       => StackColumns%head_ 
            StackColumns%head_          => StencilColumn
            return
          endif
        else
          ! get all corresponding leaves of the stencil column
          ColumnLeaves = column_get_leaves(StencilColumn)
          dA = 0._dp; MeshColumn%prop_%qx_(dim, :, iStencil) = 0._dp
          Lscale = 2**(ColumnLeaves%head_%l_ - MeshColumn%l_)
          MeshColumn%prop_%qdx_(dim, iStencil) = ColumnLeaves%head_%prop_%dx_(dim)*Lscale
          do while (associated(ColumnLeaves%head_))
            LeafColumn => ColumnLeaves%head_
            ColumnLeaves%head_ => LeafColumn%nextLeaf_
            LeafColumn%nextLeaf_ => null()
            MeshColumn%prop_%qx_(dim, :, iStencil) = MeshColumn%prop_%qx_(dim, :, iStencil) &
              + LeafColumn%prop_%qmid_(:, 3 - dimLV, InOrOut - 1)*column_dA(LeafColumn)
            dA = dA + column_dA(LeafColumn)
          enddo
          MeshColumn%prop_%qx_(dim, :, iStencil) = MeshColumn%prop_%qx_(dim, :, iStencil)/dA
        endif
      enddo
      ! the dx in longitudinal direction must be scaled 
      if (dim == dimX) MeshColumn%prop_%qdx_(dim, :) = MeshColumn%prop_%qdx_(dim, :)/MeshColumn%prop_%cosp_
      ! compute the 1D coefficient from PPM interpolation
      if (StencilEnd == 2) then
        call PPM1D(MeshColumn%prop_%a_(dim, :), MeshColumn%prop_%b_(dim, :), &
          MeshColumn%prop_%qx_(dim, :, -2:), MeshColumn%prop_%qdx_(dim, -2:), -2, ntracer)
        MeshColumn%prop_%has_dim(dim) = .true.
      endif
      MeshColumn%prop_%qx_(dimX:, :, -2:) = 0._dp
      MeshColumn%prop_%qdx_(dim, :) = 0._dp
    enddo
    completed = .true.
  end subroutine GetCoefficients

  subroutine get_Ghost(MeshColumn, ntracer, InOrOut, dimLV)
  ! remap the coarse grid onto fine grid using 2-D PPM
    use class_kind,                only: dp
    use class_AMRTypes,            only: column
    use class_parameter,           only: dimX, dimY
    use class_column,              only: column_dA
    integer,               intent(in) :: ntracer, dimLV, InOrOut
    type(column), pointer, intent(in) :: MeshColumn
    real(dp)                          :: ax(ntracer + 1), bx(ntracer + 1), ay(ntracer + 1), by(ntracer + 1)
    integer                           :: i, j, nGhost
    real(dp)                          :: dx(dimX:dimY), x(dimX:dimY, 2), area, dA

    nGhost = MeshColumn%prop_%nGhost_
    allocate(MeshColumn%prop_%qghost_(0:2**nGhost - 1, 0:2**nGhost-1, ntracer + 1) )
    dx = MeshColumn%prop_%dx_/(2**nGhost)

    ax = MeshColumn%prop_%a_(dimX, :)
    bx = MeshColumn%prop_%b_(dimX, :)
    ay = MeshColumn%prop_%a_(dimY, :)
    by = MeshColumn%prop_%b_(dimY, :)
    area = column_dA(MeshColumn)

    ! remapping the refined cell
    do i = 0, ishft(1, nGhost) - 1
      do j = 0, ishft(1, nGhost) - 1
        x(:, 1) = MeshColumn%prop_%x_ + [i, j]*dx
        x(:, 2) = x(:, 1) + dx
        x(dimY, 1) = sin(x(dimY, 1))
        x(dimY, 2) = sin(x(dimY, 2))
        dA = dx(dimX)*(x(dimY, 2) - x(dimY, 1))
        x(dimX, :) = (x(dimX, :) - MeshColumn%prop_%x_(dimX))/MeshColumn%prop_%dx_(dimX) - 0.5_dp
        x(dimY, :) = (x(dimY, :) - sin(MeshColumn%prop_%x_(dimY)))/(sin(MeshColumn%prop_%x_(dimY) &
          + MeshColumn%prop_%dx_(dimY)) - sin(MeshColumn%prop_%x_(dimY))) - 0.5_dp
        MeshColumn%prop_%qghost_(i, j, :) = (PPM2D_remap(x(dimX, 2), x(dimY, :), ax, ay, &
            bx, by, MeshColumn%prop_%qmid_(:, 3 - dimLV, InOrOut - 1), ntracer) &
              - PPM2D_remap(x(dimX, 1), x(dimY, :), ax, ay, bx, by, MeshColumn%prop_%qmid_(:, 3-dimLV, InOrOut-1), ntracer))*area/dA       
      enddo
    enddo
  end subroutine get_Ghost

  subroutine PPM1D(ax, bx, q, dx, minStencil, ntracer)
  ! get 1D PPM coefficients
    use class_kind, only: dp
    integer, intent(in)          :: ntracer, minStencil
    real(dp), intent(in)         :: q(:, minStencil:)
    real(dp), intent(in)         :: dx(minStencil:)
    real(dp), intent(inout)        :: ax(1:), bx(1:)

    integer                      :: k, i
    real(dp)                     :: dm(-1:1), q_l, q_r

    do i = 1, ntracer + 1
      do k = -1, 1
        dm(k) = dm1(dx(-1+k:1+k), q(i, -1+k:1+k))
      enddo
      q_l = al1(dm(-1:0), dx(-2:1), q(i, -1:0))
      q_r = al1(dm(0:1), dx(-1:2), q(i, 0:1))

      ! call limiter(q(i, 0), dm(0), q_l, q_r)
      ax(i) = q_r - q_l
      bx(i) = 3.0_dp*(2.0_dp*q(i, 0) - (q_l + q_r))
    enddo
  end subroutine PPM1D

  function ppm2D_remap(x, y, ax, ay, bx, by, q0, ntracer)
  ! 2D PPM integration
    use class_kind, only: dp
    integer, intent(in)  :: ntracer
    real(dp), intent(in) :: ax(:), ay(:), bx(:), by(:), x, y(:), q0(:)
    real(dp)             :: dy, dy2, dy3
    real(dp)             :: ppm2D_remap(ntracer+1)
    dy = y(2) - y(1); dy2 = y(2)*y(2) - y(1)*y(1); dy3 = y(2)**3 - y(1)**3
    ppm2D_remap(:) = (q0 + (bx + by)/12._dp)*x*dy + 0.5*ax*dy*x*x + &
      0.5_dp*ay*x*dy2 - (bx/3._dp)*dy*x**3 - (by/3._dp)*x*dy3 
  end function ppm2D_remap

  function ppm1D_remap(x, ax, bx, q0, ntracer)
  ! 1D PPM integration
    use class_kind,     only: dp
    integer, intent(in)  :: ntracer
    real(dp), intent(in) :: ax(:), bx(:), x, q0(:)
    real(dp)             :: ppm1D_remap(ntracer+1)
    ppm1D_remap = (q0 + bx/12._dp)*x + 0.5_dp*ax*x*x - (bx/3._dp)*x*x*x
  end function ppm1D_remap

  function dm1(q_d, q_adv)
  ! slope calculation on unequidistant grid
    use class_kind, only: dp
    real(dp),intent(in)  :: q_d(-1:), q_adv(-1:)
    real(dp)             :: dm, dm1
    ! real(dp)             :: dm_max, dm_min
    dm = 0.5_dp*(q_d(0)/(q_d(-1) + q_d(0) + q_d(1))) * &
      (  ( (2*q_d(-1) + q_d(0)) / (q_d(1) + q_d(0)) )*(q_adv(1) - q_adv(0)) &
        + ( (2*q_d(1) + q_d(0)) / (q_d(-1) + q_d(0)) )*(q_adv(0) - q_adv(-1))  )
    dm1 = dm
    ! dm_max = max(q_adv(-1), q_adv(0), q_adv(1)) - q_adv(0)
    ! dm_min = q_adv(0) - min(q_adv(-1), q_adv(0), q_adv(1))
    ! dm1 = sign(min(abs(dm), dm_min, dm_max), dm)
  end function dm1

  function al1(dm, q_d, q_adv)
  ! quartic polynomial interpolation from ppn
    use class_kind, only: dp
    real(dp),intent(in)  :: dm(-1:), q_d(-2:), q_adv(-1:)
    real(dp)             :: al1
    al1 = q_adv(-1) &
      + (q_d(-1)/sum(q_d(-1:0)))*(q_adv(0) - q_adv(-1))&
      + (2._dp/sum(q_d(-2:1))) &
      *((q_d(0)*q_d(-1)/sum(q_d(-1:0))) &
        *( sum(q_d(-2:-1))/(2*q_d(-1) + q_d(0)) &
        - sum(q_d(0:1))/(2*q_d(0) + q_d(-1)) )*(q_adv(0) - q_adv(-1)) &
      - q_d(-1)*dm(0)*sum(q_d(-2:-1))/(2*q_d(-1) + q_d(0)) &
      + q_d(0)*dm(-1)*sum(q_d(0:1))/(q_d(-1) + 2*q_d(0)))
  end function al1

  function d2q(q_d, q_adv)
  ! gradient of slope
    use class_kind, only: dp
    real(dp),intent(in)  :: q_d(-1:1), q_adv(-1:1)
    real(dp)             :: d2q
    d2q = (1._dp/(q_d(-1) + q_d(0) + q_d(1))) *((q_adv(1) - &
      q_adv(0))/(q_d(1) + q_d(0)) - (q_adv(0) - &
        q_adv(-1))/(q_d(0) + q_d(-1)))
  end function d2q

  ! subroutine limiter(q0, da, q_l, q_r)
  ! ! L04 limiter appendix B3 AND B4
  !   use class_kind, only: dp
  !   real(dp),intent(in)     :: q0, da
  !   real(dp), intent(inout) :: q_l, q_r
  !   q_l = q0 - sign(min(abs(2*da), abs(q_l - q0)), 2*da)
  !   q_r = q0 + sign(min(abs(2*da), abs(q_r - q0)), 2*da)
  ! end subroutine limiter

  function boundary_check(i, l, nx, ny, l0) result(val)
  ! shifting the index based on boundary condition
    use class_parameter,            only: dimX, dimY
    integer, intent(in)                :: i(:), nx, ny
    integer, intent(in)                :: l
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

  function get_range(i, l, x, l0) result(val)
  ! get the coordinate value based on index and level
    use class_kind, only: dp
    integer,  intent(in) :: i, l, l0
    real(dp), intent(in) :: x(0:)
    integer              :: Lscale, i0
    real(dp)             :: dx
    real(dp)             :: val(2)
    Lscale = 2**(l - l0)
    i0 = i/Lscale
    dx = (x(i0 + 1) - x(i0))/Lscale
    val(1) = x(i0) + (i - i0*Lscale)*dx
    val(2) = x(i0) + (i + 1 - i0*Lscale)*dx
  end function get_range
end module class_CISL