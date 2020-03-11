module class_refine
implicit none
private
public :: mesh_refine, mesh_refineInit
contains
  subroutine mesh_refineInit(self, nx, ny, nz, ntracer, MaxRef, theta_r, theta_c, dt, uInput, vInput, nLimit)
  ! initial refinement of the mesh: no coarsening and predictive step
    use class_mesh,                only: mesh
    use class_kind,                only: dp
    use class_AMRTypes,            only: column, face
    use class_column,              only: ColumnList_AMR_Mark, ColumnList_refine, column_stability
    use class_column,              only: ColumnList_q
    use class_setup,               only: setup_init
    use class_parameter,           only: dimX, dimY
    use class_setupWind,           only: setup_interp    
    type(mesh),            intent(inout) :: self(1:)
    real(dp),              intent(in)    :: dt, theta_r, theta_c
    real(dp),              intent(in)    :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)
    integer,               intent(in)    :: nx, ny, nz, ntracer, MaxRef, nLimit
    type(column), pointer                :: MeshColumn, MeshColumnNext
    type(face),   pointer                :: MeshFace
    integer                              :: i, nchange, nCycle, nRefine, dim, Lscale
    nRefine = 1; nCycle = 0
    ! update velocity
    do dim = dimX, dimY
      MeshFace => self(2)%faces_(dim)%head_
      do i = 1, self(2)%faces_(dim)%n_
        call setup_interp(MeshFace, self(2)%coord_%x_, self(2)%coord_%y_, &
          self(2)%coord_%l0_, nx, ny, nz, uInput, vInput, dim)
        MeshFace => MeshFace%next_
      enddo
    enddo    
    do while ((nRefine /= 0) .and. (nCycle < nLimit))
      nchange = self(2)%columns_%n_
      nRefine = nchange
      ! mark the number of refinement
      call ColumnList_AMR_Mark(self(2)%columns_, self(2)%faces_, self(2)%coord_%l0_, &
        MaxRef, ntracer, nx, ny, nz, theta_r, theta_c, dt, self(2)%coord_%x_, &
        self(2)%coord_%y_, self(2)%coord_%dxMax_, uInput, vInput, .true., nCycle)
      ! refine the meshes
      MeshColumn => self(2)%columns_%head_
      do while(associated(MeshColumn))
        MeshColumnNext => MeshColumn%next_
        if (MeshColumn%prop_%i_(dimY) == 0) MeshColumn%prop_%NumChange_ = 0
        if (MeshColumn%prop_%i_(dimY) == ny - 1) MeshColumn%prop_%NumChange_ = 0
        Lscale = 2**(MeshColumn%l_ - self(2)%coord_%l0_)
        if ( (MeshColumn%prop_%i_(dimY)/Lscale == 1) .and. &
          (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
          MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_
        if ( (MeshColumn%prop_%i_(dimY)/Lscale == ny - 2) .and. &
          (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
          MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_
        if (MeshColumn%prop_%NumChange_ > 0) &
          call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
            self(2)%coord_%x_, self(2)%coord_%y_, self(2)%coord_%l0_, nx, ny, nz, uInput, vInput)
        MeshColumn => MeshColumnNext
      enddo
      ! ensure there is no stability problem due to the mesh refinement
      nchange = self(2)%columns_%n_ - nchange
      do while (nchange /= 0)
        MeshColumn => self(2)%columns_%head_
        do i = 1, self(2)%columns_%n_
          call column_stability(MeshColumn, dt)
          MeshColumn => MeshColumn%next_
        enddo
        nchange = self(2)%columns_%n_
        MeshColumn => self(2)%columns_%head_
        do while(associated(MeshColumn))
          MeshColumnNext => MeshColumn%next_
          if (MeshColumn%prop_%i_(dimY) == 0) MeshColumn%prop_%NumChange_ = 0
          if (MeshColumn%prop_%i_(dimY) == ny - 1) MeshColumn%prop_%NumChange_ = 0
          Lscale = 2**(MeshColumn%l_ - self(2)%coord_%l0_)
          if ( (MeshColumn%prop_%i_(dimY)/Lscale == 1) .and. &
            (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
            MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_
          if ( (MeshColumn%prop_%i_(dimY)/Lscale == ny - 2) .and. &
            (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
            MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_            
          if (MeshColumn%prop_%NumChange_ > 0) &
            call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
              self(2)%coord_%x_, self(2)%coord_%y_, self(2)%coord_%l0_, nx, ny, nz, uInput, vInput)
          MeshColumn => MeshColumnNext
        enddo
        nchange = self(2)%columns_%n_ - nchange
      enddo
      nCycle = nCycle + 1
      ! reinitialize the tracer concentration
      call ColumnList_q(self(2)%columns_)
      ! call setup_init(self(2)%columns_, nz, trim(adjustl(TestName)))
      nRefine = self(2)%columns_%n_ - nRefine
    enddo
    call ColumnList_q(self(2)%columns_)
    MeshColumn => self(2)%columns_%head_
    do i = 1, self(2)%columns_%n_
      MeshColumn%prop_%q_(1:, 1:) = MeshColumn%prop_%q0_(1:, 1:)
      MeshColumn => MeshColumn%next_
    enddo    
  end subroutine mesh_refineInit

  subroutine mesh_refine(self, nx, ny, nz, ntracer, MaxRef, theta_r, theta_c, dt, uInput, vInput, nLimit)
  ! refine the mesh at each time step
    use class_mesh,                only: mesh
    use class_kind,                only: dp
    use class_AMRTypes,            only: column, face
    use class_column,              only: ColumnList_AMR_Mark, ColumnList_refine, column_stability
    use class_column,              only: ColumnList_q
    use class_predict,             only: predict
    use class_parameter,           only: dimX, dimY
    use class_setupWind,           only: setup_interp
    type(mesh),            intent(inout) :: self(1:)
    real(dp),              intent(in)    :: dt, theta_r, theta_c
    real(dp),              intent(in)    :: uInput(1:, 1:, 1:), vInput(1:, 0:, 1:)    
    integer,               intent(in)    :: nx, ny, nz, ntracer, MaxRef, nLimit
    type(column), pointer                :: MeshColumn, MeshColumnNext
    type(face),   pointer                :: MeshFace
    integer                              :: i, nchange, nCycle, nRefine, dim, Lscale

    ! update velocity
    do dim = dimX, dimY
      MeshFace => self(2)%faces_(dim)%head_
      do i = 1, self(2)%faces_(dim)%n_
        call setup_interp(MeshFace, self(2)%coord_%x_, self(2)%coord_%y_, &
          self(2)%coord_%l0_, nx, ny, nz, uInput, vInput, dim)
        MeshFace => MeshFace%next_
      enddo
    enddo
    nRefine = 1; nCycle = 0
    do while ((nRefine /= 0) .and. (nCycle < nLimit))
      nchange = self(2)%columns_%n_
      nRefine = nchange
      ! predict new time step
      call predict(self, self(2)%coord_%l0_, nx, ny, nz, ntracer, dt)
      ! mark the refined meshes
      call ColumnList_AMR_Mark(self(2)%columns_, self(2)%faces_, self(2)%coord_%l0_, &
        MaxRef, ntracer, nx, ny, nz, theta_r, theta_c, dt, self(2)%coord_%x_, &
        self(2)%coord_%y_, self(2)%coord_%dxMax_, uInput, vInput, .false., nCycle)
      ! refine the mesh
      MeshColumn => self(2)%columns_%head_
      do while(associated(MeshColumn))
        MeshColumnNext => MeshColumn%next_
          if (MeshColumn%prop_%i_(dimY) == 0) MeshColumn%prop_%NumChange_ = 0
          if (MeshColumn%prop_%i_(dimY) == ny - 1) MeshColumn%prop_%NumChange_ = 0
          Lscale = 2**(MeshColumn%l_ - self(2)%coord_%l0_)
          if ( (MeshColumn%prop_%i_(dimY)/Lscale == 1) .and. &
            (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
            MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_
          if ( (MeshColumn%prop_%i_(dimY)/Lscale == ny - 2) .and. &
            (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
            MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_            
        if (MeshColumn%prop_%NumChange_ > 0) &
          call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
            self(2)%coord_%x_, self(2)%coord_%y_, self(2)%coord_%l0_, nx, ny, nz, uInput, vInput)
        MeshColumn => MeshColumnNext
      enddo
      ! ensure the stability of the computation for next step
      nchange = self(2)%columns_%n_ - nchange
      do while (nchange /= 0)
        MeshColumn => self(2)%columns_%head_
        do i = 1, self(2)%columns_%n_
          call column_stability(MeshColumn, dt)
          MeshColumn => MeshColumn%next_
        enddo
        nchange = self(2)%columns_%n_
        MeshColumn => self(2)%columns_%head_
        do while(associated(MeshColumn))
          MeshColumnNext => MeshColumn%next_
          if (MeshColumn%prop_%i_(dimY) == 0) MeshColumn%prop_%NumChange_ = 0
          if (MeshColumn%prop_%i_(dimY) == ny - 1) MeshColumn%prop_%NumChange_ = 0
          Lscale = 2**(MeshColumn%l_ - self(2)%coord_%l0_)
          if ( (MeshColumn%prop_%i_(dimY)/Lscale == 1) .and. &
            (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
            MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_
          if ( (MeshColumn%prop_%i_(dimY)/Lscale == ny - 2) .and. &
            (MeshColumn%l_ + MeshColumn%prop_%NumChange_ > self(2)%coord_%l0_ + 1) ) &
            MeshColumn%prop_%NumChange_ = self(2)%coord_%l0_ + 1 - MeshColumn%l_            
          if (MeshColumn%prop_%NumChange_ > 0) &
            call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
              self(2)%coord_%x_, self(2)%coord_%y_, self(2)%coord_%l0_, nx, ny, nz, uInput, vInput)
          MeshColumn => MeshColumnNext
        enddo
        nchange = self(2)%columns_%n_ - nchange
      enddo
      nCycle = nCycle + 1
      ! update the columns for next predictive step
      call ColumnList_q(self(2)%columns_)
      nRefine = self(2)%columns_%n_ - nRefine
    enddo
    call ColumnList_q(self(2)%columns_)
  end subroutine mesh_refine
end module class_refine
