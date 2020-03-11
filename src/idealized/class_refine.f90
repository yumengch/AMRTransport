module class_refine
implicit none
private
public :: mesh_refine, mesh_refineInit
contains
  subroutine mesh_refineInit(self, nx, ny, nz, t, ntracer, MaxRef, theta_r, theta_c, dt, beta, TestName, nLimit)
  ! initial refinement of the mesh: no coarsening and predictive step
    use class_mesh,                only: mesh
    use class_kind,                only: dp
    use class_AMRTypes,            only: column, face
    use class_column,              only: ColumnList_AMR_Mark, ColumnList_refine, column_stability
    use class_column,              only: ColumnList_q
    use class_setup,               only: setup_init
    use class_parameter,           only: dimX, dimY
    use class_setupWind,           only: setup_ideal    
    type(mesh),            intent(inout) :: self(1:)
    character(len =*),     intent(in)    :: TestName
    real(dp),              intent(in)    :: dt, theta_r, theta_c, beta
    integer,               intent(in)    :: nx, ny, nz, t, ntracer, MaxRef, nLimit
    type(column), pointer                :: MeshColumn, MeshColumnNext
    type(face),   pointer                :: MeshFace
    integer                              :: i, nchange, nCycle, nRefine, dim
    nRefine = 1; nCycle = 0
    ! update velocity
    do dim = dimX, dimY
      MeshFace => self(2)%faces_(dim)%head_
      do i = 1, self(2)%faces_(dim)%n_
        call setup_ideal(MeshFace, TestName, nz, t, dt, beta, dim)
        MeshFace => MeshFace%next_
      enddo
    enddo

    do while ((nRefine /= 0) .and. (nCycle < nLimit))
      nchange = self(2)%columns_%n_
      nRefine = nchange
      ! mark the number of refinement
      call ColumnList_AMR_Mark(self(2)%columns_, self(2)%faces_, self(2)%coord_%l0_, &
        MaxRef, ntracer, nx, ny, nz, theta_r, theta_c, dt, self(2)%coord_%x_, &
        self(2)%coord_%y_, self(2)%coord_%dxMax_, beta, t, "Old", .true., nCycle)
      ! refine the meshes
      MeshColumn => self(2)%columns_%head_
      do while(associated(MeshColumn))
        MeshColumnNext => MeshColumn%next_
        if (MeshColumn%prop_%NumChange_ > 0) &
          call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
            nz, TestName, dt, beta, t)
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
          if (MeshColumn%prop_%NumChange_ > 0) &
            call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
              nz, TestName, dt, beta, t)
          MeshColumn => MeshColumnNext
        enddo
        nchange = self(2)%columns_%n_ - nchange
      enddo
      nCycle = nCycle + 1
      ! reinitialize the tracer concentration
      call setup_init(self(2)%columns_, nz, trim(adjustl(TestName)))
      nRefine = self(2)%columns_%n_ - nRefine
    enddo
  end subroutine mesh_refineInit

  subroutine mesh_refine(self, nx, ny, nz, t, ntracer, MaxRef, theta_r, theta_c, dt, beta, TestName, nLimit)
  ! refine the mesh at each time step
    use class_mesh,                only: mesh
    use class_kind,                only: dp
    use class_AMRTypes,            only: column, face
    use class_column,              only: ColumnList_AMR_Mark, ColumnList_refine, column_stability
    use class_column,              only: ColumnList_q
    use class_predict,             only: predict
    use class_parameter,           only: dimX, dimY
    use class_setupWind,           only: setup_ideal
    type(mesh),            intent(inout) :: self(1:)
    character(len =*),     intent(in)    :: TestName
    real(dp),              intent(in)    :: dt, theta_r, theta_c, beta
    integer,               intent(in)    :: nx, ny, nz, t, ntracer, MaxRef, nLimit
    type(column), pointer                :: MeshColumn, MeshColumnNext
    type(face),   pointer                :: MeshFace
    integer                              :: i, nchange, nCycle, nRefine, dim

    ! update velocity
    do dim = dimX, dimY
      MeshFace => self(2)%faces_(dim)%head_
      do i = 1, self(2)%faces_(dim)%n_
        call setup_ideal(MeshFace, TestName, nz, t, dt, beta, dim)
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
        self(2)%coord_%y_, self(2)%coord_%dxMax_, beta, t, TestName, .false., nCycle)
      ! refine the mesh
      MeshColumn => self(2)%columns_%head_
      do while(associated(MeshColumn))
        MeshColumnNext => MeshColumn%next_

        if (MeshColumn%prop_%NumChange_ > 0) &
          call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
            nz, TestName, dt, beta, t)
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
          if (MeshColumn%prop_%NumChange_ > 0) &
            call ColumnList_refine(self(2)%columns_, MeshColumn, self(2)%faces_, &
              nz, TestName, dt, beta, t)
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
