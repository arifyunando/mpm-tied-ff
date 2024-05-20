PROGRAM Tied_Free_Field_MPM
  USE MPM_CORE
  USE FEM_CORE
  USE IO
  USE LOGGING
  USE JSON_MODULE
  IMPLICIT NONE
  !************************** VARIABLES DECLARATION ***************************!
  INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
  REAL(iwp)          :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp, three=3.0_iwp
  CHARACTER(128)     :: input_directory, input_filename
  CHARACTER(:), ALLOCATABLE :: output_directory

  ! Data Structure
  TYPE(json_file) :: input_json           ! Main input file
  TYPE(mpm_body), ALLOCATABLE  :: mbod(:) ! may have multiple mp bodies
  TYPE(fem_body), ALLOCATABLE  :: fbod(:) ! may have multiple fe bodies
  TYPE(mpm_grid)  :: mesh                 ! main computaion mesh, singleton

  !--- Time Stepping and Iteration Variables
  ! Time Stepping
  INTEGER   :: step, nsteps, output_steps, print_steps
  REAL(iwp) :: dt
  ! Plastic iteration
  INTEGER   :: iteration, iteration_limit
  LOGICAL   :: converged, found

  ! Dummies
  INTEGER :: i,j,k,index,dummy
  CHARACTER(1) :: i_char

  !**************************** CONSTANT PARAMETER ****************************!
  !                                                                            !
  INTEGER :: nst=4    ! Number of Stress/Strain Terms 
                      ! (3: Plane Stress/Strain, 4: Axisymmetric, 6: 3-D)
                      ! With appropriate Linear Derivative Operator (Bee),
                      ! it is also possible to have Plane Stress/Strain
                      ! conditions with nst = 4
  INTEGER :: nodof=2  ! Number of Degree of Freedom
  INTEGER :: nip=4    ! Number of Gaussian Integration Point for Double Mapping

  !**************************** PREPARATION PHASE *****************************!
  !                                                                            !
  
  !--- Initiate Program Settings  
  ! Load Input JSON and Get Input/Output Path
  CALL IO_GET_INPUT_DIRECTORY(input_directory, input_filename)
  CALL IO_LOAD_JSON(input_directory, input_filename, input_json)
  CALL input_json%GET('post_processing.path', output_directory)
 
  ! Analysis and Post-Processing Parameters
  CALL input_json%GET('post_processing.output_steps', output_steps)
  CALL input_json%GET('analysis.nsteps', nsteps)
  CALL input_json%GET('analysis.dt', dt)
  CALL input_json%GET('analysis.max_iter', iteration_limit)

  !***************************** MODEL INITIATION *****************************!
  !                                                                            !

  !--- Load Material Points / Particle Bodies
  CALL input_json%INFO('particles', found, dummy, index)
  ALLOCATE(mbod(index))
  IF ( found ) THEN
    DO i=1, index
      CALL mbod(i)%LOAD_PARTICLES(i, input_directory, input_json)
    END DO
  END IF

  !--- Initiate Mesh
  CALL mesh%LOAD_MESH(input_directory, input_json)
  mesh%mesh_ind=4
  
  !--- Insert Material Points to Background Mesh
  CALL mbod(1)%SET_MESH(mesh)
  
  !--- Load Material Constitutive Model Parameters
  CALL mbod(1)%SET_MATERIAL(input_json)


  !**************************** INITIAL CONDITIONS ****************************!
  !                                                                            !
  
  !--- Activate/Flag Initial Compulational Grid
  CALL mbod(1)%FLAG_ELEMENTS
  CALL mbod(1)%ACTIVATE_ELEMENTS

  !--- Calculate Initial Material Points State (Mass, Velocity, Stresses, Loads)
  
  !--- Calculate Initial Conditions and Properties of Material Points
  ! Compute Dee Matrix (Constitutive Model)
  
  !--- Print Initial Conditions Visualization (VTK)

  CALL IO_POINT_VIZ(                         &
    input=0,                                 &
    coord=mbod(1)%gm_coord,                  &
    a_ins=mbod(1)%a_ins,                     &
    evpt=mbod(1)%epsinvacum,                 &
    m_stress=mbod(1)%m_stress,               &
    m_stress_inc=mbod(1)%m_stress_change,    &
    acc=mbod(1)%m_acc,                       &
    velocity=mbod(1)%m_velocity,             &
    cohesion=mbod(1)%mpcp,                   &
    devstress=mbod(1)%devstress,             &
    meanstress=mbod(1)%mean_stress,          &
    mpyield=mbod(1)%mpyield,                 &
    directory=output_directory               &
  )
  
  CALL IO_PARAVIEW(                          &
    input=0,                                 &
    node_type=4,                             &
    coord=mesh%g_coord,                      &
    num=mesh%num,                            &
    directory=output_directory               &
  )
    
  !****************************** TIME STEPPING *******************************!
  !                                                                            !

  !--- Initiate Time Stepping Variables
  print_steps = 0
  TIME_STEPS : DO step = 1, nsteps

    !**************************** Pre-Calculations ****************************!
    PRINT '(A, I8, A, I8)', "MPM Calculation Step :", step, "/", nsteps

    !--- Body Solution
    ! Reset Variables
    converged = .TRUE.
    CALL mbod(1)%RESET_MP()
    
    !--- Compute Domain Stiffness Matrix
    CALL mbod(1)%CONSTRUCT_DMMPM_STIFFNESS(nip, nst)

    ! Compute Mesh Momentum and Mass (GIMP / CMPM)

    !--- Compute Domain Mass Matrix

    !--- Combined Stiffness

    !************************ Plastic Strain Iteration ************************!
    iteration = 0
    PLASTIC_ITERATION : DO WHILE (iteration < iteration_limit)
      iteration = iteration + 1
      ! PRINT '(A, I8)', "  Iteration", iteration
      EXIT
    END DO PLASTIC_ITERATION
    
    !*************************** Post-Calculations ****************************!
    !--- Reverse Mapping (Node -> MP)
    !--- Save Material Point Data (VTK)
    print_steps = print_steps + 1 
    IF (print_steps == output_steps .and. .false.) THEN
      CALL IO_POINT_VIZ(                         &
        input=step,                              &
        coord=mbod(1)%gm_coord,                  &
        a_ins=mbod(1)%a_ins,                     &
        evpt=mbod(1)%epsinvacum,                 &
        m_stress=mbod(1)%m_stress,               &
        m_stress_inc=mbod(1)%m_stress_change,    &
        acc=mbod(1)%m_acc,                       &
        velocity=mbod(1)%m_velocity,             &
        cohesion=mbod(1)%mpcp,                   &
        devstress=mbod(1)%devstress,             &
        meanstress=mbod(1)%mean_stress,          &
        mpyield=mbod(1)%mpyield,                 &
        directory=output_directory               &
      )
      
      CALL IO_PARAVIEW(                          &
        input=step,                              &
        node_type=4,                             &
        coord=mesh%g_coord,                      &
        num=mesh%num,                            &
        directory=output_directory               &
      )

      print_steps = 0
    END IF
    !--- Reallocate Memory
    !--- Determine and Activate New Element for Next Time Step
  END DO TIME_STEPS
END PROGRAM Tied_Free_Field_MPM