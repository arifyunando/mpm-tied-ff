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
  TYPE(json_file) :: input_json ! Main input file
  TYPE(mpm_body)  :: mbod       ! Simulated body
  TYPE(fem_body)  :: ffbod(2)   ! Free-Field bodies
  TYPE(mpm_grid)  :: mesh       ! Computaional background mesh

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
  INTEGER :: ndim=2   ! Number of dimensions (x, y, z)
  INTEGER :: nst=4    ! Number of Stress/Strain Terms 
                      ! (3: Plane Stress/Strain, 4: Axisymmetric, 6: 3-D)
                      ! With appropriate Linear Derivative Operator (Bee),
                      ! it is also possible to have Plane Stress/Strain
                      ! conditions with nst = 4
  INTEGER :: nodof=2  ! Number of Degree of Freedom (ndim + 1 if there's water)
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
  CALL mbod%LOAD_PARTICLES(1, input_directory, input_json)

  !--- Initiate Mesh
  CALL mesh%LOAD_MESH(input_directory, input_json)
  
  !--- Insert Material Points to Background Mesh
  CALL mbod%SET_MESH(mesh)
  
  !--- Load Material Constitutive Model Parameters
  CALL mbod%SET_MATERIAL(input_json)


  !**************************** INITIAL CONDITIONS ****************************!
  !                                                                            !
  
  !--- Activate/Flag Initial Compulational Grid
  CALL mbod%FLAG_ELEMENTS()
  CALL mbod%GET_PARTICLE_SUPPORT_SIZE()
  CALL mbod%GET_SUPPORT_ELEMENTS()

  !--- Calculate Initial Material Points State (Mass, Velocity, Stresses, Loads)
  
  !--- Calculate Initial Conditions and Properties of Material Points
  ! Compute Dee Matrix (Constitutive Model)
  
  !--- Print Initial Conditions Visualization (VTK)

  CALL IO_POINT_VIZ(                         &
    input=0,                                 &
    coord=mbod%gm_coord,                     &
    a_ins=mbod%a_ins,                        &
    evpt=mbod%epsinvacum,                    &
    m_stress=mbod%m_stress,                  &
    m_stress_inc=mbod%m_stress_change,       &
    acc=mbod%m_acc,                          &
    velocity=mbod%m_velocity,                &
    cohesion=mbod%mpcp,                      &
    devstress=mbod%devstress,                &
    meanstress=mbod%mean_stress,             &
    mpyield=mbod%mpyield,                    &
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
    CALL mbod%RESET_MP()
    
    !--- Compute Domain Stiffness Matrix
    CALL mbod%CONSTRUCT_DMMPM_STIFFNESS(nip, nst)

    !--- Compute Domain Mass Matrix

    !--- Combined Stiffness
    
    ! FOR TESTING
    mbod%gm_coord(1,:) = mbod%gm_coord(1,:) + 0.0002

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
        coord=mbod%gm_coord,                     &
        a_ins=mbod%a_ins,                        &
        evpt=mbod%epsinvacum,                    &
        m_stress=mbod%m_stress,                  &
        m_stress_inc=mbod%m_stress_change,       &
        acc=mbod%m_acc,                          &
        velocity=mbod%m_velocity,                &
        cohesion=mbod%mpcp,                      &
        devstress=mbod%devstress,                &
        meanstress=mbod%mean_stress,             &
        mpyield=mbod%mpyield,                    &
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
    CALL mbod%FLAG_ELEMENTS()
    CALL mbod%GET_SUPPORT_ELEMENTS()
  END DO TIME_STEPS
END PROGRAM Tied_Free_Field_MPM