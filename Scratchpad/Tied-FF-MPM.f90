PROGRAM Tied_Free_Field_MPM
  USE JSON_MODULE
  USE MPM_CORE
  USE CLASS_MESH
  USE LOGGING
  USE IO
  IMPLICIT NONE
  !************************** VARIABLES DECLARATION ***************************!
  INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
  REAL(iwp)          :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp, three=3.0_iwp
  CHARACTER(128)     :: input_directory, input_filename
  CHARACTER(:), ALLOCATABLE :: output_directory

  ! Data Structure
  TYPE(json_file) :: input_json ! Main input file
  TYPE(mpm_body)  :: mbod,fbod(2) ! Simulated body and free-field bodies
  TYPE(mesh)      :: grid       ! Computaional background mesh

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

  ! Geometry
  INTEGER :: nx,ny    ! Number of cell in each directions
  REAL(iwp):: w1,h1   ! Model width and height 

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
  ! Define main body geometry
  nx=40; ny=20
  w1=20.0; h1=10.0

  !--- Load Material Points and Free-Field Boundaries
  ALLOCATE(mbod%particles,fbod(1)%particles,fbod(2)%particles)
  CALL mbod%particles%GENERATE(nx,ny,w1,h1,node=nip,                     &
    offsetx=5,offsety=2,nst=nst,ndim=ndim,nodof=nodof)
  CALL fbod(1)%particles%GENERATE(nx=2,ny=ny,w1=1.0_iwp,h1=h1,node=nip,  &
    offsetx=0,offsety=2,nst=nst,ndim=ndim,nodof=nodof)
  CALL fbod(2)%particles%GENERATE(nx=2,ny=ny,w1=1.0_iwp,h1=h1,node=nip,  &
    offsetx=nx+5+3,offsety=2,nst=nst,ndim=ndim,nodof=nodof)

  !--- Load/Generate MPM Mesh
  CALL grid%GENERATE(nx=nx+4,ny=ny+5,w1=22.0_iwp,h1=12.5_iwp,node=nip,offsetx=3,&
    nst=nst,ndim=ndim,nodof=nodof)


  !--- Insert Material Points to Background Mesh
  CALL mbod%SET_MESH(grid)
  
  !--- Load Material Constitutive Model Parameters
  ! CALL mbod%SET_MATERIAL(input_json)


  !**************************** INITIAL CONDITIONS ****************************!
  !                                                                            !
  !--- Determine Global Boundary Conditions
  
  !--- Activate/Flag Initial Compulational Grid
  CALL mbod%LOCATE_PARTICLES()
  CALL mbod%GET_PARTICLE_SUPPORT_SIZE()
  CALL mbod%GET_SUPPORT_ELEMENTS()
  CALL mbod%GET_SUPPORT_NODES()

  !--- Calculate Initial Material Points State (Mass, Velocity, Stresses, Loads)
  
  !--- Calculate Initial Conditions and Properties of Material Points
  ! Compute Dee Matrix (Constitutive Model)
  
  !--- Print Initial Conditions Visualization (VTK)

  CALL PLOT(0)
    
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
    CALL mbod%particles%RESET_MP()

    !--- Compute Nodal Mass
    
    !--- Compute Domain Stiffness Matrix
    ! CALL mbod%CONSTRUCT_DMMPM_STIFFNESS(nip, nst)

    !--- Compute Domain Mass Matrix

    !--- Combined Stiffness
    
    ! FOR TESTING
    ! mbod%gm_coord(1,:) = mbod%gm_coord(1,:) + 0.0002

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
    IF (print_steps == output_steps .and. .True.) THEN
      CALL PLOT(step)
      print_steps = 0
    END IF
    !--- Reallocate Memory
    !--- Determine and Activate New Element for Next Time Step
    ! CALL mbod%FLAG_ELEMENTS()
    ! CALL mbod%GET_SUPPORT_ELEMENTS()
  END DO TIME_STEPS

  CONTAINS

  SUBROUTINE PLOT(input)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::input
    CALL IO_POINT_VIZ(                         &
      input=input,                             &
      coord=mbod%particles%gm_coord,                     &
      a_ins=mbod%particles%a_ins,                        &
      evpt=mbod%particles%epsinvacum,                    &
      m_stress=mbod%particles%m_stress,                  &
      m_stress_inc=mbod%particles%m_stress_change,       &
      acc=mbod%particles%m_acc,                          &
      velocity=mbod%particles%m_velocity,                &
      cohesion=mbod%particles%mpcp,                      &
      devstress=mbod%particles%devstress,                &
      meanstress=mbod%particles%mean_stress,             &
      mpyield=mbod%particles%mpyield,                    &
      directory=output_directory,              &
      argv="MPM_Particles"                     &
    )
    CALL IO_POINT_VIZ(                         &
      input=input,                             &
      coord=fbod(1)%particles%gm_coord,                  &
      a_ins=fbod(1)%particles%a_ins,                     &
      evpt=fbod(1)%particles%epsinvacum,                 &
      m_stress=fbod(1)%particles%m_stress,               &
      m_stress_inc=fbod(1)%particles%m_stress_change,    &
      acc=fbod(1)%particles%m_acc,                       &
      velocity=fbod(1)%particles%m_velocity,             &
      cohesion=fbod(1)%particles%mpcp,                   &
      devstress=fbod(1)%particles%devstress,             &
      meanstress=fbod(1)%particles%mean_stress,          &
      mpyield=fbod(1)%particles%mpyield,                 &
      directory=output_directory,              &
      argv="FEM_1_particles"                   &
    )
    CALL IO_POINT_VIZ(                         &
      input=input,                             &
      coord=fbod(2)%particles%gm_coord,                  &
      a_ins=fbod(2)%particles%a_ins,                     &
      evpt=fbod(2)%particles%epsinvacum,                 &
      m_stress=fbod(2)%particles%m_stress,               &
      m_stress_inc=fbod(2)%particles%m_stress_change,    &
      acc=fbod(2)%particles%m_acc,                       &
      velocity=fbod(2)%particles%m_velocity,             &
      cohesion=fbod(2)%particles%mpcp,                   &
      devstress=fbod(2)%particles%devstress,             &
      meanstress=fbod(2)%particles%mean_stress,          &
      mpyield=fbod(2)%particles%mpyield,                 &
      directory=output_directory,              &
      argv="FEM_2_particles"                   &
    )
    CALL IO_PARAVIEW(                          &
      input=input,                             &
      node_type=4,                             &
      coord=grid%g_coord,                      &
      num=grid%g_num,                          &
      directory=output_directory,              &
      argv="MPM_mesh"                          &
    )
    CALL IO_PARAVIEW(                          &
      input=input,                             &
      node_type=4,                             &
      coord=fbod(1)%particles%g_coord,                   &
      num=fbod(1)%particles%g_num,                       &
      directory=output_directory,              &
      argv="FEM_1_mesh"                        &
    )
    CALL IO_PARAVIEW(                          &
      input=input,                             &
      node_type=4,                             &
      coord=fbod(2)%particles%g_coord,                   &
      num=fbod(2)%particles%g_num,                       &
      directory=output_directory,              &
      argv="FEM_2_mesh"                        &
    )
  END SUBROUTINE
END PROGRAM Tied_Free_Field_MPM