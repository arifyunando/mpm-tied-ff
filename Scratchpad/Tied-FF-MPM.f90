PROGRAM Tied_Free_Field_MPM
  USE MPM_Core
  USE FEM_Core
  USE MPM_IO
  USE LOGGING
  USE JSON_MODULE
  IMPLICIT NONE
  !************************** VARIABLES DECLARATION ***************************!
  INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
  REAL(iwp)          :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp, three=3.0_iwp
  CHARACTER(128)     :: input_directory, input_filename
  CHARACTER(:), ALLOCATABLE :: output_directory
  CHARACTER(:), ALLOCATABLE :: mesh_filename, particles_filename, entity_filename

  ! Data Structure
  TYPE(mpm_body)  :: mbod
  TYPE(fem_body)  :: fbod(2)
  TYPE(mpm_grid)  :: mesh
  TYPE(json_file) :: input_json, entity_json

  !--- Time Stepping and Iteration Variables
  ! Time Stepping
  INTEGER   :: step, nsteps, output_steps, print_steps
  REAL(iwp) :: dt
  ! Plastic iteration
  INTEGER   :: iteration, iteration_limit

  !**************************** CONSTANT PARAMETER ****************************!
  INTEGER :: nst=4    ! Number of Stress/Strain Terms (3: Plane Stress/Strain, 4: Axisymmetric, 6: 3-Dimension)
  INTEGER :: nodof=2  ! Number of Degree of Freedom

  !**************************** PREPARATION PHASE *****************************!
  !--- Initiate Program Settings  

  !--- Load Input Data and Get Input/Output Path
  CALL GET_INPUT_DIRECTORY(input_directory, input_filename)
  CALL LOAD_JSON_DATA(input_directory, input_filename, input_json)
  CALL input_json%GET('post_processing.path', output_directory)
  
  !--- Filenames
  CALL input_json%GET('mesh.mesh', mesh_filename)
  CALL input_json%GET('mesh.entity_sets', entity_filename)
  CALL input_json%GET('particles(1).generator.location', particles_filename)

  !--- Mesh and Material Points Setups
  CALL MPMCORE_LOAD_MESH(input_directory, trim(mesh_filename), mesh)
  CALL MPMCORE_LOAD_PARTICLES(input_directory, trim(particles_filename), mbod) 
  CALL LOAD_JSON_DATA(input_directory, trim(entity_filename), entity_json)

  !--- Analysis and Post-Processing Parameter
  CALL input_json%GET('post_processing.output_steps', output_steps)
  CALL input_json%GET('analysis.nsteps', nsteps)
  CALL input_json%GET('analysis.dt', dt)
  CALL input_json%GET('analysis.max_iter', iteration_limit)

  !***************************** MODEL INITIATION *****************************!

  !--- Setup Global Boundary Conditions
  CALL MPMCORE_FORM_GLOBAL_NF(input_json, entity_json, mesh)

  !--- Insert Material Points
  CALL MPMCORE_INSERT_MPS(mbod, mesh)
  mbod%a_ele = 0 ! element id for each material point
  mbod%c_ele = 0 ! total of MP inside each element
  mbod%k_ele = 0 ! total of MP in the domain
  mbod%d_ele = 0 ! activated element array (1 active/0 deactive)

  
  !************************* VARIABLES INITIALIZATION *************************!
  !--- Activate Initial Compulational Grid
  !--- Get Initial Volume and Mass for Each Material Points
  !--- Calculate Initial Conditions and Properties of Material Points
  !--- Print Initial Conditions Visualization (VTK)

  CALL POINT_VIZ(                         &
    input=0,                              &
    coord=mbod%gm_coord,                  &
    a_ins=mbod%a_ins,                     &
    evpt=mbod%epsinvacum,                 &
    m_stress=mbod%m_stress,               &
    m_stress_inc=mbod%m_stress_change,    &
    acc=mbod%m_acc,                       &
    velocity=mbod%m_velocity,             &
    cohesion=mbod%mpcp,                   &
    devstress=mbod%devstress,             &
    meanstress=mbod%mean_stress,          &
    mpyield=mbod%mpyield,                 &
    directory=output_directory            &
  )
  
  CALL PARAVIEW(                          &
    input=0,                              &
    node_type=4,                          &
    coord=mesh%n_coords,                  &
    num=mesh%el_codes,                    &
    directory=output_directory            &
  )
    
  !****************************** TIME STEPPING *******************************!
  !--- Initiate Time Stepping Variables
  print_steps = 0
  TIME_STEPS : DO step = 1, nsteps
    !**************************** Pre-Calculations ****************************!
    PRINT '(A, I8, A, I8)', "MPM Calculation Step :", step, "/", nsteps
    !--- Body Solution
    !--- Stiffness
    !--- Combined Stiffness
    mbod%gm_coord(1, :) = mbod%gm_coord(1, :) + 0.000001_iwp
    mbod%gm_coord(2, :) = mbod%gm_coord(2, :) + 0.000005_iwp
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
    IF (print_steps == output_steps) THEN
      CALL POINT_VIZ(                         &
        input=step,                           &
        coord=mbod%gm_coord,                  &
        a_ins=mbod%a_ins,                     &
        evpt=mbod%epsinvacum,                 &
        m_stress=mbod%m_stress,               &
        m_stress_inc=mbod%m_stress_change,    &
        acc=mbod%m_acc,                       &
        velocity=mbod%m_velocity,             &
        cohesion=mbod%mpcp,                   &
        devstress=mbod%devstress,             &
        meanstress=mbod%mean_stress,          &
        mpyield=mbod%mpyield,                 &
        directory=output_directory            &
      )
      
      CALL PARAVIEW(                          &
        input=step,                           &
        node_type=4,                          &
        coord=mesh%n_coords,                  &
        num=mesh%el_codes,                    &
        directory=output_directory            &
      )

      print_steps = 0
    END IF
    !--- Reallocate Memory
    !--- Determine and Activate New Element for Next Time Step
  END DO TIME_STEPS
END PROGRAM Tied_Free_Field_MPM