MODULE IO
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp

  CONTAINS

  SUBROUTINE IO_GET_INPUT_DIRECTORY(directory, filename)
    !
    ! This subroutine get the input filename and directory from 
    ! command-line arguments. Directory is placed after `-f` flag and 
    ! the filename is after `-i` flag
    !
    IMPLICIT NONE
    CHARACTER(*), INTENT(INOUT) :: directory, filename
    CHARACTER(16) :: arg
    INTEGER   :: i
    DO i = 1, command_argument_count()
      CALL get_command_argument(i, arg)
      IF (trim(arg) == "-f") THEN
        CALL get_command_argument(i+1, directory)
      END IF
      IF (trim(arg) == "-i") THEN
        CALL get_command_argument(i+1, filename)
      END IF
    END DO
  END SUBROUTINE IO_GET_INPUT_DIRECTORY


  SUBROUTINE IO_LOAD_JSON(directory, filename, json_output)
    USE JSON_MODULE
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: directory, filename
    TYPE(json_file), INTENT(OUT) :: json_output
    CALL json_initialize()
    CALL json_output%load_file(filename=trim(directory)//trim(filename))
  END SUBROUTINE IO_LOAD_JSON


  SUBROUTINE IO_LOAD_MESH(directory, filename, g_coords, g_elements)
    !
    ! Read formatted mesh file and obtain both node coordinates and 
    ! element connectivity matrix
    !
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: directory, filename
    REAL(iwp), ALLOCATABLE, INTENT(OUT) :: g_coords(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: g_elements(:,:)
    CHARACTER(512) :: readline
    INTEGER :: startline, i
    INTEGER :: n_nodes, n_elements, node_types
    ! check starting point
    OPEN(10, FILE=trim(directory)//trim(filename), STATUS='OLD')
    startline = 1
    DO
      READ(10, '(A)') readline
      IF (readline(1:1) /= "!") EXIT
      startline = startline + 1
    END DO

    ! Skipping commented rows
    rewind(10)
    DO i = 1, startline - 1
      READ(10, *)
    END DO

    ! Read 'number of nodes' and 'number of elements' and allocate memory
    READ(10, *) n_nodes, n_elements, node_types
    ALLOCATE(g_coords(2,n_nodes), g_elements(node_types, n_elements))

    ! Read nodes coordinate
    DO i = 1, n_nodes
      READ(10, *) g_coords(:, i)
    END DO

    ! Read elements connectivity number
    DO i = 1, n_elements
      READ(10, *) g_elements(:, i)
    END DO
    CLOSE(10)
  END SUBROUTINE IO_LOAD_MESH


  SUBROUTINE IO_LOAD_PARTICLE(directory, filename, g_mpcoords)
    !
    !
    !
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: directory, filename
    REAL(iwp), ALLOCATABLE, INTENT(OUT) :: g_mpcoords(:,:)
    CHARACTER(512) :: readline
    INTEGER :: startline, i
    INTEGER :: n_particles
    ! check starting point
    OPEN(10, FILE=trim(directory)//trim(filename), STATUS='OLD')
    startline = 1
    DO
      READ(10, '(A)') readline
      IF (readline(1:1) /= "!") EXIT
      startline = startline + 1
    END DO

    ! Skipping commented rows
    rewind(10)
    DO i = 1, startline - 1
      READ(10, *)
    END DO

    ! Read 'number of particles'
    READ(10, *) n_particles
    ALLOCATE(g_mpcoords(2,n_particles))

    ! Read material points coordinate
    DO i = 1, n_particles
      READ(10, *) g_mpcoords(:, i)
    END DO
  END SUBROUTINE IO_LOAD_PARTICLE


  SUBROUTINE IO_POINT_VIZ(input,coord,a_ins,evpt,m_stress,m_stress_inc,acc,     &
    velocity,cohesion,devstress,meanstress,mpyield,directory,argv)
    !
    ! SUBROUTINE used to save visualise outputs to Paraview format
    ! https://dav.lbl.gov/archive/NERSC/Software/express/help6.1/help/reference/dvmac/UCD_Form.htm
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: input
    REAL(iwp), INTENT(IN) :: coord(:,:)
    REAL(iwp), INTENT(IN) :: a_ins(:,:), evpt(:), m_stress(:,:),                &
      m_stress_inc(:,:), acc(:,:), velocity(:,:), cohesion(:), devstress(:),    &
      meanstress(:), mpyield(:)
    CHARACTER(*), OPTIONAL, INTENT(IN):: argv, directory
    
    CHARACTER(128) :: cnumber, s_directory="Results/", s_argv="Particles"
    INTEGER :: i, nmps

    !--- File setup
    if (present(directory)) s_directory = directory
    if (present(argv)) s_argv = argv
    CALL EXECUTE_COMMAND_LINE("mkdir -p "//trim(s_directory))
    WRITE(cnumber, '(i8.8)') input
    OPEN(10,FILE=trim(s_directory)//trim(s_argv)//"_"//trim(cnumber)//'.inp')

    !--- Variables
    nmps = ubound(coord, 2)

    !--- File Comments
    write(10,*) '#'
    write(10,*) '# Simple AVS UCD File'
    write(10,*) '#'

    !--- General Structure
    ! Number of nodes, number of cells, number of data items per node, 
    ! number of data items per cell, model number (always 0)
    write(10,'(2i7.1,a)') nmps, nmps, ' 3 0 0' 

    !--- Node Coordinates
    do i=1,nmps ! coordinates
      write(10,'(i5.1,3g13.5)') i, coord(:,i), zero
    end do
    
    !--- Cell Descriptions
    do i=1,nmps ! points
        write(10,'(2i5.1,a,i5.1)') i, 0, '1 pt ', i
    end do

    !--- Node Based Data Descriptions
    write(10,*) '10 2 1 4 4 2 2 1 1 1 1'
    write(10,*) 'Displacement, meters'
    write(10,*) 'Plastic Strains, meters'
    write(10,*) 'Stress, kPa'
    write(10,*) 'Stress Increase, kPa'
    write(10,*) 'Acceleration, m2/s'
    write(10,*) 'Velocity, m/s'
    write(10,*) 'Cohesion, kPa'
    write(10,*) 'Dev-stress, no-unit'
    write(10,*) 'Mean-stress, no-unit'
    write(10,*) 'Yield-value, no-unit'
    
    !--- Node Based Data
    DO i=1,nmps
      WRITE(10,'(i5.1,19f13.3)') i, a_ins(:,i), evpt(i), m_stress(:,i),         &
        m_stress_inc(:,i), acc(:,i), velocity(:,i), cohesion(i), devstress(i),  &
        meanstress(i), mpyield(i)
    END DO

    CLOSE(10)
  END SUBROUTINE IO_POINT_VIZ


  SUBROUTINE IO_PARAVIEW(input,node_type,coord,num,nf,kdiag,diag,loads,         &
    ddylds,gravlo,d1x1,d2x1,vcm,fcont,normals,f_fint,kv,mv,directory,argv)
    !
    ! Subroutine used to save visualization output of computational mesh
    ! https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: input, node_type
    INTEGER, INTENT(IN)    :: num(:,:)
    REAL(iwp), INTENT(IN)  :: coord(:,:)
    CHARACTER(*), OPTIONAL :: directory, argv
    INTEGER, OPTIONAL, INTENT(IN)   :: nf(:,:), kdiag(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: diag(:), loads(:), ddylds(:), gravlo(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: d1x1(:), d2x1(:), vcm(:), fcont(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: normals(:,:), f_fint(:), kv(:), mv(:)

    CHARACTER(128) :: cnumber, s_directory="Results/", s_argv="Paraview"
    INTEGER :: i, iel, nels, nn

    !--- File setup
    if (present(directory)) s_directory = directory
    if (present(argv)) s_argv = argv
    CALL EXECUTE_COMMAND_LINE("mkdir -p "//trim(s_directory))
    WRITE(cnumber, '(i8.8)') input
    OPEN(10,FILE=trim(s_directory)//trim(s_argv)//"_"//trim(cnumber)//'.vtk')

    !--- Variables
    nn = ubound(coord, 2)
    nels = ubound(num, 2)
    
    !--- File Description Parts (cf. VTK File Format)
    WRITE(10,'(a)')'# vtk DataFile Version 3.0' ! Header
    WRITE(10,'(a)')"vtk output"                 ! Title
    WRITE(10,'(a)')"ASCII"                      ! Data type
    WRITE(10,'(a)')""                           !
    WRITE(10,'(a)')"DATASET UNSTRUCTURED_GRID"  ! Geometry/Topology

    !--- Datasets: Geometry
    ! Node Coordinates
    WRITE(10,'(1A6,1I5,1A7)')"POINTS", nn, "float"
    DO i=1, nn
        WRITE(10,'(3f9.4)') coord(:,i), zero
    END DO
    ! Cell Constructs
    WRITE(10,'(a)')""
    SELECT CASE(node_type)
      CASE(4)
        WRITE(10,'(1A6,2I6)')"CELLS", nels, nels*(1+4)
        DO iel=1, nels
            WRITE(10,'(5I5)') node_type,  &
              num(1,iel)-1, num(4,iel)-1, num(3,iel)-1, num(2,iel)-1
        END DO

        WRITE(10,'(a)')""
        WRITE(10,'(1A10,1I5)')"CELL_TYPES", nels
        DO iel = 1 , nels
            WRITE(10,*)"9"
        END DO

      CASE DEFAULT
        WRITE(*,*)"wrong number of nodes input in paraview"
    END SELECT

    IF (.not. present(nf)) RETURN
    !--- Datasets: Node Values
    WRITE(10,'(a)')""
    WRITE(10,'(1A10,1I5)')"POINT_DATA", nn
    NODAL_MASS: IF (present(diag)) THEN
      WRITE(10,'(a)')"vectors mass float "
      DO i=1, nn
          IF(nf(1,i)==0 .and. nf(2,i)==0) WRITE(10,'(3f15.4)')zero, zero, zero
          IF(nf(1,i)>0  .and. nf(2,i)==0) WRITE(10,'(3f15.4)')diag(nf(1,i)+1), zero, zero
          IF(nf(1,i)==0 .and. nf(2,i)>0 ) WRITE(10,'(3f15.4)')zero, diag(nf(2,i)+1), zero
          IF(nf(1,i)>0  .and. nf(2,i)>0 ) WRITE(10,'(3f15.4)')diag(nf(1,i)+1), diag(nf(2,i)+1), zero
      END DO
    END IF NODAL_MASS

    F_TOTAL: IF (present(loads)) THEN
      WRITE(10,'(a)')"vectors F_Total float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')loads(nf(:,i)+1), zero
      END DO
    END IF F_TOTAL

    F_INTERNAL: IF (present(ddylds)) THEN
      WRITE(10,'(a)')"vectors F_Internal float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')ddylds(nf(:,i)+1), zero
      END DO
    END IF F_INTERNAL

    F_EXTERNAL: IF (present(gravlo)) THEN
      WRITE(10,'(a)')"vectors F_External float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')gravlo(nf(:,i)+1), zero
      END DO
    END IF F_EXTERNAL

    F_KINETIC: IF (present(vcm)) THEN
      WRITE(10,'(a)')"vectors F_Kinetic float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')vcm(nf(:,i)+1), zero
      END DO
    END IF F_KINETIC

    IF (present(f_fint)) THEN
      WRITE(10,'(a)')"vectors ffint float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')f_fint(nf(:,i)+1), zero
      END DO
    END IF

    IF (present(normals)) THEN
      WRITE(10,'(a)')"vectors normals float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')normals(1,i), normals(2,i), zero
      END DO 
    END IF

    IF (present(fcont)) THEN
      WRITE(10,'(a)')"vectors Fcont float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')fcont(nf(:,i)+1), zero
      END DO
    END IF

    NODAL_VELOCITY: IF (present(d1x1)) THEN
      WRITE(10,'(a)')"vectors Velocity float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')d1x1(nf(:,i)+1), zero
      END DO
    END IF NODAL_VELOCITY
    
    NODAL_ACC: if (present(d2x1)) THEN
      WRITE(10,'(a)')"vectors Acceleration float "
      DO i=1, nn
          WRITE(10,'(3f15.6)')d2x1(nf(:,i)+1), zero
      END DO 
    END IF NODAL_ACC

    STIFFNESS_MATRIX: if (present(kv) .and. present(kdiag)) THEN
      WRITE(10,'(a)')"vectors KM_Stiffness float "
      DO i=1, nn
          IF(nf(1,i)==0 .and. nf(2,i)==0) WRITE(10,'(3f15.4)') zero, zero, zero
          IF(nf(1,i)>0  .and. nf(2,i)==0) WRITE(10,'(3f15.4)') kv(kdiag(nf(1,i))), zero, zero
          IF(nf(1,i)==0 .and. nf(2,i)>0 ) WRITE(10,'(3f15.4)') zero, kv(kdiag(nf(2,i))), zero
          IF(nf(1,i)>0  .and. nf(2,i)>0 ) WRITE(10,'(3f15.4)') kv(kdiag(nf(1,i))), kv(kdiag(nf(2,i))), zero
      END DO
    END IF STIFFNESS_MATRIX

    MASS_MATRIX: if (present(mv) .and. present(kdiag)) THEN
      WRITE(10,'(a)')"vectors Mv_Mass float "
      DO i=1, nn
          IF(nf(1,i)==0 .and. nf(2,i)==0) WRITE(10,'(3f15.4)') zero, zero, zero
          IF(nf(1,i)>0  .and. nf(2,i)==0) WRITE(10,'(3f15.4)') mv(kdiag(nf(1,i))), zero, zero
          IF(nf(1,i)==0 .and. nf(2,i)>0 ) WRITE(10,'(3f15.4)') zero, mv(kdiag(nf(2,i))), zero
          IF(nf(1,i)>0  .and. nf(2,i)>0 ) WRITE(10,'(3f15.4)') mv(kdiag(nf(1,i))),mv(kdiag(nf(2,i))), zero
      END DO
    END IF MASS_MATRIX

    ! Close File
    CLOSE(10)
  END SUBROUTINE IO_PARAVIEW

END MODULE IO