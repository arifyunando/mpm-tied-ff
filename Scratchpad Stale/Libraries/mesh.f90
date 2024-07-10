MODULE CLASS_MESH
  USE JSON_MODULE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)

  TYPE::mpm_grid
    INTEGER:: id
    INTEGER:: mesh_ind  ! Indicates mesh type (+:rectilinear, -:isoparametric)
                        ! For rectilinear (0-50: Cartesian; 51-100: Non-Cartesian)
    ! Constants
    INTEGER:: ndim          ! number of problem dimension
    INTEGER:: ndof, nodof   ! number of element degree of freedom, and node dof
    INTEGER:: nn, nels, nod ! number of nodes, elements, and node/element
    INTEGER:: nip           ! number of integration points
    INTEGER:: nomps         ! number of material points per cell initially

    ! Tracking variables
    REAL(iwp),ALLOCATABLE:: g_coord(:,:) ! Node Global Coordinates
    INTEGER,ALLOCATABLE:: num(:,:) ! Node index that build up an element
    INTEGER,ALLOCATABLE:: g(:,:)   ! global steering vector of the mesh
    INTEGER,ALLOCATABLE:: nf(:,:)  ! Track degree of freedoms
    
    ! Mesh Geometries
    REAL(iwp),ALLOCATABLE:: cellsize(:,:) ! cellsize of each cell block
    INTEGER,ALLOCATABLE:: neighb(:,:)     ! list of neighbouring cell indexes

    ! Indexes
    INTEGER,ALLOCATABLE:: MPPE(:)  ! Initial Material Point per Element

    ! Sparse Matrix variables
    INTEGER,ALLOCATABLE:: kdiag(:) ! Indexes in the diagonal of skyline matrix

    ! Material Point Tracers
    INTEGER,ALLOCATABLE:: c_ele(:) ! number of MPs inside an element
    INTEGER,ALLOCATABLE:: d_ele(:) ! activated element array (1:active; 0:deactive)
    INTEGER,ALLOCATABLE:: k_ele(:) ! total of MP in the domain / Accum. of c_ele
    
    ! Nodal Properties (printable to VTK)
    REAL(iwp),ALLOCATABLE:: eld(:) ! element nodal displacement 
    REAL(iwp),ALLOCATABLE:: vel(:) ! element nodal velocity 
    REAL(iwp),ALLOCATABLE:: acc(:) ! element nodal acceleration
    REAL(iwp),ALLOCATABLE:: f_int(:) ! nodal internal force (ddylds)
    REAL(iwp),ALLOCATABLE:: f_ext(:) ! nodal external force (gravlo)
    REAL(iwp),ALLOCATABLE:: f_kin(:) ! nodal kinetic force (vcm)
    REAL(iwp),ALLOCATABLE:: kv(:) ! Global 'effective' stiffness matrix
    REAL(iwp),ALLOCATABLE:: mv(:) ! Global 'effective' mass matrix
    CONTAINS

    PROCEDURE :: LOAD_MESH => CLASS_MESH_LOAD_MESH
    PROCEDURE :: FORM_GLOBAL_NF => CLASS_MESH_FORM_GLOBAL_NF
  END TYPE

  CONTAINS

  SUBROUTINE CLASS_MESH_LOAD_MESH(this, directory, input_json, nodof, nst)
    !
    ! Constructor of The MPM_GRID Class Object
    !
    USE FUNCTIONS
    USE JSON_MODULE
    USE IO
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN)       :: directory
    TYPE(json_file), INTENT(INOUT) :: input_json
    CLASS(mpm_grid), INTENT(INOUT) :: this
    INTEGER, OPTIONAL, INTENT(INOUT) :: nodof, nst
    CHARACTER(:), ALLOCATABLE :: filename, cell_type
    INTEGER :: def_nodof=2, def_nst=4, ind

    ! Default Arguments
    if (present(nodof)) def_nodof = nodof
    if (present(nst)) def_nst = nst
    
    ! Get Mesh Filename
    CALL input_json%GET('mesh.mesh', filename)

    ! Get Cell Type
    CALL input_json%GET('mesh.cell_type', cell_type)
    if (cell_type == "ED2Q4G") this%mesh_ind = 4

    ! Get mesh nodal locations
    CALL IO_LOAD_MESH(trim(directory), trim(filename), this%g_coord, this%num)
    
    ! Determine Mesh Variable
    this%ndim = ubound(this%g_coord, 1)
    this%nn = ubound(this%g_coord, 2)
    this%nels = ubound(this%num, 2)
    this%nod = ubound(this%num, 1)
    this%nodof = def_nodof
    this%ndof = this%nod * this%nodof

    ! Allocate Memory for variables with known shape
    ALLOCATE(                       &
      this%nf(def_nodof, this%nn),  &
      this%g(this%ndof,this%nels),  &
      this%c_ele(this%nels),        &
      this%d_ele(this%nels),        &
      this%k_ele(0:this%nels),      &
      this%neighb(this%nels,8)      &
    )

    ! Apply Global Boundary Condition
    CALL CLASS_MESH_FORM_GLOBAL_NF(this, directory, input_json)

    ! Create Global Steering Vector
    CALL CLASS_MESH_FORM_STEERING_VECTOR(this)

    ! Determmine Cellsize
    CALL CLASS_MESH_CALCULATE_CELLSIZE(this)
  END SUBROUTINE CLASS_MESH_LOAD_MESH


  SUBROUTINE CLASS_MESH_FORM_GLOBAL_NF(this, directory, input_json)
    USE JSON_MODULE
    USE FUNCTIONS
    USE IO
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN)       :: directory
    TYPE(json_file), INTENT(INOUT) :: input_json
    CLASS(mpm_grid), INTENT(INOUT) :: this
    TYPE(json_file) :: entity_json
    LOGICAL :: found=.true.
    INTEGER :: i, nset_id, dir
    INTEGER, ALLOCATABLE :: nodes(:)
    CHARACTER(:), ALLOCATABLE :: entity_filename
    CHARACTER(1) :: ic

    CALL input_json%GET('mesh.entity_sets', entity_filename)
    CALL IO_LOAD_JSON(trim(directory), trim(entity_filename), entity_json)

    this%nf = 1
    i = 1
    DO WHILE(found)
      write(ic, '(I1)') i
      CALL input_json%GET(                                                      &
      'mesh.boundary_conditions.displacement_constraints('//ic//').nset_id',    &
      nset_id, found                                                            &
      )
      CALL input_json%GET(                                                      &
      'mesh.boundary_conditions.displacement_constraints('//ic//').dir',        &
      dir, found                                                                &
      )
      IF (found) THEN
        write(ic, '(I1)') nset_id
        CALL entity_json%GET('node_sets('//ic//').set', nodes)
        this%nf(dir,nodes) = 0
      END IF
      i = i+1
    END DO
    CALL FORMNF(this%nf)
  END SUBROUTINE CLASS_MESH_FORM_GLOBAL_NF

  
  SUBROUTINE CLASS_MESH_FORM_STEERING_VECTOR(this)
    USE FUNCTIONS
    IMPLICIT NONE
    CLASS(mpm_grid), INTENT(INOUT) :: this
    INTEGER,ALLOCATABLE::num(:),g(:)
    INTEGER::iel
    ALLOCATE(num(this%nod), g(this%nodof))
    DO iel=1, this%nels
      num = this%num(:,iel)
      CALL NUM_TO_G(num, this%nf, g)
      this%g(:,iel) = g
    END DO
  END SUBROUTINE CLASS_MESH_FORM_STEERING_VECTOR


  SUBROUTINE CLASS_MESH_CALCULATE_CELLSIZE(this)
    !
    ! Returns the cell size of each cell that build up the mesh 
    ! Note:
    !   Currently only assume constant cellsize
    !
    IMPLICIT NONE
    CLASS(mpm_grid), INTENT(INOUT) :: this
    SELECT CASE(this%mesh_ind)
    CASE (1:50) ! Cartesian Grid / Constant Cellsize
      ALLOCATE(this%cellsize(this%ndim, this%nels))
      this%cellsize = ABS(this%g_coord(1, 1) - this%g_coord(1, 2))
    CASE (51:) ! Rectilinear Grid / Constant Cellsize in 1-dir
      ! get all the node coordinates of the stencils
      ! nstencil = UBOUND(stencils, 1)
      ! ALLOCATE(coords(ndim,nstencil))
      ! coords = g_coord(:,stencils(:,s))
      ! ! Accumulate cellsize and average over all stencils
      ! count=0
      ! cellsize = 0.0_iwp
      ! DO i=2,nstencil
      !   IF (stencils(1,s) == stencils(i,s)) CYCLE
      !   cellsize = cellsize + SQRT(SUM((coords(:,1) - coords(:,i))**2))
      !   count = count + 1
      ! END DO
      ! cellsize = cellsize/count
    CASE (:-1) ! Isoparametric
    CASE DEFAULT
      WRITE(*, *) "Mesh shape index is not valid"
    END SELECT
  END SUBROUTINE CLASS_MESH_CALCULATE_CELLSIZE


  ! TODOs
  ! SUBROUTINE CLASS_MESH_LOOK_NEIGHBOUR(neighb)
  !   IMPLICIT NONE
  !   INTEGER,INTENT(OUT)::neighb(:,:) ! shape(max_neighbour, n_elements)
  ! END SUBROUTINE CLASS_MESH_LOOK_NEIGHBOUR


  ! SUBROUTINE CLASS_MESH_CONSTRUCT_STENCILS(this)
  !   IMPLICIT NONE
  !   CLASS(mpm_grid), INTENT(INOUT) :: this
  ! END SUBROUTINE CLASS_MESH_CONSTRUCT_STENCILS


  ! SUBROUTINE CLASS_MESH_CALCULATE_PARTICLE_SUPPORT_LENGTH(this)
  !   !
  !   ! Calculate lp based on the cellsize and the number of particles
  !   ! in direction
  !   !
  !   IMPLICIT NONE
  !   CLASS(mpm_grid), INTENT(INOUT) :: this
  
  ! END SUBROUTINE CLASS_MESH_CALCULATE_PARTICLE_SUPPORT_LENGTH



END MODULE CLASS_MESH