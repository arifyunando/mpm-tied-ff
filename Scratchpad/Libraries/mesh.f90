MODULE CLASS_MESH
  USE JSON_MODULE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)

  TYPE::mesh
    INTEGER:: id
    INTEGER:: mesh_ind  ! Indicates mesh type (+:rectilinear, -:isoparametric)
                        ! For rectilinear (0-50: Cartesian; 51-100: Non-Cartesian)
    LOGICAL:: is_initiated=.false.
    
    ! Constants
    INTEGER:: ndim          ! number of problem dimension
    INTEGER:: ndof, nodof   ! number of element degree of freedom, and node dof
    INTEGER:: nn, nels, nod ! number of nodes, elements, and node/element
    INTEGER:: nip           ! number of integration points
    INTEGER:: nomps         ! number of material points per cell initially

    ! Tracking variables
    REAL(iwp),ALLOCATABLE:: g_coord(:,:) ! Node Global Coordinates
    INTEGER,ALLOCATABLE:: g_num(:,:) ! Node index that build up an element
    
    ! Mesh Geometries
    REAL(iwp),ALLOCATABLE:: cellsize(:,:)      ! cellsize of each cell
    INTEGER,ALLOCATABLE:: v_neighbour_id(:,:)  ! list of neighbouring cell indexes
    INTEGER,ALLOCATABLE:: n_neighbour_id(:)    

    ! Computation Indexes
    INTEGER,ALLOCATABLE:: g(:,:)   ! global steering vector of the mesh
    INTEGER,ALLOCATABLE:: nf(:,:)  ! Track degree of freedoms
    
    ! Sparse Matrix variables
    INTEGER,ALLOCATABLE:: kdiag(:) ! Indexes in the diagonal of skyline matrix
    
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

    PROCEDURE :: LOAD => p_LOAD_MESH
    PROCEDURE :: GENERATE => p_GENERATE_MESH
    ! PROCEDURE :: FORM_GLOBAL_NF => CLASS_MESH_FORM_GLOBAL_NF
  END TYPE

  CONTAINS

!****************************** PUBLIC FUNCTIONS ******************************!
!                                                                              !
  
  SUBROUTINE p_LOAD_MESH(this, directory, input_json, nodof, nst)
    !
    ! Constructor of The MPM_GRID Class Object
    !
    USE FUNCTIONS
    USE JSON_MODULE
    USE IO
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN)       :: directory
    TYPE(json_file), INTENT(INOUT) :: input_json
    CLASS(mesh), INTENT(INOUT) :: this
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
    CALL IO_LOAD_MESH(trim(directory),trim(filename),this%g_coord,this%g_num)
    
    ! Initiate mesh variables
    CALL m_INITIATE_MESH(this,this%g_coord,this%g_num,nst,nodof)

    ! Apply Global Boundary Condition
    ! CALL CLASS_MESH_FORM_GLOBAL_NF(this, directory, input_json)
  END SUBROUTINE p_LOAD_MESH


  SUBROUTINE p_GENERATE_MESH(this,nx,ny,w1,h1,node,  &
    offsetx,offsety,nst,ndim,nodof)
    USE GEOMETRY_GENERATOR
    USE FUNCTIONS
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nx,ny,node
    INTEGER,INTENT(IN)::nst,ndim,nodof
    REAL(iwp),INTENT(IN)::w1,h1
    CLASS(mesh),INTENT(INOUT)::this
    INTEGER,OPTIONAL,INTENT(IN)::offsetx,offsety
    INTEGER::offx=0,offy=0
    ! Generate Rectangular FE Mesh
    IF(PRESENT(offsetx)) offx = offsetx
    IF(PRESENT(offsety)) offx = offsety
    CALL RECTANGULAR_2D(this%g_coord,this%g_num,this%nn,this%nels,&
                        nx,ny,w1,h1,node,offx,offy)
    this%mesh_ind=4
    CALL m_INITIATE_MESH(this,this%g_coord,this%g_num,nst,nodof)
  END SUBROUTINE p_GENERATE_MESH


  SUBROUTINE p_APPLY_DISPLACEMENT_BOUNDARY()
    IMPLICIT NONE

  END SUBROUTINE p_APPLY_DISPLACEMENT_BOUNDARY

!                                                                              !
!****************************** PRIVATE FUNCTIONS *****************************!
!                                                                              !
  
  SUBROUTINE m_INITIATE_MESH(this,g_coord,g_num,nst,nodof)
    IMPLICIT NONE
    CLASS(mesh),INTENT(INOUT)::this
    REAL(iwp),ALLOCATABLE,INTENT(IN)::g_coord(:,:)
    INTEGER,ALLOCATABLE,INTENT(IN)::g_num(:,:)
    INTEGER,INTENT(IN)::nodof,nst
    ! Local Variable
    INTEGER::ndim,nn,nels,nod,ndof
    !--- Determine and Assign Mesh Constants
    ndim = ubound(g_coord,1)
    nn = ubound(g_coord,2)
    nels = ubound(g_num,2)
    nod = ubound(g_num,1)
    ndof = nod*nodof
    this%ndim=ndim; this%nn=nn; this%nels=nels
    this%nod=nod; this%nodof=nodof; this%ndof=ndof
    this%g_coord=g_coord; this%g_num=g_num
    !--- Allocate Memory for variables with known shape
    ALLOCATE(                    &
      this%nf(nodof,nn),         &
      this%g(ndof,nels)          &
    )
    !--- Locate Cell Neighbours
    CALL m_LOCATE_NEIGHBOUR(this)
    !--- Determmine Cellsize
    CALL m_CALCULATE_CELLSIZE(this)

    ! Mark initiation
    this%is_initiated=.true.
  END SUBROUTINE m_INITIATE_MESH

  
  SUBROUTINE m_FORM_STEERING_VECTOR(this)
    USE FUNCTIONS
    IMPLICIT NONE
    CLASS(mesh), INTENT(INOUT) :: this
    INTEGER,ALLOCATABLE::num(:),g(:)
    INTEGER::iel
    ALLOCATE(num(this%nod), g(this%nodof))
    DO iel=1, this%nels
      num = this%g_num(:,iel)
      CALL NUM_TO_G(num, this%nf, g)
      this%g(:,iel) = g
    END DO
  END SUBROUTINE m_FORM_STEERING_VECTOR


  SUBROUTINE m_CALCULATE_CELLSIZE(this)
    !
    ! Returns the cell size of each cell that build up the mesh 
    ! Note:
    !   Currently only assume constant cellsize
    !
    IMPLICIT NONE
    CLASS(mesh), INTENT(INOUT) :: this
    SELECT CASE(this%mesh_ind)
    CASE (1:50) ! Cartesian Grid / Constant Cellsize
      ALLOCATE(this%cellsize(this%ndim,this%nels))
      this%cellsize = ABS(this%g_coord(2, 1) - this%g_coord(2, 2))
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
  END SUBROUTINE m_CALCULATE_CELLSIZE


  SUBROUTINE m_LOCATE_NEIGHBOUR(this)
    IMPLICIT NONE
    CLASS(mesh),INTENT(INOUT)::this
    ! Local Variable
    INTEGER,ALLOCATABLE::neighbour_id(:,:),neighbour_count(:)
    INTEGER,ALLOCATABLE::this_element_nodes(:),that_element_nodes(:)
    INTEGER::this_node,that_node
    INTEGER::this_ele,that_ele,p,q,i
    
    ALLOCATE(                     &
      neighbour_id(8,this%nels),  & ! assume max neighbour is 8 cells
      neighbour_count(this%nels)  &
    )

    !--- Determine neighbour_id and neighbour_count
    neighbour_count=1
    DO this_ele=1,this%nels
      this_element_nodes=this%g_num(:,this_ele)
      DO that_ele=this_ele+1,this%nels
        that_element_nodes=this%g_num(:,that_ele)
        ! Check whether that_ele is already a neighbour
        DO i=1,neighbour_count(this_ele)
          IF (neighbour_id(i,this_ele) == that_ele) GO TO 100
        END DO
        ! Check whether that_ele is a neighbour
        DO p=1,this%nod
          this_node=this_element_nodes(p)
          DO q=1,this%nod
            that_node=that_element_nodes(q)
            IF (this_node == that_node) THEN
              neighbour_id(neighbour_count(this_ele),this_ele)=that_ele
              neighbour_id(neighbour_count(that_ele),that_ele)=this_ele
              neighbour_count(this_ele) = neighbour_count(this_ele) + 1
              neighbour_count(that_ele) = neighbour_count(that_ele) + 1
              GO TO 100 ! Check next element
            END IF
          END DO
        END DO
        100 CONTINUE
      END DO
    END DO

    !--- Assign neighbours to collections
    this%n_neighbour_id = neighbour_count - 1
    this%v_neighbour_id = neighbour_id
  END SUBROUTINE m_LOCATE_NEIGHBOUR


  SUBROUTINE m_GET_DISPLACEMENT_BOUNDARY_NODES()
  END SUBROUTINE m_GET_DISPLACEMENT_BOUNDARY_NODES

  ! SUBROUTINE m_FORM_GLOBAL_NF(this,nodes_per_direction,node_sets)
  !   USE FUNCTIONS
  !   IMPLICIT NONE
  !   CLASS(mesh), INTENT(INOUT) ::this
  !   INTEGER,ALLOCATABLE,INTENT(IN)::nodes_per_direction(:),node_sets(:,:)
  !   INTEGER :: node,dir
  !   this%nf=1
  !   DO dir=1,size(nodes_per_direction)
  !     DO node=1,nodes_per_direction(dir)
  !       this%nf(dir,node)=0
  !     END DO
  !   END DO
  !   CALL FORMNF(this%nf)
  ! END SUBROUTINE m_FORM_GLOBAL_NF

  ! SUBROUTINE CLASS_MESH_CONSTRUCT_STENCILS(this)
  !   IMPLICIT NONE
  !   CLASS(mesh), INTENT(INOUT) :: this
  ! END SUBROUTINE CLASS_MESH_CONSTRUCT_STENCILS


  ! SUBROUTINE CLASS_MESH_CALCULATE_PARTICLE_SUPPORT_LENGTH(this)
  !   !
  !   ! Calculate lp based on the cellsize and the number of particles
  !   ! in direction
  !   !
  !   IMPLICIT NONE
  !   CLASS(mesh), INTENT(INOUT) :: this
  
  ! END SUBROUTINE CLASS_MESH_CALCULATE_PARTICLE_SUPPORT_LENGTH



END MODULE CLASS_MESH