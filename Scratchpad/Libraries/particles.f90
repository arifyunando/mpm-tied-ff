MODULE CLASS_PARTICLE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp

  TYPE::particles
    INTEGER::id
    CHARACTER(16)::name
    LOGICAL::is_initiated=.false.
    
    ! Constants
    INTEGER::nn    ! Number of nodes
    INTEGER::nmps  ! Number of material points
    INTEGER::ndim  ! Dimensions
    INTEGER::nst   ! Number of stress/strain terms
    INTEGER::nels  ! Number of FE Elements
    INTEGER::nodof ! Number of DOF per nodes
    INTEGER::ndof  ! Number of DOF per element

    ! Commons Variable
    ! INTEGER,ALLOCATABLE::nf(:,:)        ! Node-Equation indexing
    ! INTEGER,ALLOCATABLE::neighbour(:,:) ! list of neighbouring cell_id
    ! INTEGER,ALLOCATABLE::tied_nn(:,:)   ! list of tied nodes
    ! INTEGER,ALLOCATABLE::base_nn(:)     ! list of base nodes
    
    ! Variables to track material points
    REAL(iwp),ALLOCATABLE::gm_coord(:,:) ! MP global coordinates

    ! GIMP support tracker
    INTEGER,ALLOCATABLE::n_support_nodes(:)   ! List of Node ids inside mp support domain
    INTEGER,ALLOCATABLE::v_support_nodes(:,:) 
    INTEGER,ALLOCATABLE::member_elements(:,:) ! list of Elements ids inside mp support domain
    
    ! Varibales to track Finite Element Grid
    REAL(iwp),ALLOCATABLE::g_coord(:,:) ! Mesh node global coordinates
    INTEGER,ALLOCATABLE::g_num(:,:)     ! global cell-nodes ids 
    INTEGER,ALLOCATABLE::g_g(:,:)       ! global equations number vector
    INTEGER,ALLOCATABLE::num(:)         ! Cell-node indexes
    INTEGER,ALLOCATABLE::g(:)           ! Cell equations number vector

    ! Numerical Integration
    REAL(iwp),ALLOCATABLE::mweights(:)   ! Initial Particle Weights (Initial Particle Volume)
    REAL(iwp),ALLOCATABLE::lm_coord(:,:) ! Local Particle Coordinates
    REAL(iwp),ALLOCATABLE::lp(:,:)       ! Local Particle Size (GIMP)

    ! External/Boundary Conditions
    ! REAL(iwp),ALLOCATABLE::ground_acc(:)  ! Ground Acceleration
    ! REAL(iwp),ALLOCATABLE::g_matrix(:) ! Gravitational loading
    
    ! Particle Material Properties
    ! REAL(iwp),ALLOCATABLE::mp_dens(:) ! Particle densitiy
    REAL(iwp),ALLOCATABLE::mpcp(:)    ! Particle cohesion
    
    ! ! Kinetics and Kinematics
    ! REAL(iwp),ALLOCATABLE::eps_acum(:,:)        ! Accumulative particle strain
    REAL(iwp),ALLOCATABLE::m_stress(:,:)        ! Current Particle Stress
    ! REAL(iwp),ALLOCATABLE::m_stress_ini(:,:)    ! Initial Particle Stress
    REAL(iwp),ALLOCATABLE::m_stress_change(:,:) ! Stress Increment
    
    ! ! Constitutive Model
    REAL(iwp),ALLOCATABLE::dee(:,:)       ! Elastic Stiffness Matrix
    REAL(iwp),ALLOCATABLE::mean_stress(:) ! Mean Stress (p-q diagram)
    REAL(iwp),ALLOCATABLE::Devstress(:)   ! Deviatoric Stress (p-q diagram)
    REAL(iwp),ALLOCATABLE::epsinvacum(:)  ! Particle Plastic Strain (?)
    REAL(iwp),ALLOCATABLE::mpyield(:)     ! Yield Function Value (F=0:Plastic, F<0:Elastic)
    
    ! Reverse Mapping Particles
    ! REAL(iwp),ALLOCATABLE::accp(:,:) ! MP present step acceleration
    ! REAL(iwp),ALLOCATABLE::accb(:,:) ! MP future step acceleration
    ! REAL(iwp),ALLOCATABLE::vccp(:,:) ! MP present step velocity
    ! REAL(iwp),ALLOCATABLE::ins(:,:)  ! MP displacement vector
    REAL(iwp),ALLOCATABLE::a_ins(:,:)      ! MP accumulated displacement
    REAL(iwp),ALLOCATABLE::m_acc(:,:)      ! MP acceleration
    REAL(iwp),ALLOCATABLE::m_velocity(:,:) ! MP velocity

    ! Particle Mass and Volume. Constant throughout time stepping
    REAL(iwp),ALLOCATABLE::m_mass(:)       ! MP mass
    REAL(iwp),ALLOCATABLE::m_volume(:)     ! Current particle volume
    
    ! ! Global Body Field
    REAL(iwp),ALLOCATABLE::gravlo(:)  ! Nodal Gravity loading vector
    REAL(iwp),ALLOCATABLE::loads(:)   ! Nodal total forces (before solve), displacement (after solve)
    REAL(iwp),ALLOCATABLE::kv(:)      ! Constructed stiffness matrix
    REAL(iwp),ALLOCATABLE::mv(:)      ! Constructed mass matrix
    REAL(iwp),ALLOCATABLE::mass(:)    ! Diagonal vector of lumped mass matrix (old: diag)
    REAL(iwp),ALLOCATABLE::ddylds(:)  ! Nodal internal Forces (integral of stresses)
    REAL(iwp),ALLOCATABLE::d1x1(:)    ! Nodal velocity vector
    REAL(iwp),ALLOCATABLE::d2x1(:)    ! Nodal acceleration vector
    REAL(iwp),ALLOCATABLE::kinup_d2x1(:) ! Nodal velocity update vector
    REAL(iwp),ALLOCATABLE::kinup_d1x1(:) ! Nodal acceleration update vector
    REAL(iwp),ALLOCATABLE::kinup_Ground_d2x1(:) ! Nodal acceleration update vector due to ground movement
    REAL(iwp),ALLOCATABLE::x1(:)      ! Nodal Displacement Vector
    REAL(iwp),ALLOCATABLE::x1_ini(:)  ! Nodal initical displacement vector
    REAL(iwp),ALLOCATABLE::x1_change(:) ! Difference between x1 and x1_ini

    ! ! Data structure - NEED TO BE REFACTORED TO ALLOW BETTER 3-ARRAY STRUCT
    ! INTEGER::n_nonzero ! number of non-zero elements in the sparse matrix
    INTEGER,ALLOCATABLE::kdiag(:)    ! skyline matrix diagonal index
    ! INTEGER,ALLOCATABLE::ia(:)       ! CSR Form: Row Indexes
    ! INTEGER,ALLOCATABLE::ja(:)       ! CSR Form: Column Indexes
    ! REAL(iwp),ALLOCATABLE::kv_csr(:) ! CSR Form: Constructed Stiffness
    
    ! ! Contact-specific variables
    ! LOGICAL  ::activate=.false.
    ! REAL(iwp)::phi, gimptol

    CONTAINS

    ! Subroutines and Functions
    ! PROCEDURE::LOAD => p_LOAD_PARTICLES
    PROCEDURE::GENERATE => p_GENERATE_PARTICLES
    PROCEDURE::RESET_MP => p_RESET_MP
  END TYPE

  CONTAINS

!****************************** PUBLIC FUNCTIONS ******************************!
!                                                                              !
  
  SUBROUTINE p_GENERATE_PARTICLES(this,nx,ny,w1,h1,node,  &
    offsetx,offsety,nst,ndim,nodof)
    !
    ! Generate Quadrilateral FEM Body with the number of integration points 
    ! similar to its nodes per element number. 
    !
    USE GEOMETRY_GENERATOR
    USE FUNCTIONS
    IMPLICIT NONE
    ! Dummy Variables
    INTEGER,INTENT(IN)::nx,ny,node
    INTEGER,INTENT(IN)::nst,ndim,nodof
    REAL(iwp),INTENT(IN)::w1,h1
    CLASS(particles),INTENT(INOUT)::this
    INTEGER,OPTIONAL,INTENT(IN)::offsetx,offsety
    ! Local Variables
    REAL(iwp),ALLOCATABLE::g_coord(:,:),gm_coord(:,:),m_coord(:,:),coord(:,:)
    REAL(iwp),ALLOCATABLE::points(:,:),weights(:),fun(:)
    INTEGER,ALLOCATABLE::MPPE(:),num(:),g_num(:,:)
    INTEGER::nmps,i,iel,nip,nn,nels
    INTEGER::offx=0,offy=0
    ! Generate Rectangular FE Mesh
    IF(PRESENT(offsetx)) offx = offsetx
    IF(PRESENT(offsety)) offy = offsety
    CALL RECTANGULAR_2D(g_coord,g_num,nn,nels,nx,ny,w1,h1,node,offx,offy)
    ! Determine number of generated particles
    nmps=nels*node
    nip=node
    ALLOCATE(gm_coord(ndim,nmps))
    ! Allocate local variables
    ALLOCATE(MPPE(nip),num(nip),coord(node,ndim),fun(node),&
             points(node,ndim),weights(nip),m_coord(nip,ndim))
    ! Generate Particles from Gaussian Integration Points
    DO i=1,node
      MPPE(i)=i
    END DO
    CALL SAMPLE_GAUSS_IP('quadrilateral',points,weights)
    DO iel=1,nels
      num=g_num(:,iel)
      coord=TRANSPOSE(g_coord(:,num))
      DO i=1,nip
        CALL SHAPE_FUN(fun,points,i)
        m_coord(i,:)=MATMUL(fun,coord)
      END DO
      gm_coord(:,MPPE)=TRANSPOSE(m_coord)
      MPPE(:)=MPPE(:)+nip
    END DO
    ! Initiate particle class
    CALL m_INITIATE_PARTICLES(this,gm_coord,g_coord,nst,nodof)
  END SUBROUTINE p_GENERATE_PARTICLES


  ! SUBROUTINE p_LOAD_PARTICLES(this,index,directory,input_json,    &
  !   ndim,nodof,nst)
  !   USE JSON_MODULE
  !   USE IO
  !   IMPLICIT NONE
  !   CLASS(particles),INTENT(INOUT)::this
  !   INTEGER,INTENT(IN)::index ! particles array index in the input json file
  !   CHARACTER(*),INTENT(IN)::directory
  !   TYPE(json_file),INTENT(INOUT)::input_json
  !   INTEGER,OPTIONAL,INTENT(INOUT)::nodof,nst,ndim
  !   INTEGER::nmps,def_nodof=2,def_nst=4,def_ndim=2
  !   CHARACTER(1)::i_char
  !   CHARACTER(:),ALLOCATABLE::filename,name
  !   ! get particle generator data
  !   WRITE(i_char, '(I1)') index
  !   CALL input_json%GET('particles('//i_char//').generator.name', name)
  !   CALL input_json%GET('particles('//i_char//').generator.location', filename)
  !   this%name=name
  !   ! get particle locations
  !   CALL IO_LOAD_PARTICLE(trim(directory), trim(filename), this%gm_coord)
  !   ! Determine Particles Variable
  !   if (present(ndim)) def_ndim = ndim
  !   if (present(nodof)) def_nodof = nodof
  !   if (present(nst)) def_nst = nst
  !   nmps=ubound(this%gm_coord,2)
  !   CALL m_INITIATE_PARTICLES(this,def_ndim,nmps,def_nst,def_nodof)
  ! END SUBROUTINE p_LOAD_PARTICLES


  SUBROUTINE p_RESET_MP(this)
    !
    ! Reset time-stepping variables
    !
    IMPLICIT NONE
    CLASS(particles), INTENT(INOUT) :: this
    ! Displacement
    this%x1                = zero
    this%x1_ini            = zero
    this%x1_change         = zero
    ! Loads
    this%gravlo            = zero
    this%loads             = zero
    ! Kinematics
    this%d1x1              = zero
    this%d2x1              = zero
    this%kinup_Ground_d2x1 = zero
  END SUBROUTINE p_RESET_MP

!                                                                              !
!****************************** PRIVATE FUNCTIONS *****************************!
!                                                                              !
  
  SUBROUTINE m_INITIATE_PARTICLES(this,gm_coord,g_coord,nst,nodof)
    !
    ! Initiate Particle and Allocate Variables
    !
    CLASS(particles),INTENT(INOUT)::this
    REAL(iwp),ALLOCATABLE,INTENT(IN)::gm_coord(:,:),g_coord(:,:)
    INTEGER,INTENT(IN)::nst,nodof
    ! Local Variable
    INTEGER::nmps,ndim
    !--- Determine and assign particle contants
    nmps=ubound(gm_coord,2)
    ndim=ubound(g_coord,1)
    this%nmps=nmps; this%ndim=ndim; this%nst=nst; this%nodof=nodof
    this%gm_coord=gm_coord; this%g_coord=g_coord
    !--- Allocate Memory for variables that depends on nn, nels, nmps, and nst
    ALLOCATE(this%lm_coord(ndim,nmps))      
    ALLOCATE(this%mweights(nmps))           
    ! ALLOCATE(this%m_mass(nmps))             
    ! ALLOCATE(this%m_volume(nmps))           
    ALLOCATE(this%m_stress(nst,nmps))       
    ALLOCATE(this%m_velocity(nodof,nmps))   
    ALLOCATE(this%a_ins(ndim,nmps))         
    ! ALLOCATE(this%vccp(ndim,nmps))          
    ! ALLOCATE(this%accp(ndim,nmps))          
    ALLOCATE(this%epsinvacum(nmps))         
    ALLOCATE(this%Devstress(nmps))          
    ! ALLOCATE(this%ins(ndim,nmps))           
    ! ALLOCATE(this%eps_acum(nst,nmps))       
    ! ALLOCATE(this%support_nodes(9,nmps))    
    ! ALLOCATE(this%m_stress_ini(nst,nmps))   
    ALLOCATE(this%m_stress_change(nst,nmps))
    ALLOCATE(this%m_acc(nodof,nmps))        
    ALLOCATE(this%mpyield(nmps))            
    ALLOCATE(this%mean_stress(nmps))        
    ALLOCATE(this%lp(ndim,nmps))            
    ALLOCATE(this%member_elements(nmps,4))  
    ! ALLOCATE(this%accb(ndim,nmps))          
    ! ALLOCATE(this%mp_dens(nmps))            
    ALLOCATE(this%mpcp(nmps))
    
    ! Mark initiation
    this%is_initiated=.true.
  END SUBROUTINE m_INITIATE_PARTICLES


END MODULE CLASS_PARTICLE