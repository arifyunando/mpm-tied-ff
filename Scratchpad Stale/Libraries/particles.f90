MODULE CLASS_PARTICLE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp

  TYPE::particles
    INTEGER::id
    CHARACTER(16)::name
    
    ! Constants
    INTEGER::nn    ! Number of nodes
    INTEGER::nmps  ! Number of material points
    INTEGER::ndim  ! Dimensions
    INTEGER::nst   ! Number of stress/strain terms
    INTEGER::nels  ! Number of FE Elements
    INTEGER::nodof ! Number of DOF per nodes
    INTEGER::ndof  ! Number of DOF per element

    ! Commons Variable
    REAL(iwp),ALLOCATABLE::gm_coord(:,:) ! MP global coordinates
    INTEGER,ALLOCATABLE::nf(:,:)        ! Node-Equation indexing
    INTEGER,ALLOCATABLE::neighbour(:,:) ! list of neighbouring cell_id
    INTEGER,ALLOCATABLE::tied_nn(:,:)   ! list of tied nodes
    INTEGER,ALLOCATABLE::base_nn(:)     ! list of base nodes
    
    ! Variables to track material points
    INTEGER::n_active_ele
    INTEGER,ALLOCATABLE::a_ele(:) ! Element index each MP resides 
    INTEGER,ALLOCATABLE::c_ele(:) ! Number of MPs inside a given element
    INTEGER,ALLOCATABLE::d_ele(:) ! Active element (1:active; 0:deactive)
    INTEGER,ALLOCATABLE::k_ele(:) ! Accumulation of c_ele
    INTEGER,ALLOCATABLE::support_nodes(:,:)   ! List of Node ids inside mp support domain
    INTEGER,ALLOCATABLE::member_elements(:,:) ! list of Elements ids inside mp support domain
    
    ! Varibales to track Finite Element Grid
    REAL(iwp),ALLOCATABLE::g_coord(:,:) ! Mesh node global coordinates
    INTEGER,ALLOCATABLE::g_num(:,:)     ! global element node ids 
    INTEGER,ALLOCATABLE::num(:)         ! Mesh node indexes
    INTEGER,ALLOCATABLE::g_g(:,:)       ! global steering factor
    INTEGER,ALLOCATABLE::g(:)           ! local steering factor

    ! Numerical Integration
    REAL(iwp),ALLOCATABLE::mweights(:)   ! Initial Particle Weights (Initial Particle Volume)
    REAL(iwp),ALLOCATABLE::lm_coord(:,:) ! Local Particle Coordinates
    REAL(iwp),ALLOCATABLE::lp(:,:)       ! Local Particle Size (GIMP)

    ! External/Boundary Conditions
    REAL(iwp),ALLOCATABLE::ground_acc(:)  ! Ground Acceleration
    REAL(iwp),ALLOCATABLE::g_matrix(:) ! Gravitational loading
    
    ! Particle Material Properties
    REAL(iwp),ALLOCATABLE::mp_dens(:) ! Particle densitiy
    REAL(iwp),ALLOCATABLE::mpcp(:)    ! Particle cohesion
    
    ! Kinetics and Kinematics
    REAL(iwp),ALLOCATABLE::eps_acum(:,:)        ! Accumulative particle strain
    REAL(iwp),ALLOCATABLE::m_stress(:,:)        ! Current Particle Stress
    REAL(iwp),ALLOCATABLE::m_stress_ini(:,:)    ! Initial Particle Stress
    REAL(iwp),ALLOCATABLE::m_stress_change(:,:) ! Stress Increment
    
    ! Constitutive Model
    REAL(iwp),ALLOCATABLE::dee(:,:)       ! Elastic Stiffness Matrix
    REAL(iwp),ALLOCATABLE::mean_stress(:) ! Mean Stress (p-q diagram)
    REAL(iwp),ALLOCATABLE::Devstress(:)   ! Deviatoric Stress (p-q diagram)
    REAL(iwp),ALLOCATABLE::epsinvacum(:)  ! Particle Plastic Strain (?)
    REAL(iwp),ALLOCATABLE::mpyield(:)     ! Yield Function Value (F=0:Plastic, F<0:Elastic)
    
    ! Reverse Mapping Particles
    REAL(iwp),ALLOCATABLE::accp(:,:) ! MP present step acceleration
    REAL(iwp),ALLOCATABLE::accb(:,:) ! MP future step acceleration
    REAL(iwp),ALLOCATABLE::vccp(:,:) ! MP present step velocity
    REAL(iwp),ALLOCATABLE::ins(:,:)  ! MP displacement vector
    REAL(iwp),ALLOCATABLE::a_ins(:,:)      ! MP accumulated displacement
    REAL(iwp),ALLOCATABLE::m_mass(:)       ! MP mass
    REAL(iwp),ALLOCATABLE::m_volume(:)     ! Current particle volume
    REAL(iwp),ALLOCATABLE::m_acc(:,:)      ! MP acceleration
    REAL(iwp),ALLOCATABLE::m_velocity(:,:) ! MP velocity
    
    ! Global Body Field
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

    ! Data structure - NEED TO BE REFACTORED TO ALLOW BETTER 3-ARRAY STRUCT
    INTEGER::n_nonzero ! number of non-zero elements in the sparse matrix
    INTEGER,ALLOCATABLE::kdiag(:)    ! skyline matrix diagonal index
    INTEGER,ALLOCATABLE::ia(:)       ! CSR Form: Row Indexes
    INTEGER,ALLOCATABLE::ja(:)       ! CSR Form: Column Indexes
    REAL(iwp),ALLOCATABLE::kv_csr(:) ! CSR Form: Constructed Stiffness
    
    ! Contact-specific variables
    LOGICAL  ::activate=.false.
    REAL(iwp)::phi, gimptol

    CONTAINS

    ! Subroutines and Functions
    PROCEDURE::LOAD_PARTICLES => CLASS_PARTICLE_LOAD_PARTICLES
    PROCEDURE::RESET_MP => CLASS_PARTICLE_RESET_MP
    PROCEDURE::GENERATE_PARTICLE_IN_CELL => CLASS_PARTICLE_GENERATE_PARTICLE_IN_CELL
    PROCEDURE::ALLOCATE_VARIABLES => CLASS_PARTICLE_ALLOCATE_VARIABLES
  END TYPE

  ! Make procedures (subroutines / functions) private
  PRIVATE :: CLASS_PARTICLE_LOAD_PARTICLES, CLASS_PARTICLE_RESET_MP

  CONTAINS

  SUBROUTINE CLASS_PARTICLE_GENERATE_PARTICLE_IN_CELL(this,nx,ny,w1,h1,node,  &
    offsetx,offsety,nst,ndim,nodof)
    !
    ! Generate Quadrilateral FEM Body with the number of integration points 
    ! similar to its nodes per element number. 
    !
    USE GEOMETRY_GENERATOR
    USE FUNCTIONS
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nx,ny,node
    INTEGER,INTENT(IN)::nst,ndim,nodof
    REAL(iwp),INTENT(IN)::w1,h1
    CLASS(particles),INTENT(INOUT)::this
    INTEGER,OPTIONAL,INTENT(IN)::offsetx,offsety
    REAL(iwp),ALLOCATABLE::points(:,:),weights(:),fun(:),m_coord(:,:)
    INTEGER,ALLOCATABLE::MPPE(:),num(:),coord(:,:)
    INTEGER::nmps,i,iel,nip
    INTEGER::offx=0,offy=0
    ! Generate Rectangular FE Mesh
    IF(PRESENT(offsetx)) offx = offsetx
    IF(PRESENT(offsety)) offx = offsety
    CALL RECTANGULAR_2D(this%g_coord,this%g_num,this%nn,this%nels,&
                        nx,ny,w1,h1,node,offx,offy)
    ! Determine number of generated particles
    nmps=this%nels*node
    nip=node
    ! Allocate local variables
    ALLOCATE(MPPE(nip),num(nip),coord(node,ndim),fun(node),&
             points(node,ndim),weights(nip),m_coord(nip,ndim))
    ! Generate Particles from Gaussian Integration Points
    CALL this%ALLOCATE_VARIABLES(ndim,nmps,nst,nodof)
    DO i=1,node
      MPPE(i)=i
    END DO
    CALL SAMPLE_GAUSS_IP('quadrilateral',points,weights)
    DO iel=1,this%nels
      num=this%g_num(:,iel)
      coord=TRANSPOSE(this%g_coord(:,num))
      DO i=1,nip
        CALL SHAPE_FUN(fun,points,i)
        m_coord(i,:)=MATMUL(fun,coord)
      END DO
      this%gm_coord(:,MPPE)=TRANSPOSE(m_coord)
      this%a_ele(MPPE)=iel
      DO i=1,nip
        MPPE(i)=MPPE(i)+nip
      END DO
    END DO
  END SUBROUTINE CLASS_PARTICLE_GENERATE_PARTICLE_IN_CELL


  SUBROUTINE CLASS_PARTICLE_LOAD_PARTICLES(this,index,directory,input_json,    &
    ndim,nodof,nst)
    USE JSON_MODULE
    USE IO
    IMPLICIT NONE
    CLASS(particles),INTENT(INOUT)::this
    INTEGER,INTENT(IN)::index ! particles array index in the input json file
    CHARACTER(*),INTENT(IN)::directory
    TYPE(json_file),INTENT(INOUT)::input_json
    INTEGER,OPTIONAL,INTENT(INOUT)::nodof,nst,ndim
    INTEGER::nmps,def_nodof=2,def_nst=4,def_ndim=2
    CHARACTER(1)::i_char
    CHARACTER(:),ALLOCATABLE::filename,name

    ! get particle generator data
    WRITE(i_char, '(I1)') index
    CALL input_json%GET('particles('//i_char//').generator.name', name)
    CALL input_json%GET('particles('//i_char//').generator.location', filename)
    
    ! get particle locations
    CALL IO_LOAD_PARTICLE(trim(directory), trim(filename), this%gm_coord)
    
    ! Determine Particles Variable
    this%name=name
    nmps = ubound(this%gm_coord, 2)
    if (present(ndim)) def_ndim = ndim
    if (present(nodof)) def_nodof = nodof
    if (present(nst)) def_nst = nst

    CALL this%ALLOCATE_VARIABLES(def_ndim,nmps,def_nst,def_nodof)

    ! Set Material Properties to Zero
    IF(this%nmps>1)THEN 
      this%m_acc           = zero
      this%m_stress        = zero
      this%m_stress_ini    = zero
      this%m_stress_change = zero
      this%accp            = zero
      this%a_ins           = zero
      this%d1x1            = zero
      this%d2x1            = zero
      this%accb            = zero
      this%vccp            = zero
      this%epsinvacum      = zero
      this%Devstress       = zero
      this%ins             = zero
      this%mass            = zero
      this%loads           = zero
      this%gravlo          = zero
      this%ddylds          = zero
      this%m_velocity      = zero
      this%mean_stress     = zero
      this%mpyield         = zero
      this%kv              = zero
      this%mv              = zero
      this%kinup_d1x1      = zero
      this%kinup_d2x1      = zero
    END IF

  END SUBROUTINE CLASS_PARTICLE_LOAD_PARTICLES


  SUBROUTINE CLASS_PARTICLE_ALLOCATE_VARIABLES(this,ndim,nmps,nst,nodof)
    !
    ! Allocate Particle Material Properties
    !
    CLASS(particles),INTENT(INOUT)::this
    INTEGER,INTENT(IN)::ndim,nmps,nst,nodof
    this%nmps=nmps
    this%ndim=ndim
    this%nst=nst
    this%nodof=nodof
    ALLOCATE(this%lm_coord(ndim,nmps))      
    ALLOCATE(this%mweights(nmps))           
    ALLOCATE(this%m_mass(nmps))             
    ALLOCATE(this%m_volume(nmps))           
    ALLOCATE(this%m_stress(nst,nmps))       
    ALLOCATE(this%m_velocity(nodof,nmps))   
    ALLOCATE(this%a_ins(ndim,nmps))         
    ALLOCATE(this%a_ele(nmps))              
    ALLOCATE(this%vccp(ndim,nmps))          
    ALLOCATE(this%accp(ndim,nmps))          
    ALLOCATE(this%epsinvacum(nmps))         
    ALLOCATE(this%Devstress(nmps))          
    ALLOCATE(this%ins(ndim,nmps))           
    ALLOCATE(this%eps_acum(nst,nmps))       
    ALLOCATE(this%support_nodes(9,nmps))    
    ALLOCATE(this%m_stress_ini(nst,nmps))   
    ALLOCATE(this%m_stress_change(nst,nmps))
    ALLOCATE(this%m_acc(nodof,nmps))        
    ALLOCATE(this%mpyield(nmps))            
    ALLOCATE(this%mean_stress(nmps))        
    ALLOCATE(this%lp(ndim,nmps))            
    ALLOCATE(this%member_elements(nmps,4))  
    ALLOCATE(this%accb(ndim,nmps))          
    ALLOCATE(this%mp_dens(nmps))            
    ALLOCATE(this%mpcp(nmps))               
  END SUBROUTINE CLASS_PARTICLE_ALLOCATE_VARIABLES


  SUBROUTINE CLASS_PARTICLE_RESET_MP(this)
    !
    ! Reset time-stepping variables
    !
    IMPLICIT NONE
    CLASS(particles), INTENT(INOUT) :: this
    this%gravlo            = zero
    this%d1x1              = zero
    this%d2x1              = zero
    this%support_nodes     = zero
    this%x1                = zero
    this%loads             = zero
    this%x1_ini            = zero
    this%x1_change         = zero
    this%kinup_Ground_d2x1 = zero
  END SUBROUTINE CLASS_PARTICLE_RESET_MP


END MODULE CLASS_PARTICLE