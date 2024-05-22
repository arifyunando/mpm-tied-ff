MODULE CLASS_PARTICLE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp

  TYPE::mpm_particles
    INTEGER :: id
    CHARACTER(16):: name
    
    ! Constants
    INTEGER:: nmps ! number of material points

    ! Commons Variable
    REAL(iwp),ALLOCATABLE:: gm_coord(:,:) ! MP global coordinates
    INTEGER,ALLOCATABLE:: b(:) 
    INTEGER,ALLOCATABLE:: nf(:,:)
    INTEGER,ALLOCATABLE:: MPPE(:)
    INTEGER,ALLOCATABLE:: neighb(:,:)
    INTEGER,ALLOCATABLE:: nodecont(:)
    INTEGER,ALLOCATABLE:: newels(:)
    INTEGER,ALLOCATABLE:: valuesg(:)
    INTEGER,ALLOCATABLE:: boundnod(:)
    INTEGER,ALLOCATABLE:: kdiag(:)
    INTEGER,ALLOCATABLE:: kconst(:)
    INTEGER,ALLOCATABLE:: nodecont_f(:)
    INTEGER,ALLOCATABLE:: tied_nn(:,:)
    INTEGER,ALLOCATABLE:: base_nn(:)
    
    ! Variables to track material points
    INTEGER:: n_active_ele
    INTEGER,ALLOCATABLE:: a_ele(:) ! Element index each MP resides 
    INTEGER,ALLOCATABLE:: c_ele(:) ! Number of MPs inside a given element
    INTEGER,ALLOCATABLE:: d_ele(:) ! Active element (1:active; 0:deactive)
    INTEGER,ALLOCATABLE:: k_ele(:) ! Accumulation of c_ele
    INTEGER,ALLOCATABLE:: support_nodes(:,:) ! List of Node ids inside mp support domain
    INTEGER,ALLOCATABLE:: member_elements(:,:) ! list of Elements ids inside mp support domain
    
    ! Active Elements Variable
    REAL(iwp),ALLOCATABLE:: g_coord(:,:) ! Mesh node global coordinates
    INTEGER,ALLOCATABLE:: g_num(:,:) ! global element node ids 
    INTEGER,ALLOCATABLE:: num(:)     ! Mesh node indexes
    INTEGER,ALLOCATABLE:: g_g(:,:)   ! global steering factor
    INTEGER,ALLOCATABLE:: g(:)       ! local steering factor

    ! Variables to track active MPM elements
  
    ! Numerical Integration
    REAL(iwp),ALLOCATABLE:: mweights(:)   ! Initial Particle Weights (Initial Particle Volume)
    REAL(iwp),ALLOCATABLE:: lm_coord(:,:) ! Local Particle Coordinates
    REAL(iwp),ALLOCATABLE:: lp(:,:)       ! Local Particle Size (GIMP)

    ! Particle Properties (Print to VTK)
    REAL(iwp),ALLOCATABLE:: a_ins(:,:)    ! Particle Displacement
    REAL(iwp),ALLOCATABLE:: epsinvacum(:) ! Particle Plastic Strain (?)
    REAL(iwp),ALLOCATABLE:: m_mass(:)     ! Particle Mass
    REAL(iwp),ALLOCATABLE:: m_acc(:,:)    ! Particle Acceleration
    REAL(iwp),ALLOCATABLE:: m_velocity(:,:) ! Particle Velocity

    REAL(iwp),ALLOCATABLE:: m_stress_ini(:,:) ! Initial Particle Stress
    REAL(iwp),ALLOCATABLE:: m_stress(:,:)     ! Current Particle Stress
    REAL(iwp),ALLOCATABLE:: m_stress_change(:,:) ! Stress Increment
    REAL(iwp),ALLOCATABLE:: mean_stress(:)    ! Mean Stress (p-q diagram)
    REAL(iwp),ALLOCATABLE:: Devstress(:)      ! Deviatoric Stress (p-q diagram)
    REAL(iwp),ALLOCATABLE:: mpyield(:)        ! Yield Criterion Function Value

    ! Particle Material Properties
    REAL(iwp),ALLOCATABLE:: mp_dens(:) ! Particle densitiy
    REAL(iwp),ALLOCATABLE:: mpcp(:)    ! Particle cohesion

    ! Kinetics and Kinematics
    REAL(iwp),ALLOCATABLE:: dee(:,:)    ! Elastic Stiffness Matrix

    REAL(iwp),ALLOCATABLE:: g_matrix(:) ! Gravitational loading

    REAL(iwp),ALLOCATABLE:: eps_acum(:,:)   ! Accumulative particle strain
    REAL(iwp),ALLOCATABLE:: ini_density(:)  ! Initial density
    REAL(iwp),ALLOCATABLE:: ini_volume(:)   ! Initial particle volume
    REAL(iwp),ALLOCATABLE:: m_volume(:)     ! Current particle volume
    REAL(iwp),ALLOCATABLE:: ground_acc(:)   ! Ground Acceleration
    REAL(iwp),ALLOCATABLE:: ins(:,:)        ! 
    
    ! Single body field
    REAL(iwp),ALLOCATABLE:: accp(:,:)
    REAL(iwp),ALLOCATABLE:: accb(:,:)
    REAL(iwp),ALLOCATABLE:: vccp(:,:)
    REAL(iwp),ALLOCATABLE:: a_field(:,:)

    ! Body force
    REAL(iwp),ALLOCATABLE:: ddylds(:)
    REAL(iwp),ALLOCATABLE:: diag(:)
    REAL(iwp),ALLOCATABLE:: d1x1(:)
    REAL(iwp),ALLOCATABLE:: d2x1(:)
    REAL(iwp),ALLOCATABLE:: eps_m(:,:)
    REAL(iwp),ALLOCATABLE:: eps_1(:,:)
    REAL(iwp),ALLOCATABLE:: eps_2(:,:)
    REAL(iwp),ALLOCATABLE:: fnorm(:)
    REAL(iwp),ALLOCATABLE:: fcont(:)
    REAL(iwp),ALLOCATABLE:: fdamp(:)
    REAL(iwp),ALLOCATABLE:: Freact(:)
    REAL(iwp),ALLOCATABLE:: f_fint(:)
    REAL(iwp),ALLOCATABLE:: gravlo(:)
    REAL(iwp),ALLOCATABLE:: kv(:)
    REAL(iwp),ALLOCATABLE:: kp(:)
    REAL(iwp),ALLOCATABLE:: kp_2(:)
    REAL(iwp),ALLOCATABLE:: kinup_d2x1(:)
    REAL(iwp),ALLOCATABLE:: kinup_d1x1(:)
    REAL(iwp),ALLOCATABLE:: loads(:)
    REAL(iwp),ALLOCATABLE:: mv(:)
    REAL(iwp),ALLOCATABLE:: mvkv(:)
    REAL(iwp),ALLOCATABLE:: m_mvp(:)
    REAL(iwp),ALLOCATABLE:: m_mva(:)
    REAL(iwp),ALLOCATABLE:: m_field(:)
    REAL(iwp),ALLOCATABLE:: m_phi(:)
    REAL(iwp),ALLOCATABLE:: normal(:,:)
    REAL(iwp),ALLOCATABLE:: temp_d1x1(:)
    REAL(iwp),ALLOCATABLE:: temp_d2x1(:)
    REAL(iwp),ALLOCATABLE:: tangent(:,:)
    REAL(iwp),ALLOCATABLE:: v_field(:,:)
    REAL(iwp),ALLOCATABLE:: vcm(:)
    REAL(iwp),ALLOCATABLE:: vel_change(:)
    REAL(iwp),ALLOCATABLE:: x1(:)
    REAL(iwp),ALLOCATABLE:: x1_orig(:)
    REAL(iwp),ALLOCATABLE:: x1_change(:)
    REAL(iwp),ALLOCATABLE:: x1_ini(:)
    REAL(iwp),ALLOCATABLE:: f_earth(:)
    REAL(iwp),ALLOCATABLE:: kinup_Ground_d2x1(:)
    REAL(iwp),ALLOCATABLE:: cdamp(:)

    ! Contact-specific variables
    LOGICAL   :: activate=.false.
    REAL(iwp) :: phi, gimptol

    CONTAINS

    ! Subroutines and Functions
    PROCEDURE :: LOAD_PARTICLES => CLASS_PARTICLE_LOAD_PARTICLES
    PROCEDURE :: RESET_MP => CLASS_PARTICLE_RESET_MP
  END TYPE

  ! Make procedures (subroutines / functions) private
  PRIVATE :: CLASS_PARTICLE_LOAD_PARTICLES, CLASS_PARTICLE_RESET_MP

  CONTAINS

  SUBROUTINE CLASS_PARTICLE_LOAD_PARTICLES(this, index, directory, input_json,  &
    ndim, nodof, nst)
    USE JSON_MODULE
    USE IO
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: index
    CHARACTER(*), INTENT(IN)       :: directory
    CLASS(mpm_particles), INTENT(INOUT) :: this
    TYPE(json_file), INTENT(INOUT) :: input_json
    INTEGER, OPTIONAL, INTENT(INOUT) :: nodof, nst, ndim
    INTEGER :: def_nodof=2, def_nst=4, def_ndim=2
    CHARACTER(1) :: i_char
    CHARACTER(:), ALLOCATABLE :: filename, name

    ! get particle generator data
    WRITE(i_char, '(I1)') index
    CALL input_json%GET('particles('//i_char//').generator.name', name)
    CALL input_json%GET('particles('//i_char//').generator.location', filename)
    this%name = name

    ! get particle locations
    CALL IO_LOAD_PARTICLE(trim(directory), trim(filename), this%gm_coord)
    
    ! Determine Particles Variable
    this%nmps = ubound(this%gm_coord, 2)
    if (present(ndim)) def_ndim = nst
    if (present(nodof)) def_nodof = nst
    if (present(nst)) def_nst = nst

    ! Allocate Particle Material Properties
    ALLOCATE(                                   &
      this%lm_coord(def_ndim, this%nmps),       &
      this%b(this%nmps),                        &
      this%mweights(this%nmps),                 &
      this%m_mass(this%nmps),                   &
      this%m_volume(this%nmps),                 &
      this%m_stress(def_nst, this%nmps),        &
      this%m_velocity(def_nodof, this%nmps),    &
      this%a_ins(def_ndim, this%nmps),          &
      this%a_ele(this%nmps),                    &
      this%accp(def_ndim, this%nmps),           &
      this%vccp(def_ndim, this%nmps),           &
      this%ini_density(this%nmps),              &
      this%epsinvacum(this%nmps),               &
      this%Devstress(this%nmps),                &
      this%ins(def_ndim, this%nmps),            &
      this%eps_acum(def_nst, this%nmps),        &
      this%support_nodes(9, this%nmps),         &
      this%valuesg(this%nmps),                  &
      this%m_stress_ini(def_nst, this%nmps),    &
      this%m_stress_change(def_nst, this%nmps), &
      this%m_acc(def_nodof, this%nmps),         &
      this%mpyield(this%nmps),                  &
      this%mean_stress(this%nmps),              &
      this%lp(def_ndim, this%nmps),             &
      this%member_elements(this%nmps, 4),       &
      this%accb(def_ndim, this%nmps),           &
      this%eps_m(def_nst, this%nmps),           &
      this%eps_1(def_nst, this%nmps),           &
      this%eps_2(def_nst, this%nmps),           &
      this%mp_dens(this%nmps),                  &
      this%mpcp(this%nmps)                      &
    )

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
      this%diag            = zero
      this%loads           = zero
      this%gravlo          = zero
      this%ddylds          = zero
      this%m_velocity      = zero
      this%mean_stress     = zero
      this%mpyield         = zero
      this%kv              = zero
      this%mv              = zero
      this%x1_orig         = zero
      this%nodecont_f      = zero
      this%kinup_d1x1      = zero
      this%kinup_d2x1      = zero
    END IF

  END SUBROUTINE CLASS_PARTICLE_LOAD_PARTICLES


  SUBROUTINE CLASS_PARTICLE_RESET_MP(this)
    !
    ! Reset variables for a new time step
    !
    IMPLICIT NONE
    CLASS(mpm_particles), INTENT(INOUT) :: this
    this%gravlo     = zero
    this%m_mvp      = zero
    this%m_mva      = zero
    this%diag       = zero
    this%fcont      = zero
    this%fnorm      = zero
    this%d1x1       = zero
    this%d2x1       = zero
    this%temp_d1x1  = zero
    this%temp_d2x1  = zero
    this%valuesg    = zero
    this%support_nodes = zero
    this%nodecont   = zero
    this%nodecont_f = zero
    this%normal     = zero
    this%x1         = zero
    this%f_fint     = zero
    this%vcm        = zero
    this%eps_m      = zero
    this%eps_1      = zero
    this%eps_2      = zero
    this%v_field    = zero
    this%m_field    = zero
    this%a_field    = zero
    this%vel_change = zero
    this%loads      = zero
    this%x1_ini     = zero
    this%x1_change  = zero
    this%cdamp      = zero
    this%f_earth    = zero
    this%kinup_Ground_d2x1 = zero
  END SUBROUTINE CLASS_PARTICLE_RESET_MP


END MODULE CLASS_PARTICLE