MODULE CLASS_PARTICLE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp

  TYPE::mpm_particles
    INTEGER :: id
    CHARACTER(16):: name
    INTEGER   :: A,Bi,emps,newnodes,nyp,newel,nmps,nn,slopeopt,                 &
                 slopeel,slope1,tep
              
    REAL(iwp) :: w1,s1,h1
    INTEGER   :: nex,ney,nels,locx,locy,ale,neq,np_types,nprops,dist_x,dist_y,  &
                 nx1,nx2,ny1,ny2

    INTEGER,ALLOCATABLE:: b(:),etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),       &
      num(:),MPPE(:),neighb(:,:),nodecont(:),newels(:),valuesg(:),boundnod(:),  &
      kdiag(:),kconst(:),nodecont_f(:),tied_nn(:,:),base_nn(:)
    
    ! Variables to track material points
    INTEGER,ALLOCATABLE:: a_ele(:),c_ele(:),d_ele(:),b_ele(:),dp_ele(:),        &
      flag(:),GIMP_nodes(:,:),GIMP_node_mp(:,:),k_ele(:),elemmpoints(:,:)
    
    ! Material properties
    REAL(iwp),ALLOCATABLE:: dee(:,:),deeinv(:,:),epsinvacum(:),g_matrix(:),     &
      mpcp(:),prop(:,:)
    
    ! Variables for each material point
    REAL(iwp),ALLOCATABLE:: accp(:,:),a_ins(:,:),Devstress(:),eps_acum(:,:),    &
      gm_coord(:,:),g_coord(:,:),ini_density(:),ini_volume(:),ins(:,:),         &
      mweights(:),m_volume(:),m_coord(:,:),m_stress(:,:), m_velocity(:,:),      &
      m_acc(:,:),mpoints(:,:),m_mass(:),m_stress_ini(:,:),accb(:,:),            &  
      m_stress_change(:,:),vccp(:,:),mpyield(:),mean_stress(:),lp_mp(:,:),      &
      m_stress_prev(:,:),ground_acc(:),mp_dens(:) 
    
    ! Single body field
    REAL(iwp),ALLOCATABLE:: a_field(:,:),ddylds(:),diag(:),d1x1(:),d2x1(:),     &
      eps_m(:,:),eps_1(:,:),eps_2(:,:),fnorm(:),fcont(:),fdamp(:),Freact(:),    &
      f_fint(:),gravlo(:),kv(:),kp(:),kp_2(:),kinup_d2x1(:),kinup_d1x1(:),      &
      loads(:),mv(:),m_mvp(:),m_mva(:),m_field(:),m_phi(:),                     &
      normal(:,:),temp_d1x1(:),temp_d2x1(:),tangent(:,:),                       &
      v_field(:,:),vcm(:),vel_change(:),x1(:),                                  &
      x1_orig(:),x1_change(:),x1_ini(:),f_earth(:),kinup_Ground_d2x1(:),        &
      cdamp(:),mvkv(:)

    ! Contact-specific variables
    REAL(iwp) :: phi, gimptol
    LOGICAL   :: activate=.false.

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
      this%m_volume(this%nmps),                 &
      this%mweights(this%nmps),                 &
      this%m_mass(this%nmps),                   &
      this%m_stress(def_nst, this%nmps),        &
      this%m_velocity(def_nodof, this%nmps),    &
      this%a_ins(def_ndim, this%nmps),          &
      this%flag(this%nmps),                     &
      this%b(this%nmps),                        &
      this%mpoints(this%nmps, def_ndim),        &
      this%a_ele(this%nmps),                    &
      this%accp(def_ndim, this%nmps),           &
      this%vccp(def_ndim, this%nmps),           &
      this%ini_density(this%nmps),              &
      this%epsinvacum(this%nmps),               &
      this%Devstress(this%nmps),                &
      this%ins(def_ndim, this%nmps),            &
      this%eps_acum(def_nst, this%nmps),        &
      this%GIMP_nodes(9, this%nmps),            &
      this%valuesg(this%nmps),                  &
      this%m_stress_ini(def_nst, this%nmps),    &
      this%m_stress_change(def_nst, this%nmps), &
      this%m_acc(def_nodof, this%nmps),         &
      this%mpyield(this%nmps),                  &
      this%mean_stress(this%nmps),              &
      this%lp_mp(def_ndim, this%nmps),          &
      this%elemmpoints(this%nmps, 4),           &
      this%accb(def_ndim, this%nmps),           &
      this%eps_m(def_nst, this%nmps),           &
      this%eps_1(def_nst, this%nmps),           &
      this%eps_2(def_nst, this%nmps),           &
      this%m_stress_prev(def_nst, this%nmps),   &
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
      this%flag            = 0
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
      this%m_stress_prev   = zero
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
    this%GIMP_nodes = zero
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