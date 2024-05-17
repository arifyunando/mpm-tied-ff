MODULE MPM_CORE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp
  
  TYPE::mpm_body
    CHARACTER(10):: name
    INTEGER   :: A,Bi,emps,newnodes,nyp,newel,nmps,nn,slopeopt,                 &
                 slopeel,slope1,tep
              
    REAL(iwp) :: w1,s1,h1
    REAL(iwp) :: Young,Poiss,frictfact
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
  END TYPE

  TYPE::mpm_grid
    ! Constants
    INTEGER :: nnodes, nels

    ! Tracking variables
    REAL(iwp), ALLOCATABLE :: n_coords(:,:)   ! Node Coordinates
    INTEGER, ALLOCATABLE :: el_codes(:,:)     ! Node Index that build an element
    INTEGER, ALLOCATABLE :: nf(:,:)           ! Track degree of freedoms
  END TYPE

  CONTAINS

  SUBROUTINE MPMCORE_LOAD_MESH(directory, filename, grid, nodof, nst)
    USE MPM_IO
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: directory, filename
    TYPE(mpm_grid), INTENT(INOUT) :: grid
    INTEGER, OPTIONAL, INTENT(INOUT) :: nodof, nst
    INTEGER :: def_nodof=2, def_nst=4

    ! Get mesh nodal locations
    CALL LOAD_MESH_DATA(trim(directory), trim(filename), grid%n_coords, grid%el_codes)
    if (present(nodof)) def_nodof = 2
    if (present(nst)) def_nst = 4

    ! Determine Mesh Variable
    grid%nnodes = ubound(grid%n_coords, 2)
    grid%nels = ubound(grid%el_codes, 2)

    ALLOCATE(grid%nf(def_nodof, grid%nnodes))

  END SUBROUTINE MPMCORE_LOAD_MESH


  SUBROUTINE MPMCORE_LOAD_PARTICLES(directory, filename, body, ndim, nodof, nst)
    USE MPM_IO
    IMPLICIT NONE
    TYPE(mpm_body), INTENT(INOUT) :: body
    CHARACTER(*), INTENT(IN)      :: directory, filename
    INTEGER, OPTIONAL, INTENT(INOUT) :: nodof, nst, ndim
    INTEGER :: def_nodof=2, def_nst=4, def_ndim=2

    ! get particle locations
    CALL LOAD_PARTICLE_DATA(trim(directory), trim(filename), body%gm_coord)
    
    ! Determine Particles Variable
    body%nmps = ubound(body%gm_coord, 2)
    if (present(ndim)) def_ndim = nst
    if (present(nodof)) def_nodof = nst
    if (present(nst)) def_nst = nst

    ! Allocate Particle Material Properties
    ALLOCATE(                                                                   &
      body%m_volume(body%nmps),                                                 &
      body%mweights(body%nmps),                                                 &
      body%m_mass(body%nmps),                                                   &
      body%m_stress(def_nst, body%nmps),                                        &
      body%m_velocity(def_nodof, body%nmps),                                    &
      body%a_ins(def_ndim, body%nmps),                                          &
      body%flag(body%nmps),                                                     &
      body%b(body%nmps),                                                        &
      body%mpoints(body%nmps, def_ndim),                                        &
      body%a_ele(body%nmps),                                                    &
      body%accp(def_ndim, body%nmps),                                           &
      body%vccp(def_ndim, body%nmps),                                           &
      body%ini_density(body%nmps),                                              &
      body%epsinvacum(body%nmps),                                               &
      body%Devstress(body%nmps),                                                &
      body%ins(def_ndim, body%nmps),                                            &
      body%eps_acum(def_nst, body%nmps),                                        &
      body%GIMP_nodes(9, body%nmps),                                            &
      body%valuesg(body%nmps),                                                  &
      body%m_stress_ini(def_nst, body%nmps),                                    &
      body%m_stress_change(def_nst, body%nmps),                                 &
      body%m_acc(def_nodof, body%nmps),                                         &
      body%mpyield(body%nmps),                                                  &
      body%mean_stress(body%nmps),                                              &
      body%lp_mp(def_ndim, body%nmps),                                          &
      body%elemmpoints(body%nmps, 4),                                           &
      body%accb(def_ndim, body%nmps),                                           &
      body%eps_m(def_nst, body%nmps),                                           &
      body%eps_1(def_nst, body%nmps),                                           &
      body%eps_2(def_nst, body%nmps),                                           &
      body%m_stress_prev(def_nst, body%nmps),                                   &
      body%mp_dens(body%nmps),                                                  &
      body%mpcp(body%nmps)                                                      &
    )

    ! Set Material Properties to Zero
    IF(body%nmps>1)THEN 
      body%m_acc           = zero
      body%m_stress        = zero
      body%m_stress_ini    = zero
      body%m_stress_change = zero
      body%accp            = zero
      body%a_ins           = zero
      body%d1x1            = zero
      body%d2x1            = zero
      body%accb            = zero
      body%vccp            = zero
      body%flag            = 0
      body%epsinvacum      = zero
      body%Devstress       = zero
      body%ins             = zero
      body%diag            = zero
      body%loads           = zero
      body%gravlo          = zero
      body%ddylds          = zero
      body%m_velocity      = zero
      body%mean_stress     = zero
      body%mpyield         = zero
      body%kv              = zero
      body%mv              = zero
      body%x1_orig         = zero
      body%m_stress_prev   = zero
      body%nodecont_f      = zero
      body%kinup_d1x1      = zero
      body%kinup_d2x1      = zero
    END IF

  END SUBROUTINE MPMCORE_LOAD_PARTICLES


  SUBROUTINE MPMCORE_FORM_GLOBAL_NF(input_json, entity_json, grid)
    USE JSON_MODULE
    IMPLICIT NONE
    TYPE(json_file), INTENT(INOUT) :: input_json, entity_json
    TYPE(mpm_grid), INTENT(INOUT) :: grid
    LOGICAL :: found=.true.
    INTEGER :: i, nset_id, dir
    INTEGER, ALLOCATABLE :: nodes(:)
    CHARACTER(1) :: ic

    grid%nf = 1
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
        grid%nf(dir,nodes) = 0
      END IF
      i = i+1
    END DO
  END SUBROUTINE MPMCORE_FORM_GLOBAL_NF


  SUBROUTINE MPMCORE_INSERT_MPS(body, grid)
    IMPLICIT NONE
    TYPE(mpm_body), INTENT(INOUT) :: body
    TYPE(mpm_grid), INTENT(INOUT) :: grid
    ALLOCATE(                                                                   &
      body%c_ele(1:grid%nels),                                                  &
      body%d_ele(1:grid%nels),                                                  &
      body%k_ele(0:grid%nels)                                                   &
    )

  END SUBROUTINE MPMCORE_INSERT_MPS
  
END MODULE MPM_CORE