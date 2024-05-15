MODULE MPMCore
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  
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
    
    !Variables to track material points
    INTEGER,ALLOCATABLE:: a_ele(:),c_ele(:),d_ele(:),b_ele(:),dp_ele(:),        &
      flag(:),GIMP_nodes(:,:),GIMP_node_mp(:,:),k_ele(:),elemmpoints(:,:)
    
    !Material properties
    REAL(iwp),ALLOCATABLE:: dee(:,:),deeinv(:,:),epsinvacum(:),g_matrix(:),     &
      mpcp(:),prop(:,:)
    
    !Variables for each material point
    REAL(iwp),ALLOCATABLE:: accp(:,:),a_ins(:,:),Devstress(:),eps_acum(:,:),    &
      gm_coord(:,:),g_coord(:,:),ini_density(:),ini_volume(:),ins(:,:),         &
      mweights(:),m_volume(:),m_coord(:,:),m_stress(:,:), m_velocity(:,:),      &
      m_acc(:,:),mpoints(:,:),m_mass(:),m_stress_ini(:,:),accb(:,:),            &  
      m_stress_change(:,:),vccp(:,:),mpyield(:),mean_stress(:),lp_mp(:,:),      &
      m_stress_prev(:,:),ground_acc(:),mp_dens(:) 
    
    !Single body field
    REAL(iwp),ALLOCATABLE:: a_field(:,:),ddylds(:),diag(:),d1x1(:),d2x1(:),     &
      eps_m(:,:),eps_1(:,:),eps_2(:,:),fnorm(:),fcont(:),fdamp(:),Freact(:),    &
      f_fint(:),gravlo(:),kv(:),kp(:),kp_2(:),kinup_d2x1(:),kinup_d1x1(:),      &
      loads(:),mv(:),m_mvp(:),m_mva(:),m_field(:),m_phi(:),                     &
      normal(:,:),temp_d1x1(:),temp_d2x1(:),tangent(:,:),                       &
      v_field(:,:),vcm(:),vel_change(:),x1(:),                                  &
      x1_orig(:),x1_change(:),x1_ini(:),f_earth(:),kinup_Ground_d2x1(:),        &
      cdamp(:),mvkv(:)

    !Contact-specific variables
    REAL(iwp) :: phi, gimptol
    LOGICAL   :: activate=.false.
  END TYPE

  CONTAINS
  
END MODULE MPMCore