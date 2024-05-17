MODULE FEM_Core
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  
  TYPE::fem_body
    CHARACTER(10):: name
    INTEGER:: A,Bi,emps,mpart,newnodes,npart,nyp,newel,nmps,nn,slopeopt,        &
    slopeel,tep,skylength,nn_2D,ntot
        
    INTEGER*8           :: pt(64)
    INTEGER,ALLOCATABLE :: iparm(:)
    INTEGER,ALLOCATABLE :: ia(:),ja(:),ia_aux(:),ja_aux(:)
    INTEGER,ALLOCATABLE :: neq
                  
    REAL(iwp)::w1,s1,h1,w2,h2
    REAL(iwp)::Young,Poiss
    INTEGER:: nex,ney,nels,locx,locy,ale,np_types,nprops,dist_x,dist_y,         &
    nx1,nx2,ny1,ny2
    INTEGER,ALLOCATABLE:: b(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),                &
    num(:),MPPE(:),newels(:),kdiag(:),tied_nn(:,:)
    
    !Variables to track material points
    INTEGER,ALLOCATABLE:: a_ele(:),d_ele(:),b_ele(:),flag(:),k_ele(:)
    
    !Material properties
    REAL(iwp),ALLOCATABLE:: dee(:,:),deeinv(:,:),g_matrix(:),prop(:,:),km(:,:)

    !Variables for each material point
    REAL(iwp),ALLOCATABLE:: accp(:,:),a_ins(:,:),Devstress(:),eps_acum(:,:),    &
    gm_coord(:,:),g_coord(:,:),ins(:,:),g_coord_aux(:,:),m_pore(:,:),           &
    m_volume(:),m_coord(:,:),m_stress(:,:), m_velocity(:,:),                    &
    m_acc(:,:),mpoints(:,:),m_mass(:),m_stress_ini(:,:),accb(:,:),vccp(:,:),    &  
    mean_stress(:),ground_loads(:),m_stress_efe(:,:),m_dens(:)
        
    !Single body field
    REAL(iwp),ALLOCATABLE:: c_fact(:),ddylds(:),diag(:),d1x1(:),d2x1(:),        &
    f_fint(:),gravlo(:),kv(:),kp(:),kinup_d2x1(:),kinup_d1x1(:),loads(:),mv(:), &
    vcm(:),x1(:),P_ext(:),cdamp(:),Cc(:,:),mvis(:),kinup_Ground_d2x1(:),        &
    x1_acum(:), P_matrix(:),kv_CSR(:),kv_CSR_aux(:),c_damp(:),c_force(:),       &
    residual(:), loads_ini(:),loads_end(:),c_matrix(:,:),d1p1(:),mf_matrix(:,:),&
    mf_force(:), loads_base(:),KGC(:,:),MMS(:,:),CCS(:,:),MOD_MTX(:,:)
        
    !-Small-strain formulation
    REAL(iwp),ALLOCATABLE::damping(:),fk(:)

    !Contact-specific variables
    REAL(iwp):: gimptol,val_acc
        
    REAL(iwp),ALLOCATABLE::statev(:,:)
        
    !**Free-Field matrices & variables
    INTEGER,ALLOCATABLE::penpos(:),eq_pen(:),eq_pen_v(:),penpos_v(:)
        
    INTEGER::nels_bar,nip_1D=3,ndof_1D=6,nband,bandval
  END TYPE

  CONTAINS

END MODULE FEM_Core