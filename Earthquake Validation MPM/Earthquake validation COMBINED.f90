PROGRAM Implicit_MPM_eartquake
  USE mpm_main
  USE mpm_geom
  USE mpm_mpm
  USE mpm_gimp
  USE fem
  USE sparse_lib
  IMPLICIT NONE
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

  !** Counters and placeholders
  INTEGER:: i,j,k,w,m,n,q,s,iel,iters,nsrf,ale,limit,bod,plc,printval,printcount, &
            colx,rowy,ielloc,bod2,bod3,iel2,nnx,nny,tcel,step,iters2,noddir,    &
            itersstep,accdata

  INTEGER:: nnz,nnzacum,accval,eqno,substep,tstep,   &
            sstep,v_presc,h_presc
  
  LOGICAL :: DEBUG=.false.

  !** Element/mesh properties
  INTEGER:: ndim=2,ndof=8,nels,neq,nip=4,nn,nod=4,                              &
    nodof=2,nprops=7,np_types,nst=4,nx1,nx2,nye,ny1,ny2,nlen,                   &
    row,column,slope1,slope2,newnodes,nxe,nipslo,                               &
    cont1,cont2,nnxe,nnye,eq1,eq2,dist_y,aiters,count,dist_x,elcont,mvval  
    
  REAL(iwp)::h1,h2,s1,w1,w2,lowbound,maxbound,leftbound,rightbound,             &
              cellsize,Gravf,HSURF,waverage,fm,fk,upbound
  
  LOGICAL:: shape=.false.,slope=.false.,equilibrium=.false.,newrapconv1
  
  !*** AS's Thesis
  ! These variables will become placeholders for tracking boundary elements
  INTEGER,ALLOCATABLE:: left_boundary(:), right_boundary(:)  
  INTEGER:: maxint=1000000, offset_x, offset_y
  
  INTEGER:: gravbod1,gravbod2,smethod,slopeopt,ploop
  
  REAL(iwp):: det,dt=1.0e15_iwp,tol,Tini,Tfin,Tstress,Tfint,Tcont,              &
    dtim,ddt,frictfact,damping, Maxfric,Minfric,Fricmod,                        &
    maxacc,phi,pi,Totalload,Tforce,Tnodeloc

  CHARACTER(LEN=15)::element='quadrilateral',argv
  LOGICAL::converged,activate=.false.,stable=.false.

  REAL(iwp),ALLOCATABLE:: nt(:,:),delta(:,:),delta_1(:),g_matrix(:)

  !** CMPM/GIMP variables
  INTEGER:: mid,values,bound,aver,valuesg,fixed_freedoms
  
  REAL(iwp),ALLOCATABLE:: equations(:,:),boundel(:),identidad(:,:),funCMPM(:),  &
    beeCMPM(:,:),eldCMPM(:),beeextend(:,:),funextend2(:,:),funextend(:),        &
    jac_coord(:,:),gp(:,:),scalar_dee(:),dee_scal(:,:),derextend(:,:)
  
  INTEGER,ALLOCATABLE::eldddylds(:),eqmat(:),eqmat1(:),eqmat2(:),neighb(:,:),nod_num(:,:)

  REAL(iwp):: dv(2),dvn,dvt,mu,fric=0.0,theta=1    

  !** Material parameters/constitutive model
  INTEGER:: fail,softval !Softening parameters
  REAL(iwp):: dq1,dq2,dq3,dsbar,f,fmax,phif,lode_theta,sigm,fnew,k0,            &
    fac,dlambda,dh1,dh2,h_trial,r1,r2,Dsigmanorm,Tsigmanorm,epstol=1.0d-16,     &
    stol,qtol,r_error,dt_min=1.0d-04,DTret,ft,h,Ai,epsinv,ps1,ps2,ps3,          &
    ftol=1.0d-06,mpcr,SModulus,Maxplastic,Tret,dampsign, 			                  &
    sum_scalar,grad_c1,grad_c2,grad_psi1,grad_psi2,dQ,a_slo,KT,TLoad,           &
    dFdc,cpeak,psi,psires,pt5=0.5_iwp,penalty=1.0e20_iwp

  !** Constants/multipliers
  REAL(iwp):: d180=180.0_iwp,zero=0.0_iwp,tnph,tnps,one=1.0_iwp,                &
  two=2.0_iwp,d3=3.0_iwp,alpha_2

  !** SYSTEM ARRAYS (FEM arrays for the entire mesh)
  INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),               &
    num(:),MPPE(:),b(:),ground_acc(:)
    !MPPE=material points per element
  REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),dee(:,:),                          &
    devp(:),eld(:),eload(:),eps(:),evp(:),evpt(:),elso(:),                      &
    flow(:,:),fun(:),gravlo(:),g_coord(:,:),loads(:),m1(:,:),                   &
    m2(:,:),m3(:,:),points(:,:),prop(:,:),sigma(:),                             &
    weights(:),jac(:,:),der(:,:),deriv(:,:),x_coords(:),y_coords(:),            &
    mpcp(:),mpoints(:,:),pl(:,:),derCMPM(:,:),Dsigmae1(:),                      &
    gc(:),eldg(:),sigma_1(:),sigma_2(:),vmfl(:),vmfl_t(:,:),                    &
    dvmfl(:),qmat(:,:),caflow(:),trial(:),vmfla(:),vmflq(:),                    &
    dsigma(:),Freact(:),Dsigma1(:),Dsigma2(:),Dsigmae2(:),                      &
    flowg(:),Deps(:),sigma2(:),flowf(:),epsinvacum(:),depst(:),                 &
    Devstress(:),deriv_gauss(:,:),bee_gauss(:,:),gauss_stress(:),               &
    scalar_stress(:),der_gauss(:,:),fun_gauss(:),sigma_trial(:),                &
    eldGIMP(:),sp_coord(:,:),lp_coord(:,:),km_gauss(:,:),mm_gimp(:,:),          &
    km(:,:),vel_change(:),eps_2(:),props(:)

  REAL(iwp),ALLOCATABLE:: mm(:,:),mm_s(:,:),kv_CSR(:),kv_CSR_aux(:),eps_s(:),   &
    km_mf(:,:),mm_mf(:,:),mm_mf_acum(:,:)

  !** Data structures
  TYPE::mpm_body
    CHARACTER(10):: name
    INTEGER:: A,Bi,emps,newnodes,nyp,newel,nmps,nn,slopeopt,                    &
            slopeel,slope1,tep

              
    REAL(iwp)::w1,s1,h1
    REAL(iwp)::Young,Poiss,frictfact
    INTEGER:: nex,ney,nels,locx,locy,ale,neq,np_types,nprops,dist_x,dist_y,     &
              nx1,nx2,ny1,ny2
    INTEGER,ALLOCATABLE:: b(:),etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),       &
    num(:),MPPE(:),neighb(:,:),nodecont(:),newels(:),valuesg(:),boundnod(:),    &
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
      loads(:),mv(:),m_mvp(:),m_mva(:),m_field(:),m_phi(:),normal(:,:),x1(:),   &
      temp_d1x1(:),temp_d2x1(:),tangent(:,:),v_field(:,:),vcm(:),vel_change(:), &
      x1_orig(:),x1_change(:),x1_ini(:),f_earth(:),kinup_Ground_d2x1(:),cdamp(:),mvkv(:)
    !Contact-specific variables
    REAL(iwp):: phi,gimptol
    LOGICAL::activate=.false.

    !Additional Variables
    INTEGER*8           :: pt(64)
    INTEGER,ALLOCATABLE :: iparm(:)
    INTEGER:: mpart,npart,skylength,nn_2D,ntot
    INTEGER,ALLOCATABLE :: ia(:),ja(:),ia_aux(:),ja_aux(:)
    REAL(iwp),ALLOCATABLE:: g_coord_aux(:,:),m_pore(:,:),ground_loads(:),m_stress_efe(:,:),m_dens(:)
    REAL(iwp),ALLOCATABLE:: c_fact(:),P_ext(:),Cc(:,:),mvis(:),km(:,:),         &
      x1_acum(:), P_matrix(:),kv_CSR(:),kv_CSR_aux(:),c_damp(:),c_force(:),     &
      residual(:), loads_ini(:),loads_end(:),c_matrix(:,:),d1p1(:),mf_matrix(:,:),&
      mf_force(:), loads_base(:),KGC(:,:),MMS(:,:),CCS(:,:),MOD_MTX(:,:)
    INTEGER,ALLOCATABLE::penpos(:),eq_pen(:),eq_pen_v(:),penpos_v(:)
    REAL(iwp),ALLOCATABLE::ddsddt(:),drplde(:),drpldt(:),stran(:),              &
      time(:),predef(:),dpred(:),props(:),drot(:,:),dfgrd0(:,:),dfgrd1(:,:),    &
      stressumat(:),epsumat(:),deeumat(:,:),statev(:,:)
    REAL(iwp),ALLOCATABLE::damping(:),fk(:)
  END TYPE
  
  INTEGER:: maxfct, mnum, mtype, phase, nrhs, error, msglvl
  INTEGER:: nels_bar,nip_1D=3,ndof_1D=6,nband,bandval
  REAL(iwp),ALLOCATABLE::kd(:),mm_acum(:,:),mvis(:),Cc(:,:)
  REAL(iwp) :: dparm(64)
  REAL(iwp) :: ddum,valone
  INTEGER   :: idum(1)
 
  !** Define the MPM objects to be used in the simulation
  TYPE(mpm_body):: mbod(3)

  ! Additional variables
  INTEGER::nyp,emps,nmps,tep,newel,slopeel,A,Bi
  INTEGER::nshr,nstatv=100,npt,layer,kspt,kstep,kinc,npropsum,void
  INTEGER,PARAMETER::ndi=3,nstumat=6

  INTEGER,ALLOCATABLE:: sum_vol(:),c_ele(:),m_num(:),a_ele(:),                  &
    k_ele(:),b_ele(:),d_ele(:),dp_ele(:),flag(:),g_s(:)

  REAL(iwp),ALLOCATABLE::ini_volume(:),gm_coord(:,:),m_coord(:,:),              &
    mweights(:),m_mass(:),m_volume(:),ini_density(:),m_stress(:,:),             &
    m_velocity(:,:),a_ins(:,:),nod_stress(:,:),d1x1(:,:),                       &
    diag(:,:),d2x1(:),m_mvp(:,:),m_mva(:),ddylds(:),                            &
    accp(:,:),acc(:),vcc(:),vccp(:,:),ins(:,:),stress(:),ecm_acum(:,:),ecm(:,:)

  REAL(iwp),ALLOCATABLE::ddsddt(:),drplde(:),drpldt(:),stran(:),              &
    time(:),predef(:),dpred(:),drot(:,:),dfgrd0(:,:),dfgrd1(:,:),      &
    stressumat(:),epsumat(:),deeumat(:,:),statev(:,:)

  
  !=============================================================================
  ! Input and Initialisation                              
  !=============================================================================
    
  nlen=7
  argv='Results'
  OPEN(800,FILE='Output/log.output')
  
  OPEN(10,FILE='Input/Datafound.dat',status='old')
  OPEN(300,FILE='Input/Groundacc.dat',status='old')
  
  ! read material properties
  OPEN(400,FILE='Input/parumat.dat',status='old')
  READ(400,*)npropsum
  ALLOCATE(                                                                   &
      ddsddt(nstumat),drplde(nstumat),stran(nstumat),                         &
      props(npropsum),stressumat(nstumat),epsumat(nstumat),                   &
      deeumat(nstumat,nstumat)                                                &
  )
  READ(400,*)props

  ! Read parameter for the MPM body (bod=1) [AS]
  bod=1          
  mbod(1)%nprops=6 !For Von Mises with softening == 7
  mbod(1)%ntot=ndof
  READ(10,*)mbod(bod)%name,mbod(bod)%slopeopt,mbod(bod)%w1,mbod(bod)%h1,mbod(bod)%s1, &
    mbod(bod)%nex,mbod(bod)%ney,mbod(bod)%dist_x,mbod(bod)%dist_y,mbod(bod)%np_types
  ALLOCATE(mbod(bod)%prop(mbod(bod)%nprops,mbod(bod)%np_types))
  READ(10,*)mbod(bod)%prop

  ! Copy the data from body 1 to freefields [AS]
  DO bod=2,size(mbod)
    mbod(bod)%slopeopt=mbod(1)%slopeopt
    mbod(bod)%w1=1
    mbod(bod)%h1=mbod(1)%h1
    mbod(bod)%s1=mbod(1)%s1
    mbod(bod)%nex=mbod(1)%nex
    mbod(bod)%ney=mbod(1)%ney
    mbod(bod)%dist_x=mbod(1)%dist_x
    mbod(bod)%dist_y=mbod(1)%dist_y
    mbod(bod)%np_types=mbod(1)%np_types
    mbod(bod)%nprops=mbod(1)%nprops
    mbod(bod)%ntot=mbod(1)%ntot

    ALLOCATE(mbod(bod)%prop(mbod(bod)%nprops,mbod(bod)%np_types))
    mbod(bod)%prop=mbod(1)%prop
  END DO

  ! 
  DO bod=1,size(mbod)
    mbod(bod)%nx1=mbod(bod)%nex
    mbod(bod)%ny1=mbod(bod)%ney
    w2=0.0;h2=0.0;mbod(bod)%nx2=0;mbod(bod)%ny2=0
    mbod(bod)%nx2=0;mbod(bod)%ny2=0
  END DO
  
  
  !=============================================================================
  ! Simulation Settings                              
  !=============================================================================

  !Select "smethod" which indicats the use of MPM(smethod=1), GIMP(smethod=2), or CMPM(smethod=3)
  !It is recomended to use or smethod=1 or smethod=3
  smethod=3

  !PRINT*,'Use Gravity in bod 1:Yes 2:No'
  gravbod1=1

  !-The plastic loop jis to consider plasticity. ploop = 1 consider plasticity, any other value does not
  ploop=0
  
  printval = 5
  PRINT*, 'Printing results:', printval
  
  Tini=0.0;Tfin=0.0;Tcont=0.0
  call cpu_time(Tini)
  print*,'Initial time',Tini

  !=============================================================================
  ! Generate body and mesh geometries                             
  !=============================================================================
  
  ! Define main MPM body
  bod=1          
  mbod(bod)%ney=mbod(bod)%ny1+mbod(bod)%ny2
  mbod(bod)%nex=mbod(bod)%nx1+mbod(bod)%nx2
  mbod(bod)%newel=mbod(bod)%ney            !Number of elements in the x direction in the slope section
  mbod(bod)%slopeel=0
  mbod(bod)%newnodes=(mbod(bod)%ney-1)+1
  mbod(bod)%slope1=mbod(bod)%nx1+1
  
  cellsize=mbod(bod)%w1/mbod(bod)%nx1
  mbod(bod)%gimptol=1.0e-5

  ! Define Free-Fields body from inital main MPM body
  DO bod=2,size(mbod) 
    mbod(bod)%ney=mbod(bod)%ny1+mbod(bod)%ny2
    mbod(bod)%nex=2
    mbod(bod)%newel=mbod(bod)%ney 
    mbod(bod)%slopeel=0
    mbod(bod)%newnodes=(mbod(bod)%ney-1)+1
    mbod(bod)%slope1=mbod(bod)%nx1+1
    mbod(bod)%gimptol=1.0e-5
  END DO
  
  !*** AS's THESIS
  ! allocate left and right boundary marker such that every row of elements in free-fields body as 2 correspondence column
  ! in the main mesh
  ALLOCATE(left_boundary(mbod(2)%ney), right_boundary(mbod(2)%ney))
  
  ! The section below calculates the number of elements (nels), number of nodes (nn)
  ! and number of stress/material points (nmps) based on the inputted geometry
  ALLOCATE(mpoints(1,2))
  count_nodes: DO bod=1,size(mbod)
    ! slopeopt=2 indicates that there is no slope
    IF(mbod(bod)%slopeopt==2)THEN
      s1=0
      mbod(bod)%s1=0
      mbod(bod)%newel=0
      mbod(bod)%newnodes=0
      mbod(bod)%slope1=mbod(bod)%nex
    END IF
    
    ! disregarded since geometry has no slope 
    DO i=1,mbod(bod)%newel
      mbod(bod)%slopeel=mbod(bod)%slopeel+i  !Slopeel is the number of elements in the slope!Part of the slope
    END DO
    DO i=1,mbod(bod)%newel-1
      mbod(bod)%newnodes=mbod(bod)%newnodes+(mbod(bod)%ney-i)+1   !number of nodes in the slope section
    END DO
    
    ! calculate the number of elements and number of nodes generated
    IF (mbod(bod)%slopeopt==2)THEN
      mbod(bod)%nels=mbod(bod)%nex*mbod(bod)%ney
    ELSE
      mbod(bod)%nels=mbod(bod)%nex*mbod(bod)%ney+mbod(bod)%slopeel
    END IF
    mbod(bod)%nn=(mbod(bod)%nex+1)*(mbod(bod)%ney+1)+mbod(bod)%newnodes+1
    IF(mbod(bod)%slopeopt==2) mbod(bod)%nn=mbod(bod)%nn-1

    mbod(bod)%emps=nip
    mbod(bod)%nmps=mbod(bod)%nels*nip

    ! Allocate array dimension dictated by nels and nn
    IF(bod>1) mbod(bod)%skylength=500*mbod(bod)%nels ! for freefield bodies
    IF(bod==1)THEN
      ALLOCATE(                                             &
        mbod(bod)%g_coord(ndim,mbod(bod)%nn),               &
        mbod(bod)%g_num(nod,mbod(bod)%nels),                &
        mbod(bod)%g_g(ndof,mbod(bod)%nels),                 &
        mbod(bod)%nf(nodof,mbod(bod)%nn),                   &
        mbod(bod)%deeinv(nst,nst),                          &
        mbod(bod)%dee(nst,nst),                             &
        mbod(bod)%ini_volume(mbod(bod)%nels),               &
        mbod(bod)%c_ele(mbod(bod)%nels))                    
    ELSE
      ! Note: %ntot is defined as ndof
      ALLOCATE(                                             &
        mbod(bod)%g_coord(ndim,mbod(bod)%nn),               &
        mbod(bod)%g_num(nod,mbod(bod)%nels),                &
        mbod(bod)%g_g(mbod(bod)%ntot,mbod(bod)%nels),       &
        mbod(bod)%nf(nodof,mbod(bod)%nn),                   &
        mbod(bod)%deeinv(nst,nst),                          &
        mbod(bod)%dee(nst,nst),                             &
        mbod(bod)%ini_volume(mbod(bod)%nels),               &
        mbod(bod)%c_ele(mbod(bod)%nels))
      
      ! Tied-FF Variables
      ALLOCATE(                                             &
        mbod(bod)%g_coord_aux(ndim,mbod(bod)%nn),           &  
        mbod(bod)%g(mbod(bod)%ntot),                        &  
        mbod(bod)%tied_nn(mbod(bod)%nn,2),                  & 
        mbod(bod)%ia_aux(mbod(bod)%skylength),              &  
        mbod(bod)%ja_aux(mbod(bod)%skylength),              & 
        mbod(bod)%km(mbod(bod)%ntot,mbod(bod)%ntot),        & 
        mbod(bod)%KGC(mbod(bod)%ntot,mbod(bod)%ntot),       & 
        mbod(bod)%MMS(mbod(bod)%ntot,mbod(bod)%ntot),       & 
        mbod(bod)%CCS(mbod(bod)%ntot,mbod(bod)%ntot),       & 
        mbod(bod)%MOD_MTX(mbod(bod)%ntot,mbod(bod)%ntot)    & 
      )
    END IF
    !- IF mbod(bod)%kconst=1, then the stiffness in the full domain will be constant
    !- any other value avoid using a constant stiffness
    ALLOCATE(mbod(bod)%kconst(1))
    mbod(bod)%kconst=0
  END DO count_nodes

  
  !===========================================================================AS
  ! GLOBAL VARIABLE ALLOCATIONS (Based on MPM body and mesh)
  !===========================================================================AS
  nels = mbod(1)%nels
  
  ALLOCATE(                                                                       &
    g(ndof),weights(nip),points(nip,ndim),num(nod),coord(nod,ndim),bee(nst,ndof), &
    eps(nst),eps_2(nst),sigma(nst),sigma2(nst),evp(nst),devp(nst),jac(ndim,ndim), &
    der(ndim,nod),fun(nod),deriv(ndim,nod),                                       &
    flow(nst,nst),stress(nst),MPPE(nip),pl(nst,nst),nod_num(nod,1),               &
    sigma_1(nst),sigma_2(nst),evpt(nst),sigma_trial(nst),                         &
    vmfl(nst),vmfl_t(1,nst),dvmfl(nst),                                           &
    qmat(nst,nst),caflow(nst),trial(nst),                                         &
    vmfla(nst),vmflq(nst),flowf(nst),flowg(nst),                                  &
    dsigma(nst),depst(nst),eld(ndof),                                             &
    Dsigma1(nst),Dsigma2(nst),Deps(nst),                                          &
    Dsigmae1(nst),Dsigmae2(nst),                                                  &
    m1(nst,nst),m2(nst,nst),m3(nst,nst),                                          &
    scalar_dee(nod),dee_scal(nst,nst),                                            &
    deriv_gauss(ndim,nod),bee_gauss(nst,ndof),                                    &
    fun_gauss(nod),der_gauss(ndim,nod),                                           &
    eqmat(nodof),eqmat1(nodof),eqmat2(nodof),                                     &
    sp_coord(ndim,1),lp_coord(1,ndim),                                            &
    km_gauss(ndof,ndof),km(ndof,ndof))
    
  ALLOCATE(m_num(nip),vcc(ndof),m_coord(nip,ndim),ini_volume(nels),gc(ndim),ecm(ndof,ndof),ecm_acum(ndof,ndof))
  ALLOCATE(gp(1,2))
      
  
  !===========================================================================AS
  ! Geometry Generation
  !===========================================================================AS
  
  Mesh_create: DO bod=1,size(mbod)
    row=1; column=1
    A=0; Bi=0; dist_x=0
    
    !-------------------------------------------------------------------------AS
    ! Determin Node Locations and Generate Mesh
    !-------------------------------------------------------------------------AS

    !- Set up the global nf
    mbod(bod)%nf=0
    CALL emb_2d_bc2(mbod(bod)%nex,0,mbod(bod)%ney,mbod(bod)%newel,mbod(bod)%nf)

    !- Loop the elements to determine the global node and element numbering
    DO iel=1,mbod(bod)%nels
      CALL emb_2d_geom2(iel,bod,mbod(bod)%nex,mbod(bod)%ney,mbod(bod)%s1,mbod(bod)%newel,       &
        mbod(bod)%slopeel,row,column,mbod(bod)%slope1,mbod(bod)%w1,mbod(bod)%h1,coord,num,A,Bi, &
        mbod(bod)%dist_x,mbod(bod)%dist_y,mbod(bod)%slopeopt)

      CALL num_to_g(num,mbod(bod)%nf,g)
      mbod(bod)%g_num(:,iel)=num
      mbod(bod)%g_coord(:,num)=TRANSPOSE(coord)
      mbod(bod)%g_g(:,iel)=g
    END DO

    !- Compute the initial volume of the elements
    Initialvol: DO iel=1,mbod(bod)%nels
      num=mbod(bod)%g_num(:,iel)
      coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
      mbod(bod)%ini_volume(iel)=cellsize*cellsize
    END DO Initialvol

    !- Allocate material points state variables
    IF(bod==1)THEN
      ALLOCATE(mbod(bod)%gm_coord(ndim,mbod(bod)%nmps),mbod(bod)%m_volume(mbod(bod)%nmps),              &
        mbod(bod)%mweights(mbod(bod)%nmps),mbod(bod)%m_mass(mbod(bod)%nmps),                            &
        mbod(bod)%m_stress(nst,mbod(bod)%nmps),mbod(bod)%m_velocity(nodof,mbod(bod)%nmps),              &
        mbod(bod)%a_ins(ndim,mbod(bod)%nmps),mbod(bod)%flag(mbod(bod)%nmps),                            &
        mbod(bod)%mpoints(mbod(bod)%nmps,ndim),mbod(bod)%a_ele(mbod(bod)%nmps),                         &
        mbod(bod)%accp(ndim,mbod(bod)%nmps),mbod(bod)%vccp(ndim,mbod(bod)%nmps),                        &
        mbod(bod)%ini_density(mbod(bod)%nmps),mbod(bod)%epsinvacum(mbod(bod)%nmps),                     &
        mbod(bod)%Devstress(mbod(bod)%nmps),mbod(bod)%ins(ndim,mbod(bod)%nmps),                         &
        mbod(bod)%eps_acum(nst,mbod(bod)%nmps),mbod(bod)%GIMP_nodes(9,mbod(bod)%nmps),                  &
        mbod(bod)%valuesg(mbod(bod)%nmps),mbod(bod)%m_stress_ini(nst,mbod(bod)%nmps),                   &
        mbod(bod)%m_stress_change(nst,mbod(bod)%nmps),mbod(bod)%m_acc(nodof,mbod(bod)%nmps),            &
        mbod(bod)%mpyield(mbod(bod)%nmps),mbod(bod)%b(mbod(bod)%nmps),mbod(bod)%m_dens(mbod(bod)%nmps), &
        mbod(bod)%mean_stress(mbod(bod)%nmps),mbod(bod)%lp_mp(ndim,mbod(bod)%nmps),                     &
        mbod(bod)%elemmpoints(mbod(bod)%nmps,4),mbod(bod)%accb(ndim,mbod(bod)%nmps),                    &
        mbod(bod)%eps_m(nst,mbod(bod)%nmps),mbod(bod)%eps_1(nst,mbod(bod)%nmps),mbod(bod)%eps_2(nst,mbod(bod)%nmps),&
        mbod(bod)%m_stress_prev(nst,mbod(bod)%nmps),mbod(bod)%mp_dens(mbod(bod)%nmps))
    ELSE
      ALLOCATE(                                                             &
        mbod(bod)%gm_coord(ndim,mbod(bod)%nmps),                            &
        mbod(bod)%m_volume(mbod(bod)%nmps),                                 &
        mbod(bod)%m_mass(mbod(bod)%nmps),                                   &
        mbod(bod)%mweights(mbod(bod)%nmps),                                 &
        mbod(bod)%statev(mbod(bod)%nmps,100),                               &
        mbod(bod)%m_stress(nst,mbod(bod)%nmps),                             &
        mbod(bod)%m_velocity(nodof,mbod(bod)%nmps),                         &
        mbod(bod)%a_ins(ndim,mbod(bod)%nmps),                               &
        mbod(bod)%flag(mbod(bod)%nmps),                                     &
        mbod(bod)%b(mbod(bod)%nmps),                                        &
        mbod(bod)%mpoints(mbod(bod)%nmps,ndim),                             &
        mbod(bod)%a_ele(mbod(bod)%nmps),                                    &
        mbod(bod)%accp(ndim,mbod(bod)%nmps),                                &
        mbod(bod)%vccp(ndim,mbod(bod)%nmps),                                &
        mbod(bod)%ini_density(mbod(bod)%nmps),                              &
        mbod(bod)%m_acc(nodof,mbod(bod)%nmps),                              &
        mbod(bod)%m_dens(mbod(bod)%nmps),                                   &
        mbod(bod)%Devstress(mbod(bod)%nmps),                                &
        mbod(bod)%ins(ndim,mbod(bod)%nmps),                                 &
        mbod(bod)%eps_acum(nst,mbod(bod)%nmps),                             &
        mbod(bod)%m_stress_ini(nst,mbod(bod)%nmps),                         &
        mbod(bod)%mean_stress(mbod(bod)%nmps),                              &
        mbod(bod)%accb(ndim,mbod(bod)%nmps),                                &
        mbod(bod)%Cc(2,2),                                                  &
        mbod(bod)%m_pore(1,mbod(bod)%nmps),                                 &
        mbod(bod)%m_stress_efe(nst,mbod(bod)%nmps),                         &
        mbod(bod)%damping(mbod(bod)%nmps),                                  &
        mbod(bod)%fk(mbod(bod)%nmps),                                       &
        mbod(bod)%mp_dens(mbod(bod)%nmps)                                   &
      )
    END IF
    ALLOCATE(mbod(bod)%mpcp(mbod(bod)%nmps))
    
    
    !-------------------------------------------------------------------------AS
    ! Insert Material Points / Stress Points
    !-------------------------------------------------------------------------AS
  
    CALL sample2(element,points,weights)  
    
    ALLOCATE(mbod(bod)%MPPE(nip))
    mbod(bod)%a_ele=0
    mbod(bod)%c_ele=0
    mbod(bod)%MPPE=0
    DO i=1,nip
      mbod(bod)%MPPE(i)=i
    END DO

    InsertMP: DO iel=1,mbod(bod)%nels
      num=mbod(bod)%g_num(:,iel)
      coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
      DO i=1,nip
        CALL shape_fun(fun,points,i)
        m_coord(i,:)=MATMUL(fun,coord)
      END DO
      m_num=mbod(bod)%MPPE
      mbod(bod)%gm_coord(:,m_num)=TRANSPOSE(m_coord)
      mbod(bod)%mweights(m_num)=weights
      mbod(bod)%a_ele(m_num)=iel
      mbod(bod)%c_ele(iel)=mbod(bod)%emps
      DO j=1,nip
        mbod(bod)%MPPE(j)=mbod(bod)%MPPE(j)+nip
      END DO
    END DO InsertMP

    !-------------------------------------------------------------------------AS
    ! Calculate initial volume and mass of each MPs (not elements)
    !-------------------------------------------------------------------------AS
    DO i=1,mbod(bod)%nmps
      iel=mbod(bod)%a_ele(i)
      mbod(bod)%ini_density(i)=mbod(bod)%prop(1,1)
      mbod(bod)%mp_dens(i)=(mbod(bod)%prop(1,1)+0.75)/(1.75) !-Note that this equation consider a void ratio of 0.75 (Desn= (Gs+e)/(1+e) ), this must be changed 
      mbod(bod)%m_dens(i)=(mbod(bod)%prop(1,1)+0.75)/(1.75) !-Note that this equation consider a void ratio of 0.75 (Desn= (Gs+e)/(1+e) ), this must be changed. Also the freefield use m_dens instead mp_dens
      mbod(bod)%m_volume(i)=mbod(bod)%ini_volume(iel)/mbod(bod)%emps ! for FEM this should be based on the gaussian weight. [maybe todo]
      mbod(bod)%m_mass(i)=mbod(bod)%mp_dens(i)*mbod(bod)%m_volume(i)
    END DO
  END DO Mesh_create

  !===================================================================================
  ! ---------------------------- NEW BACKGROUND GRID ---------------------------------
  ! --The new grid must cover all the space considering the deformation of the body---
  !===================================================================================
  
  !- Deallocate MPM bodies
  DEALL_BOD:DO bod=1,1
    DEALLOCATE(mbod(bod)%nf,mbod(bod)%c_ele,mbod(bod)%g_g)
  END DO DEALL_BOD
      
  ! Calculate the span of the background mesh. Since there is only one body, the span is nx1 and ny1
  DO bod=1,1
    dist_x=mbod(bod)%dist_x+mbod(bod)%nx1
    dist_y=mbod(bod)%dist_y+mbod(bod)%ny1
  END DO  

  ! Define the geometry of the background mesh
  offset_y = 1; offset_x=30
  nx1=dist_x + offset_x ! add 30 columns of elements to the right of the domain
  ny1=dist_y + offset_y ! add 1 row of elements at the bottom of the domain
  w1=cellsize*nx1 ! Size of domain in x-dir in units of distance
  h1=cellsize*ny1 ! Size of domain in y-dir in units of distance
  nx2=0; ny2=0; s1=0 ! no additional elements and no slope
  
  ! Calculate number of elements and number of nodes
  nels=nx1*ny1
  nn=(ny1+1)*(nx1+1)
  nnxe=nx1
  nnye=ny1
  nnx=nx1+1
  nny=ny1+1

  !- Reallocate Global Variable
  ALLOCATE(nf(nodof,nn),g_g(ndof,nels),g_coord(ndim,nn),c_ele(nels),    &
    g_num(nod,nels),nod_stress(nst,nn),k_ele(0:nels),d_ele(nels),       &
    etype(nels),nt(ndof,nodof),delta(ndim,ndim),g_matrix(ndim),         &
    delta_1(ndim),acc(ndof),dp_ele(nels),neighb(nels,8))

  ALLOCATE(mm(ndof,ndof),mm_s(ndof,ndof))
    
  !- Allocate MPM body variable bounded to background mesh
  ALL_BOD:DO bod=1,1
    ALLOCATE(mbod(bod)%nf(nodof,nn),mbod(bod)%c_ele(nels),              &
      mbod(bod)%k_ele(0:nels),mbod(bod)%d_ele(nels),                    &
      mbod(bod)%etype(nels),mbod(bod)%nodecont(nn),                     &
      mbod(bod)%dp_ele(nels),mbod(bod)%v_field(0:nn,ndim),              &
      mbod(bod)%a_field(0:nn,ndim),mbod(bod)%m_field(0:nn),             &
      mbod(bod)%g_g(ndof,nels),mbod(bod)%g_matrix(ndim),                &
      mbod(bod)%normal(nodof,nn),mbod(bod)%tangent(nodof,nn),           &
      mbod(bod)%GIMP_node_mp(mbod(bod)%nmps,nn),mbod(bod)%boundnod(nn), &
      mbod(bod)%nodecont_f(nn),mbod(bod)%tied_nn(nny,6),mbod(bod)%base_nn(nn))
  END DO ALL_BOD

  ALLOCATE(m_mvp(0:nn,ndim),diag(0:nn,ndim),d1x1(0:nn,ndim))

  delta=0.0
  delta(1,1)=1.0
  delta(2,2)=1.0
  delta_1=(/1.0,1.0/)

  ! Generate the mesh and put the definitions in g_num, g_coord, and g_g
  CALL emb_2d_bc(nx1,nx2,ny1,ny2,nf)
  elements_2: DO iel=1,nels
    CALL emb_2d_geom(iel,nx1,nx2,nnye,ny2,w1,s1,w2,h1,h2,coord,num)
    CALL num_to_g(num,nf,g)
    g_num(:,iel)=num
    g_coord(:,num)=TRANSPOSE(coord)
    g_g(:,iel)=g
  END DO elements_2

  ! determine neighbouring element of the background mesh
  CALL neigh_b(nx1,ny1,nels,neighb,'x') 

  
  !===========================================================================AS
  ! Locate material point locations
  !===========================================================================AS
  
  ! a_ele(iel) = element id for each material point
  ! c_ele(iel)= total of MP inside each element
  ! k_ele(iel)= total of MP in the domain
  ! b_ele(k) = list of elements with MP inside
  MPM_Flags: DO bod=1,1
    mbod(bod)%d_ele=0
    mbod(bod)%c_ele=0
    !- loop over the material points and locate the mesh through its global coordinate
    DO i=1,mbod(bod)%nmps
      ! calculate the element number column/row-wise coordinate
      colx=mbod(bod)%gm_coord(1,i)/cellsize+1       
      rowy=ABS(mbod(bod)%gm_coord(2,i))/cellsize+1
        
      ! Determine the element number based on the coordinate and list it in a_ele
      ielloc=(rowy-1.0)*nx1+colx
      mbod(bod)%a_ele(i)=ielloc

      ! Locate local position of material point 'i' in element 'ielloc'
      num=g_num(:,ielloc)
      coord=TRANSPOSE(g_coord(:,num))
      sp_coord(:,1)=mbod(bod)%gm_coord(:,i)
      lp_coord(1,:)=mbod(bod)%mpoints(i,:)
      CALL floc(coord,sp_coord,lp_coord,i)
      mbod(bod)%mpoints(i,:)=lp_coord(1,:)
      mbod(bod)%d_ele(ielloc)=1
    END DO
    
    !- count the number of material points inside an element.
    CALL couma(nels,mbod(bod)%nmps,mbod(bod)%a_ele,mbod(bod)%c_ele,mbod(bod)%k_ele,etype)

    !- list the corresponding activated elements
    k=1; mbod(bod)%b=0
    DO i=1,mbod(bod)%nmps
      iel=mbod(bod)%a_ele(i)
      DO j=1,mbod(bod)%c_ele(iel)
        IF(mbod(bod)%b(mbod(bod)%k_ele(iel-1)+j)==0) THEN
          mbod(bod)%b(mbod(bod)%k_ele(iel-1)+j)=i
          EXIT
        END IF
      END DO
    END DO
  END DO MPM_Flags

  ! This flagging process for FEM relies on the stress points are indexed in the order of element number
  ! Non ordered stress points will BREAK the simulation
  FEM_Flags: DO bod=2,size(mbod)
    ielloc=1; m=0
    DO i=1,mbod(bod)%nmps
      m = m+1
      IF(m>nip)THEN 
        ielloc=ielloc+1; m=1 ! goto next element, reset counter
      END IF
      
      ! list element number in a_ele
      mbod(bod)%a_ele(i) = ielloc

      ! get the local coordinates of each integration points
      num           = mbod(bod)%g_num(:,ielloc)
      coord         = TRANSPOSE(mbod(bod)%g_coord(:,num))
      sp_coord(:,1) = mbod(bod)%gm_coord(:,i)
      lp_coord(1,:) = mbod(bod)%mpoints(i,:)
      CALL floc(coord,sp_coord,lp_coord,i)
      mbod(bod)%mpoints(i,:)=lp_coord(1,:)
    END DO
  END DO FEM_Flags

  
  !===========================================================================AS
  ! Setup boundary conditions
  !===========================================================================AS
  
  nf=0; g_g=0
  Body_fixities: DO bod=1,1
    IF(mbod(bod)%nmps>1)THEN 
      mbod(bod)%g_g=0;mbod(bod)%nf=0
      IF(nip==4) mbod(bod)%lp_mp=cellsize/4.0   
      IF(nip==9) mbod(bod)%lp_mp=cellsize/6.0   
      
      CALL GIMP_activenode(g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp, &
        mbod(bod)%gimptol,neighb,mbod(bod)%a_ele,mbod(bod)%nf,mbod(bod)%GIMP_node_mp)
      CALL Elem_suport(mbod(bod)%c_ele,mbod(bod)%nf,g_num,cellsize,nx1,nip,nels,  &
        g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%gimptol,smethod,mbod(bod)%elemmpoints)

      lowbound=minval(g_coord)
      upbound=maxval(g_coord)

      ! Set equation marker for zero-displacement boundary, forming nf, and determine number of equations
      j=1; m=1
      DO i=1,nn
        !- Body fixities  
        IF(g_coord(2,i)<lowbound+cellsize+0.01) mbod(bod)%nf(2,i)=0 
        IF(g_coord(2,i)<lowbound+cellsize+0.01) m=m+1
        !- Global fixities
        IF(g_coord(1,i)<cellsize+0.01_iwp)      nf(1,i)=0
        IF(g_coord(1,i)<0.01_iwp)               nf(:,i)=0
        IF(g_coord(2,i)<lowbound+0.01+cellsize) nf(:,i)=0  
      END DO
      CALL formnf(mbod(bod)%nf)
      mbod(bod)%neq=MAXVAL(mbod(bod)%nf)

      ! Allocate variables that tied to number of equations
      ALLOCATE(                                 &
        mbod(bod)%diag(0:mbod(bod)%neq),        &
        mbod(bod)%kdiag(mbod(bod)%neq),         &
        mbod(bod)%d1x1(0:mbod(bod)%neq),        &
        mbod(bod)%d2x1(0:mbod(bod)%neq),        &
        mbod(bod)%temp_d2x1(0:mbod(bod)%neq),   &
        mbod(bod)%temp_d1x1(0:mbod(bod)%neq),   &
        mbod(bod)%ddylds(0:mbod(bod)%neq),      &
        mbod(bod)%m_mvp(0:mbod(bod)%neq),       &
        mbod(bod)%m_mva(0:mbod(bod)%neq),       &
        mbod(bod)%loads(0:mbod(bod)%neq),       &
        mbod(bod)%fdamp(0:mbod(bod)%neq),       &
        mbod(bod)%gravlo(0:mbod(bod)%neq),      &
        mbod(bod)%Freact(0:mbod(bod)%neq),      &
        mbod(bod)%fnorm(0:mbod(bod)%neq),       &
        mbod(bod)%fcont(0:mbod(bod)%neq),       &
        mbod(bod)%f_fint(0:mbod(bod)%neq),      &
        mbod(bod)%f_earth(0:mbod(bod)%neq),     &
        mbod(bod)%vcm(0:mbod(bod)%neq),         &
        mbod(bod)%vel_change(0:mbod(bod)%neq),  &
        mbod(bod)%kinup_d2x1(0:mbod(bod)%neq),  &                                    
        mbod(bod)%kinup_d1x1(0:mbod(bod)%neq),  &
        mbod(bod)%kinup_Ground_d2x1(0:mbod(bod)%neq), &
        mbod(bod)%x1(0:mbod(bod)%neq),          &
        mbod(bod)%x1_orig(0:mbod(bod)%neq),     &
        mbod(bod)%x1_ini(0:mbod(bod)%neq),      &
        mbod(bod)%x1_change(0:mbod(bod)%neq),   &
        mbod(bod)%cdamp(0:mbod(bod)%neq),       &
        mbod(bod)%mvkv(0:mbod(bod)%neq))

      ! create global equation numbering g_g
      mbod(bod)%kdiag=zero
      DO iel=1,nels
        num=g_num(:,iel)
        CALL num_to_g(num,mbod(bod)%nf,g)
        CALL fkdiag(mbod(bod)%kdiag,g)
        mbod(bod)%g_g(:,iel)=g  
      END DO
     
      ! define global sparse matrix diagonals
      DO i=2,mbod(bod)%neq
        mbod(bod)%kdiag(i)=mbod(bod)%kdiag(i)+mbod(bod)%kdiag(i-1)
      END DO
    
      ! Allocate sparse matrices
      ALLOCATE(                                       &
        mbod(bod)%kv(mbod(bod)%kdiag(mbod(bod)%neq)), &
        mbod(bod)%mv(mbod(bod)%kdiag(mbod(bod)%neq)), &
        mbod(bod)%kp(mbod(bod)%kdiag(mbod(bod)%neq)))
    END IF
  END DO Body_fixities

  Freefield_fixities: DO bod=2,size(mbod)
    ! Find the boundary of the geometry from g_coord
    lowbound = minval(mbod(bod)%g_coord)
    upbound  = maxval(mbod(bod)%g_coord)  ! Computed earlier to consider only main domain coordinates
    
    ! g_coord_aux save the mesh coordinates in the initial setup
    mbod(bod)%g_coord_aux = mbod(bod)%g_coord

    ! Setup vertical boundary conditions
    mbod(bod)%nf = 1 ! first mark every equations as active
    j=1; k=1
    DO i=1,mbod(bod)%nn
      IF(mbod(bod)%g_coord(2,i)<lowbound+0.01) THEN
        mbod(bod)%nf(2,i)   = 0  
      END IF
      IF(bod>1.and.mbod(bod)%g_coord(1,i)>upbound-0.01) THEN 
        mbod(bod)%nf(1:2,i) = 0
      END IF
    END DO
    
    ! Setup Left Boundary Conditions
    mbod(bod)%tied_nn=0
    j=1
    leftbound=minval(mbod(bod)%g_coord(1,:))
    DO i=1,mbod(bod)%nn
      IF(mbod(bod)%g_coord(1,i)<leftbound+0.01_iwp)THEN
        mbod(bod)%tied_nn(j,1)=i
        j=j+1
      END IF  
    END DO

    ! Setup Right Boundary Conditions
    j=1
    rightbound=maxval(mbod(bod)%g_coord(1,:))
    DO i=1,mbod(bod)%nn
        IF(mbod(bod)%g_coord(1,i)>rightbound-0.01_iwp)THEN
            mbod(bod)%tied_nn(j,2)=i
            j=j+1
        END IF   
    END DO
    
    !- Fix pore presure in the not cornered nodes
    CALL formnf(mbod(bod)%nf)
    mbod(bod)%neq=MAXVAL(mbod(bod)%nf) ! number of equations including free-fields
    
    !-----------------------------------------------------------------------
    !-Repeat equations of one side of the domain to the other side to 
    ! simulate Tied-Degrees. Note that this operation is performed to bod 2 
    ! and 3 (i.e. bod>1), since bodies 2 and 3 are the Free-Fields 
    ! (i.e. 1D columns working individually with Tied-Degrees).
    !-----------------------------------------------------------------------
    ! This part makes the Free-Field bodies have tied degree of freedom
    ! between the left boundaries and the right boundaries
    IF(bod>1)THEN
      Tied_degree: DO i=1,mbod(bod)%ney+1
        mbod(bod)%nf(1:2,mbod(bod)%tied_nn(i,2)) = mbod(bod)%nf(1:2,mbod(bod)%tied_nn(i,1))
      END DO Tied_degree
    END IF
    
    ! Now the Number of Equation has been known, allocate some more variables
    ALLOCATE(                                                             &
      mbod(bod)%d1x1(0:mbod(bod)%neq),                                    &
      mbod(bod)%d2x1(0:mbod(bod)%neq),                                    &
      mbod(bod)%diag(0:mbod(bod)%neq),                                    &
      mbod(bod)%ddylds(0:mbod(bod)%neq),                                  &
      mbod(bod)%loads(0:mbod(bod)%neq),                                   &
      mbod(bod)%gravlo(0:mbod(bod)%neq),                                  &
      mbod(bod)%kdiag(mbod(bod)%neq),                                     &
      mbod(bod)%vcm(0:mbod(bod)%neq),                                     &
      mbod(bod)%x1(0:mbod(bod)%neq),                                      &
      mbod(bod)%kinup_d2x1(0:mbod(bod)%neq),                              &
      mbod(bod)%kinup_d1x1(0:mbod(bod)%neq),                              &
      mbod(bod)%f_fint(0:mbod(bod)%neq),                                  &
      mbod(bod)%P_ext(0:mbod(bod)%neq),                                   &
      mbod(bod)%cdamp(0:mbod(bod)%neq),                                   &
      mbod(bod)%ground_loads(0:mbod(bod)%neq),                            &
      mbod(bod)%kinup_Ground_d2x1(0:mbod(bod)%neq),                       &
      mbod(bod)%x1_acum(0:mbod(bod)%neq),                                 &
      mbod(bod)%c_force(0:mbod(bod)%neq),                                 &
      mbod(bod)%mf_force(0:mbod(bod)%neq),                                &
      mbod(bod)%residual(0:mbod(bod)%neq),                                &
      mbod(bod)%loads_ini(0:mbod(bod)%neq),                               &
      mbod(bod)%loads_end(0:mbod(bod)%neq),                               &
      mbod(bod)%loads_base(0:mbod(bod)%neq),                              &
      mbod(bod)%d1p1(0:mbod(bod)%neq),                                    &
      mbod(bod)%c_damp(0:mbod(bod)%neq)                                   &
    )
    
    ALLOCATE(                                                             &
      mbod(bod)%penpos(mbod(bod)%neq),                                    &
      mbod(bod)%penpos_v(mbod(bod)%neq),                                  &
      mbod(bod)%eq_pen(mbod(bod)%neq),                                    &
      mbod(bod)%eq_pen_v(mbod(bod)%neq)                                   &
    )
    
    mbod(bod)%g_g = 0
    nband   = 0
    bandval = 0
    
    ! This part create the steering vector. Dunno why not using num_to_g tho
    ! possibly because it can not deal with the tied degree. At the end, 
    ! calculate a the length of c_matrix and mf_matrix (AS)
    DO iel=1,mbod(bod)%nels
      mbod(bod)%g = 0
      num = mbod(bod)%g_num(:,iel)
      ! remember, g_g is organized per element and nf is organized per node
      mbod(bod)%g(1:ndof:2) = mbod(bod)%nf(1,num(:)) ! this one in x-dir
      mbod(bod)%g(2:ndof:2) = mbod(bod)%nf(2,num(:)) ! this one in y-dir
      mbod(bod)%g_g(:,iel)  = mbod(bod)%g
      ! bandval is calculated by taking the difference of g vector minmax
      ! only if it is a positive number (AS)
      bandval = MAXVAL(mbod(bod)%g,1,mbod(bod)%g>0) -                     &
                MINVAL(mbod(bod)%g,1,mbod(bod)%g>0)
      IF(nband<bandval) THEN
        nband = bandval
      END IF
    END DO
    
    ! c_matrix and mf_matrix are used to determine the kinematic update
    ! of each iterations in each time steps. F=ma & F=cv
    ALLOCATE(                                                             &
      mbod(bod)%c_matrix(mbod(bod)%neq,2*(nband+1)),                      &
      mbod(bod)%mf_matrix(mbod(bod)%neq,2*(nband+1))                      &
    )
    
    ! for each element, determine the diagonal index using fkdiag
    ! for each element, determine the diagonal index using fkdiag
    ! this will govern the shape of the sparse matrix
    mbod(bod)%kdiag=zero
    DO iel=1,mbod(bod)%nels
      mbod(bod)%g=mbod(bod)%g_g(:,iel) 
      CALL fkdiag(mbod(bod)%kdiag,mbod(bod)%g)
    END DO
    
    ! the step before only valid for each element, to build the whole matrix
    ! acummulate sum the diagonal terms
    DO i=2,mbod(bod)%neq
      mbod(bod)%kdiag(i)=mbod(bod)%kdiag(i)+mbod(bod)%kdiag(i-1)
    END DO

    ! Allocate sparse matrix
    ALLOCATE(                                                             &
      mbod(bod)%kv_CSR_aux(mbod(bod)%skylength),                          &
      mbod(bod)%mv(mbod(bod)%kdiag(mbod(bod)%neq)),                       &
      mbod(bod)%mvis(mbod(bod)%kdiag(mbod(bod)%neq))                      &
    )
  END DO Freefield_fixities

  !===========================================================================AS
  ! Variable initialization
  !===========================================================================AS
  READ(10,*)dtim,k0,aiters,tol,limit,nsrf

  ! Material parameter valid also valid for Free Field elements [AS]
  DO bod=1,size(mbod)
    mbod(bod)%Young=mbod(bod)%prop(2,1)
    mbod(bod)%Poiss=mbod(bod)%prop(3,1)
    mbod(bod)%mpcp=mbod(bod)%prop(4,1)
    cpeak=mbod(bod)%prop(4,1)
    mpcr=mbod(bod)%prop(5,1) 
    phi=mbod(bod)%prop(6,1) 
  END DO
    
  Maxplastic=0.3
  SModulus=(mpcr-cpeak)/Maxplastic
  cont2=0
  evpt=zero
  elso=zero
  fmax=-10.0
  fnew=-10
  stol=0.01_iwp
  Maxfric=1.0_iwp
  Minfric=0.50
  Fricmod=(Minfric-Maxfric)/Maxplastic

  
  !===========================================================================AS
  ! Initial Conditions
  !===========================================================================AS
  
  initial_conditions: DO bod=1,size(mbod)
    IF(bod==1)THEN
      mbod(bod)%m_acc         = zero
      mbod(bod)%m_stress      = zero
      mbod(bod)%m_stress_ini  = zero
      mbod(bod)%m_stress_change = zero
      mbod(bod)%m_stress_prev = zero
      mbod(bod)%accp          = zero
      mbod(bod)%ins           = zero
      mbod(bod)%a_ins         = zero
      mbod(bod)%x1_orig       = zero
      mbod(bod)%d1x1          = zero
      mbod(bod)%d2x1          = zero
      mbod(bod)%accb          = zero
      mbod(bod)%vccp          = zero
      mbod(bod)%kinup_d1x1    = zero
      mbod(bod)%kinup_d2x1    = zero
      mbod(bod)%epsinvacum    = zero
      mbod(bod)%Devstress     = zero
      mbod(bod)%mean_stress   = zero
      mbod(bod)%diag          = zero
      mbod(bod)%loads         = zero
      mbod(bod)%gravlo        = zero
      mbod(bod)%ddylds        = zero
      mbod(bod)%m_velocity    = zero
      mbod(bod)%mpyield       = zero
      mbod(bod)%nodecont_f    = zero
      mbod(bod)%kv            = zero
      mbod(bod)%mv            = zero
      mbod(bod)%flag          = 0
    ELSE
      mbod(bod)%m_acc         = zero
      mbod(bod)%m_stress      = zero
      mbod(bod)%m_stress_ini  = zero
      mbod(bod)%m_stress_efe  = zero
      mbod(bod)%accp          = zero
      mbod(bod)%a_ins         = zero
      mbod(bod)%d1x1          = zero
      mbod(bod)%d2x1          = zero
      mbod(bod)%accb          = zero
      mbod(bod)%vccp          = zero
      mbod(bod)%flag          = 0
      mbod(bod)%m_pore        = zero
      mbod(bod)%x1            = zero
      mbod(bod)%x1_acum       = zero
      mbod(bod)%eps_acum      = zero
      mbod(bod)%Devstress     = zero
      mbod(bod)%ins           = zero
      mbod(bod)%diag          = zero
      mbod(bod)%m_velocity    = zero
      mbod(bod)%mean_stress   = zero
      mbod(bod)%mvis          = zero
      mbod(bod)%c_fact        = zero
      mbod(bod)%kinup_d1x1    = zero
      mbod(bod)%kinup_d2x1    = zero
      mbod(bod)%loads_end     = zero
      mbod(bod)%loads_base    = zero
      mbod(bod)%P_ext         = zero
      mbod(bod)%Cc            = zero
      mbod(bod)%statev        = zero
      mbod(bod)%statev(:,7)   = props(16)    
      mbod(bod)%statev(:,13)  = 1.0_iwp
      mbod(bod)%statev(:,52)  = zero
      mbod(bod)%c_matrix      = zero
      mbod(bod)%mf_matrix     = zero
      mbod(bod)%d1p1          = zero
      mbod(bod)%c_damp        = zero
      mbod(bod)%Young         = mbod(1)%prop(2,1)
      mbod(bod)%Poiss         = mbod(1)%prop(3,1)
    END IF
  END DO initial_conditions

  ! Create consitutive matrix for each body
  Constitutive: DO bod=1,size(mbod)
    CALL deemat(mbod(bod)%dee,mbod(bod)%Young,mbod(bod)%Poiss)
  END DO Constitutive

  ! Define gravity field
  Gravf=10.00_iwp
  k0=0.50_iwp
  g_matrix=(/0.0,-0.00/)  !--Gravity acting in the vertical direction
  DO bod=1,size(mbod)
    mbod(bod)%g_matrix=(/0.0,-10.0/)
  END DO
  
  ! Print to paraview
  DO bod=1,size(mbod)
    IF(bod==1)THEN
      CALL point_viz2(0,bod,argv,mbod(bod)%gm_coord,mbod(bod)%m_stress,mbod(bod)%m_stress_change,   &
        mbod(bod)%epsinvacum,mbod(bod)%a_ins,mbod(bod)%Devstress,mbod(bod)%mean_stress,             &
        mbod(bod)%mpyield,mbod(bod)%mpcp,mbod(bod)%m_velocity,mbod(bod)%m_acc,mbod(bod)%nmps,nlen)
      CALL paraview2(cont2,bod,argv,g_coord,g_num,mbod(bod)%nf,nels,nod,nn,nlen,                    &
        mbod(bod)%diag,mbod(bod)%ddylds,mbod(bod)%d1x1,mbod(bod)%d2x1,                              &
        mbod(bod)%gravlo,mbod(bod)%loads,mbod(bod)%normal,mbod(bod)%fcont,mbod(bod)%kv,mbod(bod)%mv,&
        mbod(bod)%kdiag,mbod(bod)%vcm,mbod(bod)%f_fint)     
    ElSE
        CALL paraview(                                              &
            cont2,                                                  &
            bod,                                                    &
            argv,                                                   &
            (mbod(bod)%ny1+1),                                      &
            mbod(bod)%ny1,                                          &
            mbod(bod)%g_coord_aux,                                  &
            mbod(bod)%g_num,                                        &
            mbod(bod)%nf,                                           &
            mbod(bod)%nels,                                         &
            nod,                                                    &
            mbod(bod)%nn,                                           &
            nlen,                                                   &
            mbod(bod)%diag,                                         &
            (-1.0)*mbod(bod)%ddylds,                                &
            mbod(bod)%kinup_d1x1,                                   &
            mbod(bod)%kinup_d2x1,                                   &
            mbod(bod)%gravlo,                                       &
            mbod(bod)%x1_acum,                                      &
            mbod(bod)%P_ext,                                        &
            mbod(bod)%mv,                                           &
            mbod(bod)%kdiag,                                        &
            (-1.0)*mbod(bod)%mf_force,                              &
            mbod(bod)%loads,                                        &
            -mbod(bod)%c_damp,                                      &
            mbod(bod)%kv_CSR_aux,                                   &
            mbod(bod)%penpos,1                                      &
        )     
        CALL point_viz(                                             &
            cont2,                                                  &
            bod,                                                    &
            argv,                                                   &
            mbod(bod)%gm_coord,                                     &
            mbod(bod)%m_stress_efe,                                 &
            mbod(bod)%eps_acum,                                     &
            mbod(bod)%statev(:,7),                                  &
            mbod(bod)%a_ins,mbod(bod)%Devstress,                    &
            mbod(bod)%mean_stress,                                  &
            mbod(bod)%m_pore(1,:),                                  &
            mbod(bod)%statev(:,50),                                 &
            mbod(bod)%m_velocity,                                   &
            mbod(bod)%m_acc,mbod(bod)%nmps,                         &
            nlen,                                                   &
            1                                                       &
        )
    END IF
  END DO

  ! Obtain ground movement data
  READ(300,*)accdata
  DO bod=1,size(mbod)
    ALLOCATE(mbod(bod)%ground_acc(accdata))
  END DO
  
  READ(300,*)mbod(1)%ground_acc
  DO bod=2,size(mbod)
    mbod(bod)%ground_acc=mbod(1)%ground_acc
  END DO
  
  
!===============================================================================
! ----------------------- BEGINNING OF THE TIME STEP ---------------------------
!===============================================================================

  !============================================================================AS
  ! Iteration Setups
  !============================================================================AS
  
  stable=.true.
  step=0
  time_steps: DO w=1,100000
  step=step+1 
  
  ! Reset variables
  Reset_Variables: DO bod=1,size(mbod)
    ! general variables
    mbod(bod)%m_mvp      = zero
    mbod(bod)%m_mva      = zero
    mbod(bod)%valuesg    = zero
    mbod(bod)%GIMP_nodes = zero
    mbod(bod)%diag       = zero
    mbod(bod)%vcm        = zero
    mbod(bod)%gravlo     = zero
    mbod(bod)%cdamp      = zero
    mbod(bod)%f_earth    = zero
    mbod(bod)%loads      = zero
    mbod(bod)%x1         = zero
    mbod(bod)%d1x1       = zero
    mbod(bod)%d2x1       = zero
    mbod(bod)%kinup_d1x1 = zero
    mbod(bod)%kinup_d2x1 = zero
    mbod(bod)%kinup_Ground_d2x1 = zero
    
    ! additional variables
    mbod(bod)%normal     = zero
    mbod(bod)%temp_d1x1  = zero
    mbod(bod)%temp_d2x1  = zero
    mbod(bod)%f_fint     = zero
    mbod(bod)%fcont      = zero
    mbod(bod)%fnorm      = zero
    mbod(bod)%eps_m      = zero
    mbod(bod)%eps_1      = zero
    mbod(bod)%eps_2      = zero
    mbod(bod)%v_field    = zero
    mbod(bod)%m_field    = zero
    mbod(bod)%a_field    = zero
    mbod(bod)%vel_change = zero
    mbod(bod)%x1_ini     = zero
    mbod(bod)%x1_change  = zero
    mbod(bod)%nodecont   = zero
    mbod(bod)%nodecont_f = zero
      
    ! freefield variables
    mbod(bod)%diag       = zero
    mbod(bod)%c_force    = zero
    mbod(bod)%mf_force   = zero
    mbod(bod)%loads_end  = zero
    mbod(bod)%loads_base = zero
    mbod(bod)%ground_loads = zero
    mbod(bod)%penpos     = 0
    mbod(bod)%penpos_v   = 0
    mbod(bod)%eq_pen     = 0
    mbod(bod)%eq_pen_v   = 0
  END DO Reset_Variables
  
  ! Calculate diagonal mass matrix (diag)
  Body_Solution: DO bod=1,size(mbod)
    mvval=1
    IF(mbod(bod)%nmps>1)THEN 
      !---Compute mesh momentum and mass
      IF(bod==1)THEN
        Fext_Mass_P: DO i=1,mbod(bod)%nmps
          IF(mbod(bod)%nmps>1)THEN 
            iel=mbod(bod)%a_ele(i)    
            IF(smethod==2.or.smethod==3)THEN 
              CALL GIMP_nodsup(i,g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp, &
                mbod(bod)%mpoints,mbod(bod)%valuesg,mbod(bod)%gimptol,neighb,nny,mbod(bod)%a_ele,mbod(bod)%GIMP_nodes)
              values=mbod(bod)%valuesg(i)
      
              ALLOCATE(derextend(nodof,values),jac_coord(values,nodof),funextend2(values*2,2),eldddylds(values*2),beeextend(nst,values*nodof))

              eldddylds=zero;funextend2=zero;beeextend=zero;derextend=zero;jac_coord=zero
              !--creates the matrix with the values of shape functions and derivatives
              CALL GIMP_funder2(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)  
              !--collect in eldddylds the nuber of the equations in the support domain
              CALL gimpfunform(i,iel,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g,mvval)  
                                
              mbod(bod)%diag(eldddylds) = mbod(bod)%diag(eldddylds) + MATMUL(funextend2,mbod(bod)%m_mass(i)*delta_1)
              
              mbod(bod)%m_mvp(0)=zero;mbod(bod)%diag(0)=zero;mbod(bod)%gravlo(0)=zero;mbod(bod)%m_mva(0)=zero
              DEALLOCATE(derextend,funextend2,eldddylds,jac_coord,beeextend)
            END IF
          END IF
        END DO Fext_Mass_P
      ELSE !- Free Fields bodies
        i=1; j=1; n=1
        Fext_domain_ff: DO k=1,mbod(bod)%nmps
          mbod(bod)%diag=zero
          iel=mbod(bod)%a_ele(k) 
          
          ! obtain the element steering vector.
          ! then set the degree of freedom in x-direction to 0
          ! because the body mass only go downwards
          mbod(bod)%g=mbod(bod)%g_g(:,iel) 
          mbod(bod)%g(1:mbod(bod)%ntot:2)=0

          num=mbod(bod)%g_num(:,iel)
          coord=TRANSPOSE(mbod(bod)%g_coord(:,num)) 
          CALL shape_der(der,mbod(bod)%mpoints,k)
          jac=MATMUL(der,coord)
          jac(2,1)=zero; jac(1,2)=zero ! ensure the off diagonal is zero
          det=determinant(jac)
      
          CALL shape_fun(fun,mbod(bod)%mpoints,k)
          CALL ecmat2(ecm,fun,ndof,nodof)

          ! mass lumping per element
          DO m=1,SIZE(fun)*2
            ecm_acum=zero
            DO j=1,SIZE(fun)*2
              ecm_acum(m,m)=ecm_acum(m,m)+(ecm(m,j))
            END DO
            ecm(m,:)=zero
            ecm(m,m)=ecm_acum(m,m)
          END DO
          mm_s=mm_s+ecm*mbod(bod)%m_dens(k)*det*weights(i)
      
          ! accummulate all the diagonal mass matrix from all the elements
          ! into mbod(bod)%diag
          CALL formlump(mbod(bod)%diag,mm_s,mbod(bod)%g(1:18))
          mbod(bod)%diag(0)=zero

          ! calculate the actual gravity loading (gravlo = mm*grav_acceleration)
          mbod(bod)%gravlo(mbod(bod)%g(1:ndof)) = mbod(bod)%gravlo(mbod(bod)%g(1:ndof)) + &
                                                  mbod(bod)%diag(mbod(bod)%g(1:ndof))*(-Gravf)  
      
          mbod(bod)%gravlo(0)=zero

          ! reset mm_s for next integration point
          ! the index is to cycle through the weigths
          ! Note: this will break if a_ele is not ordered correctly
          mm_s=zero
          i=i+1
          IF(i>nip)i=1

        END DO Fext_domain_ff
      END IF
    END IF
  END DO Body_Solution
   
  ! Compute domain stiffnes matrices (kv), domain mass matrices (mv), and their modified stiffness (kp)
  ! This matrices are computed using only the double mapping technique and GIMP, regular MPM is not used
  Stiffness:DO bod=1,size(mbod)
    IF(bod==1)THEN
      IF(mbod(bod)%nmps>1)THEN 
        mbod(bod)%kv=zero
        mbod(bod)%mv=zero
        mbod(bod)%kp=zero
        KM_MV:DO k=1,mbod(bod)%nmps
          ! Build Stiffness Matrix
          IF(mbod(bod)%kconst(1)==1)THEN
            CALL sample(element,points,weights)
            !-constant stiffness for each activated element
            K_const:DO iel=1,nels 
              !IF(mbod(bod)%d_ele(iel)==1)THEN
              num=g_num(:,iel)
              !- The next conditin is to find if the element is under the influence of any material point in the domain. 
              !- If is influenced, then the FEM stiffnes is computed for that element
              IF(((mbod(bod)%nf(1,num(1))>0.or.mbod(bod)%nf(2,num(1))>0).and.     &
                  (mbod(bod)%nf(1,num(2))>0.or.mbod(bod)%nf(2,num(2))>0).and.     &
                  (mbod(bod)%nf(1,num(3))>0.or.mbod(bod)%nf(2,num(3))>0).and.     &
                  (mbod(bod)%nf(1,num(4))>0.or.mbod(bod)%nf(2,num(4))>0)).or.     &
                mbod(bod)%c_ele(iel)>0)THEN

                coord=TRANSPOSE(g_coord(:,num))
                g=mbod(bod)%g_g(:,iel) 
                km=zero
                gauss_points:DO i=1,nip
                  CALL shape_der(der,points,i)
                  jac=MATMUL(der,coord)
                  det=determinant(jac)
                  CALL invert(jac)
                  deriv=MATMUL(jac,der)
                  CALL beemat(bee,deriv)
                  km=km+MATMUL(MATMUL(TRANSPOSE(bee),mbod(bod)%dee),bee)*det*weights(i)
                END DO gauss_points    
                CALL fsparv(mbod(bod)%kv,km,g,mbod(bod)%kdiag)     
              END IF
            END DO K_const   
          ELSE    
            El_stiff:DO i=1,4 !This loop goes until 4 since 4 is the maximum number of elements affected by a material point
              iel=mbod(bod)%elemmpoints(k,i)
              Full_el:IF(mbod(bod)%elemmpoints(k,i)>0)THEN !--If is true, an element is affected by a material point (material point can be outside the element)
                num=g_num(:,iel)
                nod_num(:,1)=num
                coord=TRANSPOSE(g_coord(:,num))
                g=mbod(bod)%g_g(:,iel)
                km_gauss=zero   
                waverage=4.0/nip

                Double_map:IF(smethod>1)THEN !-Double mapping technique
                  CALL iGIMP_funder3(k,mbod(bod)%mpoints,mbod(bod)%lp_mp,nip,coord, &
                    cellsize,mbod(bod)%gm_coord,mbod(bod)%gimptol,nod_num,g_coord,  &
                    mbod(bod)%a_ele,mbod(bod)%c_ele,iel,der,fun)

                  IF(fun(1)<=0.or.fun(2)<=0.or.fun(3)<=0.or.fun(4)<=0)THEN
                    PRINT*,'Stiffness problem'
                    PRINT*,'shape functions',fun(:)
                    PRINT*,''
                    PRINT*,'material point',k
                      PRINT*,''
                    PRINT*,'element',iel
                    PRINT*,''
                    PRINT*,'lp',iel
                    PAUSE
                  END IF    
                  CALL sample(element,points,weights)
                  DO s=1,nip
                    CALL shape_der(der_gauss,points,s)
                    CALL shape_fun(fun_gauss,points,s)
                    scalar_dee=fun*fun_gauss*waverage
                    sum_scalar=SUM(scalar_dee)
                    dee_scal=mbod(bod)%dee*sum_scalar
                    jac=MATMUL(der_gauss,coord)
                    det=determinant(jac)
                    CALL invert(jac)
                    deriv_gauss=MATMUL(jac,der_gauss)
                    CALL beemat(bee_gauss,deriv_gauss)
                    IF(mbod(bod)%c_ele(iel)<1) waverage=1.5_iwp
                    
                    km_gauss = km_gauss + MATMUL(MATMUL(TRANSPOSE(bee_gauss),dee_scal),bee_gauss)*det*weights(i)
                  
                  END DO
                  
                  CALL fsparv(mbod(bod)%kv,km_gauss,g,mbod(bod)%kdiag) ! STIFFNESS MATRIX (KV)
                
                END IF Double_map

              END IF Full_el
            END DO El_stiff
          END IF  
          
          ! Build diagonal mass matrix  
          values=mbod(bod)%valuesg(k)
          ALLOCATE(derextend(nodof,values),funextend2(values*2,2),eldddylds(values*2),mm_gimp(values*2,values*2))
          CALL GIMP_funder2(k,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)
          mm_gimp=zero
          mvval=2
          CALL gimpfunform(k,iel,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g,mvval) 
          
          mvval=1
          a=1
          DO j=1,values
            DO q=1,2
              mm_gimp(a,a)=funextend2(a,q)*mbod(bod)%m_mass(k) 
              a=a+1
            END DO
          END DO
          
          CALL fsparv(mbod(bod)%mv,mm_gimp,eldddylds,mbod(bod)%kdiag) ! diagonal mass matrix (MV)
          
          DEALLOCATE(derextend,funextend2,eldddylds,mm_gimp)
        END DO KM_MV
      END IF
      
    ELSE ! Freefields
      
      i=1; j=1; n=1; q=1
      mbod(bod)%ia_aux     = 0
      mbod(bod)%ja_aux     = 0
      mbod(bod)%kv_CSR_aux = zero
      nnzacum              = 0
      mbod(bod)%skylength  = 0
      KM_MV_2D:DO k=1,mbod(bod)%nmps
        IF(i==1)THEN
          mbod(bod)%km      = zero 
          mm_s              = zero
          mbod(bod)%MOD_MTX = zero
          mm_acum           = zero
          mm_mf_acum        = zero
        END IF
        
        ! this is similar with the body solution
        iel=mbod(bod)%a_ele(k)
        num=mbod(bod)%g_num(:,iel)
        mbod(bod)%g=0
        mbod(bod)%g=mbod(bod)%g_g(:,iel) ! notice now both-dir exist
        coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
        CALL shape_der(der,mbod(bod)%mpoints,k)
        jac=MATMUL(der,coord)
        jac(2,1)=zero; jac(1,2)=zero ! ensure the off diagonal is zero
        det=determinant(jac)
        CALL invert(jac)
        deriv=MATMUL(jac,der)
        CALL beemat(bee,deriv)
            
        ! Build the stiffness matrix (km)
        CALL deemat(mbod(bod)%dee,mbod(bod)%Young,mbod(bod)%Poiss)
        mbod(bod)%km = mbod(bod)%km + MATMUL(MATMUL(TRANSPOSE(bee),mbod(bod)%dee),bee)*det*weights(i)

        ! Build the diagonal mass matrix (mv)
        CALL shape_fun(fun,mbod(bod)%mpoints,k)
        CALL ecmat2(ecm,fun,ndof,nodof)
        DO m=1,SIZE(fun)*2
          ecm_acum=zero
          DO j=1,SIZE(fun)*2
            ecm_acum(m,m)=ecm_acum(m,m)+(ecm(m,j))
          END DO
          ecm(m,:)=zero
          ecm(m,m)=ecm_acum(m,m)
        END DO
        mm_s=mm_s+ecm*mbod(bod)%m_dens(k)*det*weights(i)

        ! at the end of each element form the modified stiffness matrix (kp or MOD_MTK)
        ! according to the incremental FEM iterations (MOD_MTK = 4MMS/dtim**2 + 2C/dtim + KGC) ; C = fk*KGC + fm*MMS 
        IF(i>nip-1)THEN
          mbod(bod)%KGC = zero
          mbod(bod)%MMS = zero
          mbod(bod)%KGC = mbod(bod)%km
          mbod(bod)%MMS = mm_s
          mbod(bod)%MOD_MTX = mbod(bod)%KGC +                        &
                              4.0_iwp*mbod(bod)%MMS/dtim**2.0_iwp +  &
                              2.0_iwp*(mbod(bod)%KGC*fk +            &
                              mbod(bod)%MMS*fm)/dtim
              
          CALL formspars_unsym(   &
            mbod(bod)%ntot,       &
            mbod(bod)%g,          &
            mbod(bod)%MOD_MTX,    &
            mbod(bod)%ia_aux,     &
            mbod(bod)%ja_aux,     &
            mbod(bod)%kv_CSR_aux, &
            mbod(bod)%skylength   &
          )
            
          mm_acum    = zero
          mm_mf_acum = zero
          mm_s       = zero
        END IF    
        i=i+1
        IF(i>nip)i=1
      END DO KM_MV_2D
    END IF
  END DO Stiffness
 
 
  ! For MPM, the modified stiffness is computed here (kp)
  Combined_Stiffness:DO bod=1,1  
    !fm=0.052359878_iwp;fk=0.000265258_iwp
    fm=0.0_iwp;fk=0.0_iwp
    ! Calculate modified stiffness matrix
    IF(mbod(bod)%nmps>1)THEN 
      !--The stiffness matrix and the mass matrix stays constant through the analisys 
      mbod(bod)%kp=4.0_iwp*mbod(bod)%mv/dtim**2.0_iwp + 2.0_iwp*fm*mbod(bod)%mv/dtim + 2.0_iwp*fk*mbod(bod)%kv/dtim + mbod(bod)%kv
      CALL sparin(mbod(bod)%kp,mbod(bod)%kdiag)
    END IF
  
    !*** AS's THESIS
    ! Apply penalty position for boundary element nodes [TODO]
    
    
  END DO Combined_Stiffness


  !---------------------------------------------------------------------------------------AS
  ! These loops are to preparing the solvers for solving the Freefield bodies with Pardiso
  !---------------------------------------------------------------------------------------AS
  
  ! Change the sparse data-structure from 3-array forms to CSR form
  DO bod=2,size(mbod)
    nnzacum=0
    CALL sortadd(                                                 &
      mbod(bod)%skylength,                                        &
      mbod(bod)%ia_aux,                                           &
      mbod(bod)%ja_aux,                                           &
      mbod(bod)%kv_CSR_aux,                                       &
      mbod(bod)%neq+1,                                            &
      nnzacum,                                                    &
      mbod(bod)%penpos                                            &
    )
    ALLOCATE(                                                     &
      mbod(bod)%ia(mbod(bod)%neq+1),                              &
      mbod(bod)%ja(nnzacum),                                      &
      mbod(bod)%kv_CSR(nnzacum)                                   &
    )
    mbod(bod)%kv_CSR = zero
    mbod(bod)%ja     = 0
    mbod(bod)%ia     = 0
    mbod(bod)%kv_CSR(1:nnzacum)     = mbod(bod)%kv_CSR_aux
    mbod(bod)%ja(1:nnzacum)         = mbod(bod)%ja_aux
    mbod(bod)%ia(1:mbod(bod)%neq+1) = mbod(bod)%ia_aux(1:mbod(bod)%neq+1)
  END DO

  ! Compute conventional mass (mf_matrix) and damping matrix (c_matrix) for calculating forces
  DO bod=2,size(mbod)
    mbod(bod)%mv        = zero
    mbod(bod)%c_matrix  = zero
    mbod(bod)%mvis      = zero
    mbod(bod)%mf_matrix = zero
    mbod(bod)%km        = zero
    mbod(bod)%CCS       = zero
    ecm                 = zero
    mm_s                = zero

    i=1; n=1; j=1; q=1
    Mass_matrix:DO k=1,mbod(bod)%nmps
      ! again with the routines of making the bee matrix
      iel = mbod(bod)%a_ele(k)
      num = mbod(bod)%g_num(:,iel)
      coord = TRANSPOSE(mbod(bod)%g_coord(:,num))
      mbod(bod)%g = 0
      mbod(bod)%g = mbod(bod)%g_g(:,iel)  
      CALL shape_der(der,mbod(bod)%mpoints,k)
      jac = MATMUL(der,coord)
      det = determinant(jac)
      CALL invert(jac) 
      CALL shape_fun(fun,mbod(bod)%mpoints,k)
      deriv = MATMUL(jac,der)
      CALL beemat(bee,deriv)  
      CALL deemat(mbod(bod)%dee,mbod(bod)%Young,mbod(bod)%Poiss)

      mbod(bod)%km = mbod(bod)%km + MATMUL(MATMUL(TRANSPOSE(bee),mbod(bod)%dee),bee)*det*weights(i)
      
      CALL ecmat2(ecm,fun,ndof,nodof)
      DO m=1,SIZE(fun)*2
        ecm_acum=zero
        DO j=1,SIZE(fun)*2
            ecm_acum(m,m)=ecm_acum(m,m)+(ecm(m,j))
        END DO
        ecm(m,:)=zero
        ecm(m,m)=ecm_acum(m,m)
      END DO
      mm_s=mm_s+ecm*mbod(bod)%m_dens(k)*det*weights(i)
      
      ! but now instead of making the modified stiffness
      ! lets build the c_matrix and mf_matrix instead
      ! mf_matrix is mass matrix and c matrix is the Rayleigh damping
      ! these matrix aren't used for solving the displacement
      ! but rather to get the kinetic updates / stress recovery
      IF(i>nip-1)THEN
        mbod(bod)%MMS = zero
        mbod(bod)%CCS = zero
        
        mbod(bod)%MMS = mm_s
        mbod(bod)%CCS = mbod(bod)%km*fk+mbod(bod)%MMS*fm
            
        CALL formtb(mbod(bod)%c_matrix,mbod(bod)%CCS,mbod(bod)%g)
        CALL formtb(mbod(bod)%mf_matrix,mbod(bod)%MMS,mbod(bod)%g)

        mm_s=zero
        mbod(bod)%km=zero
      END IF
      i=i+1
      IF(i>nip)i=1
    END DO Mass_matrix
  END DO

  ! LU Decomposition
  DO bod=2,size(mbod) 
    ALLOCATE(mbod(bod)%iparm(64))
    mbod(bod)%iparm=0
    mbod(bod)%pt = 0    ! pointer initialization
    maxfct = 1
    mnum = 1
    nrhs = 1
    mtype = 11          ! Unsymmetric indefinite
    error = 0  
    msglvl = 0
    mbod(bod)%iparm(33) = 1
    mbod(bod)%iparm(3) = 4
    
    phase = 11  !// check matrix consistency
    CALL pardiso(                                                 &
      mbod(bod)%pt,                                               &
      maxfct,                                                     &
      mnum,                                                       &
      mtype,                                                      &
      phase,                                                      &
      mbod(bod)%neq,                                              &
      mbod(bod)%kv_CSR,                                           & 
      mbod(bod)%ia,                                               &
      mbod(bod)%ja,                                               &
      idum,                                                       &
      nrhs,                                                       &
      mbod(bod)%iparm,                                            &
      msglvl,                                                     &
      ddum,                                                       &
      ddum,                                                       &
      error,                                                      &
      dparm                                                       &
    )

    ! pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    ! A*X = B

    phase = 22  !// LU decompose
    CALL pardiso(                                                   &
      mbod(bod)%pt,                                               &
      maxfct,                                                     &
      mnum,                                                       &
      mtype,                                                      &
      phase,                                                      &
      mbod(bod)%neq,                                              &
      mbod(bod)%kv_CSR,                                           &
      mbod(bod)%ia,                                               &
      mbod(bod)%ja,                                               &
      idum,                                                       &
      nrhs,                                                       &
      mbod(bod)%iparm,                                            &
      msglvl,                                                     &
      ddum,                                                       &
      ddum,                                                       &
      error,                                                      &
      dparm                                                       &
    )
    
  END DO


  !============================================================================AS
  ! Iteration loops 
  !============================================================================AS
  
  !--- Solve freefield bodies stresses and displacement levels
  iters=0
  newrapconv1=.false.
  FF_Newton_Rhapson: DO WHILE(iters<limit+1)        
    iters=iters+1

    !---------------------------------------------------------------------------AS
    ! Determine Nodal Kinematics
    !---------------------------------------------------------------------------AS

    ! Calculate actual nodal acceleration 
    DO bod=2,size(mbod)
      IF(tstep<=accdata.and.tstep>=1)THEN
        DO i=1,(mbod(bod)%ney+1)*(mbod(bod)%nex)
          IF(mbod(bod)%g_coord(2,i)<lowbound+0.01_iwp)THEN
            valone = ground_acc(tstep)
            mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i)) =                    &
                mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i))+                 &
                (2.0_iwp*ground_acc(tstep)-mbod(bod)%kinup_d2x1(mbod(bod)%nf(1,i)))
            mbod(bod)%kinup_Ground_d2x1(0)=zero
          END IF    
        END DO
      END IF
    END DO

    !-Kinematic update of acceleration and velocity without 
    ! considering any boundary condition (contact)

    DO bod=2,size(mbod)
      mbod(bod)%kinup_d2x1 = (4.0_iwp*mbod(bod)%x1/dtim**2.0_iwp) - &
                              (4.0_iwp*mbod(bod)%d1x1/dtim) -        &
                              mbod(bod)%d2x1         
      mbod(bod)%kinup_d1x1 = 2.0_iwp*mbod(bod)%x1/dtim -            &
                              mbod(bod)%d1x1
      mbod(bod)%kinup_d2x1(0) = zero
      mbod(bod)%kinup_d1x1(0) = zero
      mbod(bod)%d1p1(0)       = zero
    END DO 

    !---------------------------------------------------------------------------AS
    ! Determine Nodal Forces
    !---------------------------------------------------------------------------AS
    
    DO bod=2,size(mbod)
      mbod(bod)%ddylds=zero
      n=1; j=1; a=1
      Fint_load:DO i=1,mbod(bod)%nmps
        iel=mbod(bod)%a_ele(i)   
        num=mbod(bod)%g_num(:,iel)
        coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
        mbod(bod)%g=mbod(bod)%g_g(:,iel)
        CALL shape_der(der,mbod(bod)%mpoints,i)
        jac=MATMUL(der,coord)
        det=determinant(jac)
        CALL invert(jac)
        deriv=MATMUL(jac,der)
        CALL beemat(bee,deriv)
        sigma=mbod(bod)%m_stress(:,i) 
        mbod(bod)%ddylds(mbod(bod)%g(1:ndof)) =                  &
            mbod(bod)%ddylds(mbod(bod)%g(1:ndof)) +              &
            MATMUL(TRANSPOSE(bee),sigma)*det*weights(a)
        a=a+1
        IF(a>nip)a=1
        mbod(bod)%ddylds(0)=zero 
      END DO Fint_load
    END DO 
    
    Force_MA: DO bod=2,size(mbod)
      ! This is force from the damping F=Cv
      CALL bantmul(                                               &
          mbod(bod)%c_matrix,                                     &
          mbod(bod)%kinup_d1x1,                                   &
          mbod(bod)%c_force                                       &
      )
      ! This is force from the general body motion F=Ma
      CALL bantmul(                                               &
        mbod(bod)%mf_matrix,                                      &
        mbod(bod)%kinup_d2x1,                                     &
        mbod(bod)%mf_force                                        &
      )
      mbod(bod)%vcm(0)      = zero
      mbod(bod)%c_force(0)  = zero
      mbod(bod)%mf_force(0) = zero
    END DO Force_MA 

    Ground_F: DO bod=2,size(mbod)
        mbod(bod)%kinup_Ground_d2x1(0)=zero
        ! This is force from the ground motion also F=Ma
        CALL bantmul(                                               &
            mbod(bod)%mf_matrix,                                    &
            mbod(bod)%kinup_Ground_d2x1,                            &
            mbod(bod)%ground_loads                                  &
        )
        mbod(bod)%ground_loads(0)=zero
    END DO Ground_F

    !---------------------------------------------------------------------------AS
    ! Solve the equations
    !---------------------------------------------------------------------------AS
    
    FEM_Displacements: DO bod=2,size(mbod)
      ! Combine all the forces: gravity loads, internal loads, ground loads, 
      ! displacement loads, and damping load
      mbod(bod)%loads    = zero
      mbod(bod)%residual = zero
      mbod(bod)%loads    =  mbod(bod)%gravlo - mbod(bod)%ddylds +  &
                            mbod(bod)%ground_loads -               &
                            mbod(bod)%mf_force - mbod(bod)%c_force
      mbod(bod)%loads(0)=zero

      ! Check whether the calculated loads has converge (saved in newrapconv1)
      IF(iters==1) mbod(bod)%loads_ini = mbod(bod)%loads

      ! Solve the system of equations (Ma+Cv+Ku = F)
      ! pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
      ! A*X = B ; X = mbod(bod)%loads(1:mbod(bod)%neq) ; B = mbod(bod)%residual(1:mbod(bod)%neq)
      mbod(bod)%iparm(8)  = 10 ! max numbers of iterative refinement steps
      phase = 33  ! solve equations
      call pardiso(                          &
        mbod(bod)%pt,                        &
        maxfct,                              &
        mnum,                                &
        mtype,                               &
        phase,                               &
        mbod(bod)%neq,                       &
        mbod(bod)%kv_CSR,                    &
        mbod(bod)%ia,                        &
        mbod(bod)%ja,                        &
        idum,                                &
        nrhs,                                &
        mbod(bod)%iparm,                     &
        msglvl,                              &
        mbod(bod)%loads(1:mbod(bod)%neq),    &
        mbod(bod)%residual(1:mbod(bod)%neq), &
        error                                &
      )
  
      mbod(bod)%loads(0)    = zero
      mbod(bod)%residual(0) = zero 
      mbod(bod)%loads       = mbod(bod)%residual  ! assign displacement to loads (weird indeed)
      mbod(bod)%loads_base  = mbod(bod)%loads
      mbod(bod)%x1          = mbod(bod)%x1+mbod(bod)%loads ! add displacement increments
      mbod(bod)%x1(0)       = zero
  
      ! Remember that loads here is already become a displacement
      ! if the last disp. has converge, the new disp. will also be already
      ! converged. So save it in the loads_end variable
      IF(newrapconv1)mbod(bod)%loads_end=mbod(bod)%loads
    END DO FEM_Displacements

    !---------------------------------------------------------------------------AS
    ! Update Nodal Kinematics
    !---------------------------------------------------------------------------AS
    
    ! calculate the acceleration and velocity increments
    DO bod=2,size(mbod)
        ! acceleration
        mbod(bod)%kinup_d2x1    = (4.0_iwp*mbod(bod)%x1/dtim**2.0_iwp) - &
                                  (4.0_iwp*mbod(bod)%d1x1/dtim) -        &
                                  mbod(bod)%d2x1
        ! velocity
        mbod(bod)%kinup_d1x1    = 2.0_iwp*mbod(bod)%x1/dtim - mbod(bod)%d1x1
        mbod(bod)%kinup_d2x1(0) = zero
        mbod(bod)%kinup_d1x1(0) = zero
        mbod(bod)%d1p1(0)       = zero
    END DO 

    !---------------------------------------------------------------------------AS
    ! Stress recovery
    !---------------------------------------------------------------------------AS
    DO bod=2,size(mbod)
      k=1
      MatPoints: DO i=1,mbod(bod)%nmps
        iel   = mbod(bod)%a_ele(i)
        num   = mbod(bod)%g_num(:,iel)
        coord = TRANSPOSE(mbod(bod)%g_coord(:,num))
        g_s   = mbod(bod)%g_g(:,iel)

        CALL shape_der(der,mbod(bod)%mpoints,i)
        jac   = MATMUL(der,coord)
        det   = determinant(jac)
        CALL invert(jac)
        deriv = MATMUL(jac,der)
        CALL beemat(bee,deriv)
        
        k=k+1
        IF(k>nip)k=1
        
        ! calculate updated strains distribution, remember loads is displacement
        eld   = mbod(bod)%loads(g_s)
        eps_s = MATMUL(bee,eld)

        ! calculate the stress based on the elastic strain
        mbod(bod)%eps_acum(:,i) = mbod(bod)%eps_acum(:,i)+eps_s
        CALL deemat(mbod(bod)%dee,mbod(bod)%Young,mbod(bod)%Poiss)
        sigma = MATMUL(mbod(bod)%dee,eps_s)       

        mbod(bod)%m_stress_efe(:,i)=mbod(bod)%m_stress_efe(:,i)+sigma
        CALL invar(mbod(bod)%m_stress_efe(:,i),sigm,dsbar,lode_theta)

        ! calculate the stresses
        ps1 = sigm+(2.0/3.0)*dsbar*sin(lode_theta-(2.0*3.1415926/3.0))
        ps2 = sigm+(2.0/3.0)*dsbar*sin(lode_theta)
        ps3 = sigm+(2.0/3.0)*dsbar*sin(lode_theta+(2.0*3.1415926/3.0))
        mbod(bod)%Devstress(i)   = (1.0/sqrt(two))*sqrt((ps1-ps2)**2+&
                                   (ps2-ps3)**2+(ps3-ps1)**2)
        mbod(bod)%mean_stress(i) = (mbod(bod)%m_stress_efe(1,i) + &
                                   mbod(bod)%m_stress_efe(2,i) +  &
                                   mbod(bod)%m_stress_efe(4,i))/3.0_iwp
        mbod(bod)%m_stress(1,i)  = mbod(bod)%m_stress_efe(1,i)
        mbod(bod)%m_stress(2,i)  = mbod(bod)%m_stress_efe(2,i)
        mbod(bod)%m_stress(3,i)  = mbod(bod)%m_stress_efe(3,i)
        mbod(bod)%m_stress(4,i)  = mbod(bod)%m_stress_efe(4,i)
      END DO MatPoints
    END DO
    
    ! Limit iteration to 30 cycles (AS)
    IF(iters>30) EXIT 
  END DO FF_Newton_Rhapson

  
  !--- Determine freefield boundaries displacement [TODO!!]
  
  
  
  !--- Solve main MPM body stresses and displacement levels
  iters=0; iters2=0
  newrapconv1=.false.
  MPM_Newton_Rhapson: DO WHILE(iters<limit+1)   
    iters=iters+1
    
    IF(DEBUG) write(800, *) " "
    IF(DEBUG) write(800, '("Steps :" (I10) "/" (I10))'), iters, w
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"x1"',mbod(1)%x1
    
    
    !---------------------------------------------------------------------------AS
    ! Determine Nodal Kinematics
    !---------------------------------------------------------------------------AS
    
    ! calculate node velocity (d1x1) and acceleration (d2x1)
    DO bod=1,1
      IF(mbod(bod)%nmps>1)THEN 
        MP_VelAccP: DO i=1,mbod(bod)%nmps
      
          iel=mbod(bod)%a_ele(i)    
          IF(smethod==2.or.smethod==3)THEN 
            CALL GIMP_nodsup(i,g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%mpoints, &
              mbod(bod)%valuesg,mbod(bod)%gimptol,neighb,nny,mbod(bod)%a_ele,mbod(bod)%GIMP_nodes)
            values=mbod(bod)%valuesg(i)
            ALLOCATE(derextend(nodof,values),jac_coord(values,nodof),funextend2(values*2,2),	&
                  eldddylds(values*2),beeextend(nst,values*nodof))

            eldddylds=zero;funextend2=zero;beeextend=zero;derextend=zero;jac_coord=zero
   
            CALL GIMP_funder2(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)  !--creates the matrix with the values of shape functions and derivatives
            CALL gimpfunform(i,iel,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g,mvval)  !--collect in eldddylds the nuber of the equations in the support domain
                              
            mbod(bod)%m_mvp(eldddylds)=mbod(bod)%m_mvp(eldddylds)+MATMUL(MATMUL(funextend2,mbod(bod)%m_mass(i)*delta),mbod(bod)%m_velocity(:,i))
            mbod(bod)%m_mva(eldddylds)=mbod(bod)%m_mva(eldddylds)+MATMUL(MATMUL(funextend2,mbod(bod)%m_mass(i)*delta),mbod(bod)%accb(:,i))

            mbod(bod)%m_mvp(0)=zero;mbod(bod)%gravlo(0)=zero;mbod(bod)%m_mva(0)=zero
   
            DEALLOCATE(derextend,funextend2,eldddylds,jac_coord,beeextend)
          END IF
        END DO MP_VelAccP

        DO i=1,mbod(bod)%neq
          mbod(bod)%d1x1(i)=mbod(bod)%m_mvp(i)/mbod(bod)%diag(i) ! m_mvp & diag is the contribution from all bodies and at the contact d1x1 could be also zero
          mbod(bod)%d2x1(i)=mbod(bod)%m_mva(i)/mbod(bod)%diag(i)
        END DO
      END IF
    END DO 
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"d1x1"',mbod(1)%d1x1
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"d2x1"',mbod(1)%d2x1
    
    
    ! calculate ground acceleration (kinup_Ground_d2x1)
    DO bod=1,1
      IF(w<=accdata)THEN
        DO i=1,nn
          IF(g_coord(2,i)<lowbound+cellsize+0.01_iwp.and.mbod(bod)%nf(1,i)>0)THEN
            mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i))=mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i))+(mbod(bod)%ground_acc(w)-mbod(bod)%kinup_d2x1(mbod(bod)%nf(1,i)))
            mbod(bod)%kinup_Ground_d2x1(0)=zero
          END IF    
        END DO
      END IF
    END DO
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"ground_d2x1"',mbod(1)%kinup_Ground_d2x1
    
    
    ! Kinematic update of acceleration (kinup_d2x1) and velocity (kinup_d1x1) 
    ! without considering any boundary condition (contact)
    DO bod=1,1
      mbod(bod)%kinup_d2x1 = (4.0_iwp*mbod(bod)%x1/dtim**2) - (4.0_iwp*mbod(bod)%d1x1/dtim) - mbod(bod)%d2x1
      mbod(bod)%kinup_d2x1(0)=zero
    END DO 
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"kinup_d1x1"',mbod(1)%kinup_d1x1
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"kinup_d2x1"',mbod(1)%kinup_d2x1
    
    !---------------------------------------------------------------------------AS
    ! Determine Nodal Forces
    !---------------------------------------------------------------------------AS
    
    ! Determine the internal force (ddylds) from particle stresses
    DO bod=1,1
      IF(mbod(bod)%nmps>1)THEN    
        mbod(bod)%ddylds=zero

        MPM_Fint_load:DO i=1,mbod(bod)%nmps
          iel=mbod(bod)%a_ele(i)   
          num=g_num(:,iel)
          coord=TRANSPOSE(g_coord(:,num))
          g=mbod(bod)%g_g(:,iel)
          IF(smethod==2.or.smethod==3)THEN 
            values=mbod(bod)%valuesg(i)
            ALLOCATE(derextend(nodof,values),jac_coord(values,nodof),funextend2(values*2,2),eldddylds(values*2),beeextend(nst,values*nodof))
            
            eldddylds=zero;funextend2=zero;beeextend=zero;derextend=zero;jac_coord=zero
            CALL GIMP_funder2(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)  !--creates the matrix with the values of shape functions and derivatives
            CALL gimpfunform(i,iel,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g,mvval)  !--collect in eldddylds the nuber of the equations in the support domain
            CALL Sup_coord(i,mbod(bod)%GIMP_nodes,g_coord,jac_coord)
            jac=MATMUL(derextend,jac_coord)
            jac(2,1)=0.0;jac(1,2)=0.0
            CALL invert(jac)
            derextend=MATMUL(jac,derextend)
            CALL beematgimp(i,beeextend,derextend,values)

            mbod(bod)%ddylds(eldddylds)=mbod(bod)%ddylds(eldddylds) + MATMUL(mbod(bod)%m_stress(:,i), beeextend)*(4.0*mbod(bod)%lp_mp(1,i) * mbod(bod)%lp_mp(2,i))
                              
            DEALLOCATE(derextend,funextend2,eldddylds,jac_coord,beeextend)
          ELSE
            CALL shape_der(der,mbod(bod)%mpoints,i)
            jac=MATMUL(der,coord)
            det=determinant(jac)
            CALL invert(jac)
            deriv=MATMUL(jac,der)
            CALL beemat(bee,deriv)
            waverage=1.0_iwp
            
            mbod(bod)%ddylds(g)=mbod(bod)%ddylds(g) + MATMUL(mbod(bod)%m_stress(:,i),bee)*mbod(bod)%m_volume(i)
          
          END IF
          mbod(bod)%ddylds(0)=zero 
        END DO MPM_Fint_load

      END IF
    END DO 
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"ddylds"',mbod(1)%ddylds
    
  
    ! calculate f_earth = mv x kinup_Ground_d2x1
    MPM_Ground_F: DO bod=1,1
      IF(mbod(bod)%nmps>1)THEN 
        CALL linmul_sky(mbod(bod)%mv,mbod(bod)%kinup_Ground_d2x1,mbod(bod)%f_earth,mbod(bod)%kdiag)  !--Multiplication of the mass matrix per the acceleration vector vcm=Ma 
        mbod(bod)%f_earth(0)=zero
      END IF
    END DO MPM_Ground_F 
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"f_earth"',mbod(1)%f_earth


    ! calculate cdamp = mvkv x kinup_d1x1 = (fm*mv + fk*kv) x kinup_d1x1
    Damping_force: DO bod=1,1
      IF(mbod(bod)%nmps>1)THEN
        mbod(bod)%cdamp=zero
        mbod(bod)%mvkv=fm*mbod(bod)%mv+fk*mbod(bod)%kv
        CALL linmul_sky(mbod(bod)%mvkv,mbod(bod)%kinup_d1x1,mbod(bod)%cdamp,mbod(bod)%kdiag)
        mbod(bod)%cdamp(0)=zero
      END IF
    END DO Damping_force
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"cdamp"',mbod(1)%cdamp
    

    ! calculate vcm = mv x kinup d2x1
    MPM_Force_MA: DO bod=1,1
      IF(mbod(bod)%nmps>1)THEN 
        CALL linmul_sky(mbod(bod)%mv,mbod(bod)%kinup_d2x1,mbod(bod)%vcm,mbod(bod)%kdiag)  !--Multiplication of the mass matrix per the acceleration vector vcm=Ma 
        mbod(bod)%vcm(0)=zero   
      END IF
    END DO MPM_Force_MA   
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"vcm"',mbod(1)%vcm
    
    
    ! calculate gravity loading (gravlo)
    
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"gravlo"',mbod(1)%gravlo
    
    
    !---------------------------------------------------------------------------AS
    ! Solve the equations
    !---------------------------------------------------------------------------AS

    MPM_DISPLACEMENTS: DO bod=1,1
      mbod(bod)%loads=zero
      mbod(bod)%loads=mbod(bod)%gravlo-mbod(bod)%ddylds-mbod(bod)%vcm+mbod(bod)%f_earth-mbod(bod)%cdamp
      mbod(bod)%loads(0)=zero 
      IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"loads"',mbod(1)%loads
      
      CALL spabac(mbod(bod)%kp,mbod(bod)%loads,mbod(bod)%kdiag)  
      mbod(bod)%loads(0)=zero 
      mbod(bod)%x1=mbod(bod)%x1+mbod(bod)%loads
      mbod(bod)%x1(0)=zero
    END DO MPM_DISPLACEMENTS
    

    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"loads_post"',mbod(1)%loads
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"x1_post"',mbod(1)%x1
    
    !---------------------------------------------------------------------------AS
    ! Update lagrangian phase / update nodal kinematics to next iteration values
    !---------------------------------------------------------------------------AS

    DO bod=1,1
      mbod(bod)%kinup_d2x1=(4.0_iwp*mbod(bod)%x1/dtim**2)-(4.0_iwp*mbod(bod)%d1x1/dtim)-mbod(bod)%d2x1 
      mbod(bod)%kinup_d1x1=2.0_iwp*mbod(bod)%x1/dtim-mbod(bod)%d1x1     !--Upgrade lagrangian phase, in the paper (wang)pag.162: v=(2/Dt)u-v  
      mbod(bod)%kinup_d1x1(0)=zero;mbod(bod)%kinup_d2x1(0)=zero
    END DO 
    
    
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"kinup_d1x1_p"',mbod(1)%kinup_d1x1
    IF(DEBUG) write(800, '((A15":"), *(E15.5 ","))'), '"kinup_d2x1_p"',mbod(1)%kinup_d2x1
    
    ! check convergence
    DO bod=1,1
      CALL checon_1(mbod(bod)%loads,mbod(bod)%x1_change,newrapconv1,tol) 
      mbod(bod)%x1_change=mbod(bod)%x1
    END DO 

    !---------------------------------------------------------------------------AS
    ! Stress recovery
    !---------------------------------------------------------------------------AS
    
    Stress_rec: DO bod=1,1
      MPM_MatPoints: DO i=1,mbod(bod)%nmps
        iel=mbod(bod)%a_ele(i)
        num=g_num(:,iel)
        coord=TRANSPOSE(g_coord(:,num))
        g=mbod(bod)%g_g(:,iel)
     
        !- The following is the use of the CMPM technique to compute stresses using neighbouring element strains
        !- Note that "values" equal to 16 indicate that only elements surounded by elements with material points are considered
        !- therefore boundary elements are not considered, since they are not surounded by filled elements
        !- to consider also boundary elements and CMPM the IF condition should be "values>=9" and not "values==16"   
        IF(smethod>=2) CALL bound_detect(iel,neighb,bound,values,mbod(bod)%c_ele,mbod(bod)%nf,g_num,aver)
        IF((values==16).and.smethod>=2)THEN 
          average:DO m=1,aver
            IF(m==1)THEN
              ALLOCATE(derextend(nodof,values),funCMPM(values),equations(values,values), &
                identidad(values,values),beeCMPM(nst,values*nodof),eldCMPM(values*nodof))
            END IF
            derextend=zero;funCMPM=zero;equations=zero;identidad=zero
            CALL const(equations,values,bound,m)
            CALL inv(equations,values)
            CALL shapefunct(mbod(bod)%mpoints,i,equations,funCMPM)
            CALL derivatives(mbod(bod)%mpoints,i,equations,derextend)
            CALL shape_der(der,mbod(bod)%mpoints,i)
            jac=MATMUL(der,coord)
            CALL invert(jac)
            derextend=MATMUL(jac,derextend)
            CALL beemat3(beeCMPM,derextend,values,funCMPM)

            CALL eldform(eldCMPM,mbod(bod)%loads,iel,neighb,values,mbod(bod)%g_g,bound,m)

            IF(m==1) eps=MATMUL(beeCMPM,eldCMPM)
            IF(m==2) eps_2=MATMUL(beeCMPM,eldCMPM)

            !sigma=MATMUL(mbod(bod)%dee,eps_2)
            IF(m==1) sigma_1=MATMUL(mbod(bod)%dee,eps)
            !IF(m==2)sigma_2=MATMUL(mbod(bod)%dee,eps_2)
          END DO average
          !IF(aver==1)sigma=sigma_1
          IF(aver==1) eps=eps
          !IF(aver==2)sigma=sigma_1*0.5+sigma_2*0.5
          IF(aver==2) eps=0.5_iwp*eps+0.5_iwp*eps_2
          DEALLOCATE(derextend,funCMPM,equations,identidad,beeCMPM,eldCMPM)
        ELSE IF(smethod==2)THEN
          values=mbod(bod)%valuesg(i)
          ALLOCATE(derextend(nodof,values),funextend(values),beeextend(nst,values*nodof), &
            eldCMPM(values*nodof),jac_coord(values,nodof))
     
          CALL GIMP_funder(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend)
          CALL Sup_coord(i,mbod(bod)%GIMP_nodes,g_coord,jac_coord) 
          jac=MATMUL(derextend,jac_coord)
          jac(2,1)=zero;jac(1,2)=zero
          CALL invert(jac)
          derextend=MATMUL(jac,derextend)
          CALL beematgimp(i,beeextend,derextend,values)   
          CALL eldformgimp(i,eldCMPM,mbod(bod)%loads,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g)
          eps=MATMUL(beeextend,eldCMPM)
          DEALLOCATE(derextend,funextend,beeextend,eldCMPM,jac_coord) 
        ELSE  
          CALL shape_der(der,mbod(bod)%mpoints,i)
          jac=MATMUL(der,coord)
          jac(2,1)=zero;jac(1,2)=zero
          det=determinant(jac)
          CALL invert(jac)
          deriv=MATMUL(jac,der)
          eld=mbod(bod)%loads(g)
          CALL beemat(bee,deriv)
          waverage=1.0_iwp
          eps=MATMUL(bee,eld)*waverage
          sigma=MATMUL(mbod(bod)%dee,eps)      
        END IF
     
        ! Stress updated considering "eps" computed with standard FEM or CMPM
        sigma  = MATMUL(mbod(bod)%dee,eps)  
        stress = sigma + mbod(bod)%m_stress(:,i) !Find new stress level

        ! Check whether yield is violated
        CALL invar(stress,sigm,dsbar,lode_theta)
        fnew=dsbar-SQRT(d3)*mbod(bod)%mpcp(i)
        mbod(bod)%mpyield(i)=fnew
        IF(fnew>fmax) fmax=fnew
        epsinv=mbod(bod)%epsinvacum(i)
    
        ! This part is ignored due to `ploop = 0`
        Yieldf:IF(fnew>zero.and.ploop==1)THEN
          stress=mbod(bod)%m_stress(:,i)  !Stresses inside the yield function (elastic range)    
          CALL alphavalue2(mbod(bod)%mpcp(i),stress,sigma+stress,alpha_2)
          stress=stress+alpha_2*sigma    ! Stress just at the yield function boundary
          mbod(bod)%flag(i)=1

          !-----------performace substeping method (Bin Wang modifyed algoritm)----------
          !-----------Refined explicit integration of elastoplastic models with----------
          !----------------automatic error controlSloan, Abbo & Sheng 2001---------------
              
          !--NOTE: For some reason, Bi in the paper is 1 in the code
          fail=0
          Ai=zero
          Deps=eps*(1-alpha_2)  !--part of the strains that send stresses outside the yield function
          count=0
          Tret=zero
          DTret=one

          iterations1:DO WHILE(Tret<one)
            IF(count>9999) EXIT

            depst=DTret*Deps
            Dsigmae1=MATMUL(mbod(bod)%dee,DTret*Deps)                     !--Amount of stress outside the yield function d_sig
            CALL invar(stress,sigm,dsbar,lode_theta)                      !-Invariants considering the stresses at the yield function
            CALL vmflow(stress,dsbar,vmfl)                                !- Change of the yield function in terms of the change of stresses (dF/dsigma)
            flowf=vmfl
            flowg=flowf
            CALL plasticmod(epsinv,Maxplastic,SModulus,h)
            Ai=sqrt(d3)*h                                                 !--Should be negative if the soil is softening
            CALL formdlambda(flowf,flowg,mbod(bod)%dee,depst,Ai,dlambda)  !dlambda should be positive 
            dlambda=MAX(dlambda,zero)                                     !--plastic multiplier
            Dsigma1=Dsigmae1-dlambda*MATMUL(mbod(bod)%dee,flowg)          !--Dsigma1 is the stress again in the yield function (Dsigmae1 goes out and dlambda*MATMUL(dee,flowg) goes in again)
                                                                          !--Elasto-plastic stress rate dsig=(De-Dp)*deps
            dh1=dlambda 
        
            sigma2=stress+Dsigma1                                         ! sigma2 is the stress in the yield function again
            h2=epsinv+dh1                                                 ! k = epsinv
            Dsigmae2=MATMUL(mbod(bod)%dee,DTret*deps)
            CALL invar(sigma2,sigm,dsbar,lode_theta)

            CALL vmflow(sigma2,dsbar,vmfl)
            flowf=vmfl
            flowg=flowf
            CALL plasticmod(h2,Maxplastic,SModulus,h)
            Ai=sqrt(d3)*h
            CALL formdlambda(flowf,flowg,mbod(bod)%dee,depst,Ai,dlambda)
      
            dlambda=MAX(dlambda,zero)
            Dsigma2=Dsigmae2-dlambda*MATMUL(mbod(bod)%dee,flowg)
            dh2=dlambda
      
            !- compute the new stresses and hardening parameter add the rotation stress.
            sigma_trial=(Dsigma1+Dsigma2)/two+stress
            h_trial=epsinv+(dh1+dh2)/two
              
            IF(h_trial<0.0) h_trial=0.0

            !- determine the relative error for the current substep
            Dsigmanorm=zero
            Tsigmanorm=zero
            DO m=1,nst
              Dsigmanorm=Dsigmanorm+(dsigma2(m)-dsigma1(m))**2
              Tsigmanorm=Tsigmanorm+sigma_trial(m)**2
            END DO  
            Dsigmanorm=sqrt(Dsigmanorm)
            Tsigmanorm=sqrt(Tsigmanorm)
      
            IF(Dsigmanorm==0.0 .or. Tsigmanorm==0.0)THEN
              r1=0.0
            ELSE
              r1=Dsigmanorm/(two*Tsigmanorm)
            END IF
        
            IF(h_trial==0.0)THEN
              r2=0.0
            ELSE
              r2=(ABS(dh2-dh1))/(two*h_trial)
            END IF
      
            r_error=MAX(epstol,r1,r2) 

            Tol_loop:IF(r_error>stol)THEN    !the substep has failed a smaller step is needed
              qtol=MAX(0.90_iwp*SQRT(stol/r_error),0.10_iwp)
              DTret=MAX(qtol*DTret,dt_min)
              fail=1
              count=count+1
              CYCLE !--The CYCLE statement send me back again to the beguinning of the iterations1 loop
            ELSE                    
              count=count+1
              !this substep is acceped
              stress=sigma_trial
              epsinv=h_trial
              CALL invar(stress,sigm,dsbar,lode_theta)
              mbod(bod)%mpcp(i)=epsinv*SModulus+mbod(bod)%prop(4,1)
              IF(mbod(bod)%mpcp(i)<mpcr)                mbod(bod)%mpcp(i)=mpcr 
              IF(mbod(bod)%mpcp(i)>mbod(bod)%prop(4,1)) mbod(bod)%mpcp(i)=mbod(bod)%prop(4,1)
              CALL plasticmod(epsinv,Maxplastic,SModulus,h)                            
              Ai=sqrt(d3)*h
              f=dsbar-sqrt(d3)*mbod(bod)%mpcp(i)
        
              IF(mbod(bod)%mpcp(i)<mpcr)  mbod(bod)%mpcp(i)=mpcr
              IF(mbod(bod)%mpcp(i)>cpeak) mbod(bod)%mpcp(i)=cpeak

              Tret=Tret+DTret
      
              IF(ABS(f)>ftol)THEN
                CALL yieldcorrection(stress,epsinv,cpeak,0.0_iwp,cpeak,    &
                  Maxplastic,mpcr,0.0_iwp,SModulus,mbod(bod)%dee,ft)
              END IF
        
              qtol=MIN(0.90_iwp*SQRT(stol/r_error),1.10_iwp)
              IF(fail==1) THEN
                qtol=MIN(qtol,one)
                fail=0
              END IF
              DTret=qtol*DTret
              DTret=MAX(DTret,dt_min)
              DTret=MIN(DTret,one-Tret)
        
            END IF Tol_loop
          END DO iterations1

          IF(count>9999)EXIT
        END IF Yieldf
      
        ! Obtain stress invariants after plastic strain distributions
        mbod(bod)%m_stress(:,i)=stress
        CALL invar(stress,sigm,dsbar,lode_theta)
        ps1=sigm+(2.0/3.0)*dsbar*sin(lode_theta-(2.0*3.1415926/3.0))
        ps2=sigm+(2.0/3.0)*dsbar*sin(lode_theta)
        ps3=sigm+(2.0/3.0)*dsbar*sin(lode_theta+(2.0*3.1415926/3.0))
        mbod(bod)%Devstress(i)=(1.0/sqrt(two))*sqrt((ps1-ps2)**2+(ps2-ps3)**2+(ps3-ps1)**2)
        mbod(bod)%mean_stress(i)=(stress(1)+stress(2)+stress(4))/3.0_iwp
        mbod(bod)%eps_m(:,i)=mbod(bod)%eps_m(:,i)+eps
        mbod(bod)%epsinvacum(i)=epsinv
        IF((iters>limit.or.(newrapconv1)))THEN
          mbod(bod)%eps_acum(:,i)=mbod(bod)%eps_acum(:,i)+mbod(bod)%eps_m(:,i)
        END IF
      END DO MPM_MatPoints
    END DO Stress_rec   
 
    IF(iters>limit) EXIT
  END DO MPM_Newton_Rhapson

  
  !=============================================================================AS
  ! Post iteration processes
  !=============================================================================AS
  
  ! x1_acum is here just for post-simulation purpose. Not for actual calculation
  DO bod=2,size(mbod)
      mbod(bod)%x1_acum = mbod(bod)%x1_acum+mbod(bod)%x1
      mbod(bod)%d2x1    = mbod(bod)%kinup_d2x1
      mbod(bod)%d1x1    = mbod(bod)%kinup_d1x1
  END DO

  ! Update the node displacement
  !- NOTE: Changing nodes coordenates can damage simulations 
  DO bod=2,size(mbod)
      DO i=1,mbod(bod)%nn
          mbod(bod)%g_coord_aux(1,i) = mbod(bod)%g_coord_aux(1,i) +   &
                                      mbod(bod)%x1(mbod(bod)%nf(1,i))
          mbod(bod)%g_coord_aux(2,i) = mbod(bod)%g_coord_aux(2,i) +   &
                                      mbod(bod)%x1(mbod(bod)%nf(2,i))
      END DO
  END DO
  
  
  !=============================================================================AS
  ! Reverse Mapping from Nodes to Material Points
  !=============================================================================AS
  ! For MPM, this step is important as the following step will calculate 
  ! nodal velocity and nodal acceleration from the kinematics in the MPs
  ! 
  ! For FEM, this step is solely for visualization purposes
  
  cont2=cont2+1
  REVERSE_MAPPING: DO bod=1,size(mbod)
    IF(bod==1)THEN
      mbod(bod)%vccp=zero
      IF(mbod(bod)%nmps>1)THEN 
        DO i=1,mbod(bod)%nmps
          values=mbod(bod)%valuesg(i)
          ALLOCATE(derextend(nodof,values),beeextend(nst,values*nodof),funextend2(values*2,2),	&
            eldddylds(values*2),funextend(values),jac_coord(values,nodof),eldCMPM(values*2))
          CALL GIMP_funder(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend) !--creates the matrix with the values of shape functions and derivatives
        
          CALL eldformgimp(i,eldCMPM,mbod(bod)%kinup_d2x1,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g)
          mbod(bod)%accp(1,i)=DOT_PRODUCT(eldCMPM(1:ndof:2),funextend)
          mbod(bod)%accp(2,i)=DOT_PRODUCT(eldCMPM(2:ndof:2),funextend)
          mbod(bod)%m_acc(:,i)=mbod(bod)%accp(:,i)
          mbod(bod)%m_velocity(:,i)=mbod(bod)%m_velocity(:,i)+0.5_iwp*(mbod(bod)%accp(:,i)+mbod(bod)%accb(:,i))*dtim

          CALL eldformgimp(i,eldCMPM,mbod(bod)%x1,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g)
          mbod(bod)%ins(1,i)=DOT_PRODUCT(eldCMPM(1:ndof:2),funextend)
          mbod(bod)%ins(2,i)=DOT_PRODUCT(eldCMPM(2:ndof:2),funextend)
          mbod(bod)%a_ins(:,i)=mbod(bod)%a_ins(:,i)+mbod(bod)%ins(:,i)
          mbod(bod)%gm_coord(:,i)=mbod(bod)%gm_coord(:,i)+mbod(bod)%ins(:,i)
          DEALLOCATE(derextend,beeextend,funextend2,eldddylds,funextend,jac_coord,eldCMPM)
        END DO
        mbod(bod)%accb=mbod(bod)%accp
      END IF
    ELSE
      DO i=1,mbod(bod)%nmps
        iel=mbod(bod)%a_ele(i)
        mbod(bod)%g=mbod(bod)%g_g(:,iel)
        eld=mbod(bod)%kinup_d2x1(mbod(bod)%g)
        CALL shape_fun(fun,mbod(bod)%mpoints,i)
        !** Update material point velocity
        mbod(bod)%accp(1,i)  = DOT_PRODUCT(eld(1:ndof:2),fun)
        mbod(bod)%accp(2,i)  = DOT_PRODUCT(eld(2:ndof:2),fun)
        mbod(bod)%m_acc(:,i) = mbod(bod)%accp(:,i)
        mbod(bod)%m_velocity(:,i) = mbod(bod)%m_velocity(:,i) +     &
                                    0.5_iwp*(mbod(bod)%accp(:,i) +  &
                                    mbod(bod)%accb(:,i))*dtim       
        !** Convect material points
        eld=mbod(bod)%x1(mbod(bod)%g(1:ndof))
        mbod(bod)%ins(1,i)   = DOT_PRODUCT(eld(1:ndof:2),fun)  
        mbod(bod)%ins(2,i)   = DOT_PRODUCT(eld(2:ndof:2),fun) 
        mbod(bod)%a_ins(:,i) = mbod(bod)%a_ins(:,i) + mbod(bod)%ins(:,i)
        mbod(bod)%gm_coord(:,i) = mbod(bod)%gm_coord(:,i)+mbod(bod)%ins(:,i)
      END DO
      mbod(bod)%accb=mbod(bod)%accp
    END IF
  END DO REVERSE_MAPPING

  
  !=============================================================================AS
  !                          FREEFIELD BODIES CLEANUP
  !=============================================================================AS
  ! Clean up memory for Pardiso solver
  phase = -1
  DO bod=2,size(mbod)
    CALL pardiso(                                                   &
      mbod(bod)%pt,                                               &
      maxfct,                                                     &
      mnum,                                                       &
      mtype,                                                      &
      phase,                                                      &
      mbod(bod)%neq,                                              &
      mbod(bod)%kv_CSR,                                           &
      mbod(bod)%ia,                                               &
      mbod(bod)%ja,                                               &
      idum,                                                       &
      nrhs,                                                       &
      mbod(bod)%iparm,                                            &
      msglvl,                                                     &
      mbod(bod)%loads(1:mbod(bod)%neq),                           &
      b,                                                          &
      mbod(bod)%residual(1:mbod(bod)%neq),                        &
      error                                                       &
    )
    DEALLOCATE(mbod(bod)%ia,mbod(bod)%ja,mbod(bod)%kv_CSR,mbod(bod)%iparm)
  END DO
  
  
  !=============================================================================AS
  ! Simulation Output and Visualization
  !=============================================================================AS
 
  !-- Loop to save data from both bodies and print it in point_vis
  printcount = 1
  IF(cont2/printval*printval==cont2)THEN
    printcount = printcount + 1
    IF (printcount > 5) THEN
      printcount = 1
      CLOSE(800)
      OPEN(800,FILE='Output/mpcoord'//'.txt', status="replace")
    END IF
    PRINT '("Steps :" (I10) "/" (I10))', step, 100000
    DO bod=1,size(mbod)
      IF(mbod(bod)%nmps>1)THEN 
          IF(bod==1)THEN
            CALL point_viz2(cont2,bod,argv,mbod(bod)%gm_coord,mbod(bod)%m_stress,mbod(bod)%m_stress_change,   &
              mbod(bod)%epsinvacum,mbod(bod)%a_ins,mbod(bod)%Devstress,mbod(bod)%mean_stress,             &
              mbod(bod)%mpyield,mbod(bod)%mpcp,mbod(bod)%m_velocity,mbod(bod)%m_acc,mbod(bod)%nmps,nlen)
            CALL paraview2(cont2,bod,argv,g_coord,g_num,mbod(bod)%nf,nels,nod,nn,nlen,                    &
              mbod(bod)%diag,mbod(bod)%ddylds,mbod(bod)%d1x1,mbod(bod)%d2x1,                              &
              mbod(bod)%gravlo,mbod(bod)%loads,mbod(bod)%normal,mbod(bod)%fcont,mbod(bod)%kv,mbod(bod)%mv,&
              mbod(bod)%kdiag,mbod(bod)%vcm,mbod(bod)%f_fint)     
          ElSE
              CALL paraview(                                              &
                  cont2,                                                  &
                  bod,                                                    &
                  argv,                                                   &
                  (mbod(bod)%ny1+1),                                      &
                  mbod(bod)%ny1,                                          &
                  mbod(bod)%g_coord_aux,                                  &
                  mbod(bod)%g_num,                                        &
                  mbod(bod)%nf,                                           &
                  mbod(bod)%nels,                                         &
                  nod,                                                    &
                  mbod(bod)%nn,                                           &
                  nlen,                                                   &
                  mbod(bod)%diag,                                         &
                  (-1.0)*mbod(bod)%ddylds,                                &
                  mbod(bod)%kinup_d1x1,                                   &
                  mbod(bod)%kinup_d2x1,                                   &
                  mbod(bod)%gravlo,                                       &
                  mbod(bod)%x1_acum,                                      &
                  mbod(bod)%P_ext,                                        &
                  mbod(bod)%mv,                                           &
                  mbod(bod)%kdiag,                                        &
                  (-1.0)*mbod(bod)%mf_force,                              &
                  mbod(bod)%loads,                                        &
                  -mbod(bod)%c_damp,                                      &
                  mbod(bod)%kv_CSR_aux,                                   &
                  mbod(bod)%penpos,1                                      &
              )     
              CALL point_viz(                                             &
                  cont2,                                                  &
                  bod,                                                    &
                  argv,                                                   &
                  mbod(bod)%gm_coord,                                     &
                  mbod(bod)%m_stress_efe,                                 &
                  mbod(bod)%eps_acum,                                     &
                  mbod(bod)%statev(:,7),                                  &
                  mbod(bod)%a_ins,mbod(bod)%Devstress,                    &
                  mbod(bod)%mean_stress,                                  &
                  mbod(bod)%m_pore(1,:),                                  &
                  mbod(bod)%statev(:,50),                                 &
                  mbod(bod)%m_velocity,                                   &
                  mbod(bod)%m_acc,mbod(bod)%nmps,                         &
                  nlen,                                                   &
                  1                                                       &
              )
          END IF
      END IF
    END DO
  END IF

  
  !=============================================================================AS
  ! Relocate MPM body locations for the next time step
  !=============================================================================AS

  ! allocate b_ele again, considering how many MP still inside of the original grid
  right_boundary = 0
  left_boundary = MAXINT
  Flags_2: DO bod=1,1 
    mbod(bod)%k_ele=0
    mbod(bod)%b=0
    mbod(bod)%dp_ele=mbod(bod)%d_ele
    mbod(bod)%d_ele=0

    DO i=1, mbod(bod)%nmps
      ! Determine element number based on colx and rowy
      colx=mbod(bod)%gm_coord(1,i)/cellsize+1
      rowy=ABS(mbod(bod)%gm_coord(2,i))/cellsize+1.0
      ielloc=(rowy-1.0)*nx1+colx
      
      !*** AS's Thesis
      ! Check whether element is a boundary element
      IF(rowy-offset_y <= size(right_boundary))THEN
        IF(right_boundary(rowy-offset_y) < colx) right_boundary(rowy-offset_y) = colx
        IF(left_boundary(rowy-offset_y) > colx) left_boundary(rowy-offset_y) = colx
      END IF
      
      ! Determine particle local coordinates w.r.t its element ielloc
      mbod(bod)%a_ele(i)=ielloc
      num=g_num(:,ielloc)
      coord=TRANSPOSE(g_coord(:,num))
      sp_coord(:,1)=mbod(bod)%gm_coord(:,i)
      lp_coord(1,:)=mbod(bod)%mpoints(i,:)
      CALL floc(coord,sp_coord,lp_coord,i)
      mbod(bod)%mpoints(i,:)=lp_coord(1,:)
      mbod(bod)%d_ele(ielloc)=1
    END DO

    mbod(bod)%activate=.false.
    DO i=1,nels
      IF(mbod(bod)%dp_ele(i)/=mbod(bod)%d_ele(i)) THEN
        mbod(bod)%activate=.true.
        EXIT
      END IF
    END DO

    CALL couma(nels,mbod(bod)%nmps,mbod(bod)%a_ele,mbod(bod)%c_ele,mbod(bod)%k_ele,etype)
  END DO Flags_2
  
  !*** AS's Thesis
  ! Check boundary element integrity and determine actual element number
  DO i=1,size(right_boundary)
    IF( right_boundary(i)<1 ) print *, "corresponding right boundary element", i, "is not found"
    right_boundary(i) = ((i+offset_y) - 1.0)*nx1 + right_boundary(i)
  END DO
  DO i=1,size(left_boundary)
    IF( left_boundary(i)>MAXINT-1 ) print *, "corresponding left boundary element", i, "is not found"
    left_boundary(i) = ((i+offset_y) - 1.0)*nx1 + left_boundary(i)
  END DO
  
  
  !=============================================================================AS
  ! Recreate global nodes and equations numbering
  !=============================================================================AS
  g_g=0; nf=0
  ChangeEl: DO bod=1,1
    !- Free all nodes associated with filled element,freeze the freedom of the empty element
    mbod(bod)%g_g=0
    mbod(bod)%nf=0
    IF(smethod==1)THEN
      DO iel=1,nels
        IF(mbod(bod)%c_ele(iel)>0)THEN
          num=g_num(:,iel)
          DO i=1,nod
            mbod(bod)%nf(:,num(i))=1
            nf(:,num(i))=1
          END DO
        END IF 
      END DO
    ELSE
      !-------------------------------------------------------------------------AS
      ! Find GIMP active nodes/elements
      !-------------------------------------------------------------------------AS
      mbod(bod)%GIMP_node_mp=0
      CALL GIMP_activenode(g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,  &
        mbod(bod)%gimptol,neighb,mbod(bod)%a_ele,mbod(bod)%nf,mbod(bod)%GIMP_node_mp)
      CALL Elem_suport(mbod(bod)%c_ele,mbod(bod)%nf,g_num,cellsize,nx1,nip,nels,  &
        g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%gimptol,smethod,mbod(bod)%elemmpoints)
 
      lowbound=minval(g_coord)
      upbound=maxval(g_coord)

      !-------------------------------------------------------------------------AS
      ! Locate boundary nodes
      !-------------------------------------------------------------------------AS
      j=1; m=1
      mbod(bod)%boundnod=0
      IF(bod==1)THEN
        DO i=1,nn
          !-Body fixities  
          IF(g_coord(2,i)<lowbound+cellsize+0.01) mbod(bod)%nf(2,i)=0
          IF(g_coord(2,i)<lowbound+cellsize+0.01)THEN
              mbod(bod)%boundnod(m)=i
              m=m+1
          END IF  
          !-Global fixities
          IF(g_coord(1,i)<cellsize+0.01_iwp)nf(1,i)=0
          IF(g_coord(1,i)<0.01_iwp)nf(:,i)=0
          IF(g_coord(2,i)<lowbound+0.01+cellsize)nf(:,i)=0  
        END DO
      END IF

      !-------------------------------------------------------------------------AS
      ! Determine equation numbering nf and reallocate memory
      !-------------------------------------------------------------------------AS
      CALL formnf(mbod(bod)%nf)
      mbod(bod)%neq=MAXVAL(mbod(bod)%nf)

      DEALLOCATE(mbod(bod)%d1x1,mbod(bod)%d2x1,mbod(bod)%diag,mbod(bod)%m_mvp,&
        mbod(bod)%loads,mbod(bod)%m_mva,mbod(bod)%ddylds,mbod(bod)%gravlo,    &
        mbod(bod)%fdamp,mbod(bod)%Freact,mbod(bod)%temp_d2x1,                 &
        mbod(bod)%fnorm,mbod(bod)%fcont,mbod(bod)%temp_d1x1,mbod(bod)%kdiag,  &
        mbod(bod)%vcm,mbod(bod)%x1,mbod(bod)%kinup_d2x1,mbod(bod)%kinup_d1x1, &
        mbod(bod)%x1_orig,mbod(bod)%f_fint,                                   &
        mbod(bod)%vel_change,mbod(bod)%x1_ini,                                &
        mbod(bod)%x1_change,mbod(bod)%f_earth,mbod(bod)%kinup_Ground_d2x1,    &
        mbod(bod)%cdamp,mbod(bod)%mvkv)

      ALLOCATE(mbod(bod)%d1x1(0:mbod(bod)%neq),mbod(bod)%d2x1(0:mbod(bod)%neq), &
        mbod(bod)%diag(0:mbod(bod)%neq),mbod(bod)%ddylds(0:mbod(bod)%neq),      &
        mbod(bod)%m_mvp(0:mbod(bod)%neq),mbod(bod)%loads(0:mbod(bod)%neq),      &
        mbod(bod)%fdamp(0:mbod(bod)%neq),mbod(bod)%m_mva(0:mbod(bod)%neq),      &
        mbod(bod)%Freact(0:mbod(bod)%neq),mbod(bod)%gravlo(0:mbod(bod)%neq),    &
        mbod(bod)%fnorm(0:mbod(bod)%neq),mbod(bod)%temp_d2x1(0:mbod(bod)%neq),  &
        mbod(bod)%fcont(0:mbod(bod)%neq),mbod(bod)%temp_d1x1(0:mbod(bod)%neq),  &
        mbod(bod)%kdiag(mbod(bod)%neq),mbod(bod)%vcm(0:mbod(bod)%neq),          &
        mbod(bod)%x1(0:mbod(bod)%neq),mbod(bod)%kinup_d2x1(0:mbod(bod)%neq),    &
        mbod(bod)%kinup_d1x1(0:mbod(bod)%neq),                                  &
        mbod(bod)%x1_orig(0:mbod(bod)%neq),                                     &
        mbod(bod)%f_fint(0:mbod(bod)%neq),                                      &
        mbod(bod)%vel_change(0:mbod(bod)%neq),                                  & 
        mbod(bod)%x1_ini(0:mbod(bod)%neq),mbod(bod)%x1_change(0:mbod(bod)%neq), &
        mbod(bod)%f_earth(0:mbod(bod)%neq),mbod(bod)%kinup_Ground_d2x1(0:mbod(bod)%neq), &
        mbod(bod)%cdamp(0:mbod(bod)%neq),mbod(bod)%mvkv(0:mbod(bod)%neq))
    
      !-------------------------------------------------------------------------AS
      ! Determine the shape and reallocate the skyline sparse matrix
      !-------------------------------------------------------------------------AS
      mbod(bod)%kdiag=zero
      DO iel=1,nels
        num=g_num(:,iel)
        CALL num_to_g(num,mbod(bod)%nf,g)
        CALL fkdiag(mbod(bod)%kdiag,g)
        mbod(bod)%g_g(:,iel)=g  
      END DO
    
      DO i=2,mbod(bod)%neq
        mbod(bod)%kdiag(i)=mbod(bod)%kdiag(i)+mbod(bod)%kdiag(i-1)
      END DO
  
      DEALLOCATE(mbod(bod)%kv,mbod(bod)%mv,mbod(bod)%kp)
      ALLOCATE(                                       &
        mbod(bod)%kv(mbod(bod)%kdiag(mbod(bod)%neq)), &
        mbod(bod)%mv(mbod(bod)%kdiag(mbod(bod)%neq)), &
        mbod(bod)%kp(mbod(bod)%kdiag(mbod(bod)%neq)))

    END IF
  END DO ChangeEl
  

  END DO time_steps
STOP
END PROGRAM Implicit_MPM_eartquake

