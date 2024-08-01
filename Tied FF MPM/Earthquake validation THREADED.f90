PROGRAM Implicit_MPM_eartquake
  USE mpm_main
  USE mpm_geom
  USE mpm_mpm
  USE mpm_gimp
  USE fem
  USE sparse_lib
  USE utilities
  IMPLICIT NONE
  !** Counters and placeholders
  INTEGER:: i,j,k,w,m,n,q,s,iel,iters,nsrf,ale,limit,bod,plc,printval,printcount, &
            colx,rowy,ielloc,bod2,bod3,iel2,nnx,nny,tcel,step,iters2,noddir,    &
            itersstep,accdata

  INTEGER:: nnz,nnzacum,accval,eqno,substep,tstep,   &
            sstep,v_presc,h_presc
  
  LOGICAL :: DEBUG=.true.

  !** Element/mesh properties
  INTEGER:: ndim=2,ndof=8,nels,neq,nip=4,nn,nod=4,                              &
    nodof=2,nprops=7,np_types,nst=4,nx1,nx2,nye,ny1,ny2,nlen,                   &
    row,column,slope1,slope2,newnodes,nxe,nipslo,                               &
    cont1,cont2,nnxe,nnye,eq1,eq2,dist_y,aiters,count,dist_x,elcont,mvval  
    
  REAL(iwp)::h1,h2,s1,w1,w2,lowbound,maxbound,leftbound,rightbound,             &
              cellsize,Gravf,HSURF,waverage,fm,fk,upbound
  
  LOGICAL:: shape=.false.,slope=.false.,equilibrium=.false.,newrapconv1, use_damping=.true.
  
  !*** AS's Thesis
  ! These variables will become placeholders for tracking boundary elements
  INTEGER,ALLOCATABLE:: left_boundary(:), right_boundary(:), iel_boundary(:,:)
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
    flowg(:),Deps(:),sigma2(:),flowf(:),epsinv_acum(:),depst(:),                &
    Devstress(:),deriv_gauss(:,:),bee_gauss(:,:),gauss_stress(:),               &
    scalar_stress(:),der_gauss(:,:),fun_gauss(:),sigma_trial(:),                &
    eldGIMP(:),sp_coord(:,:),lp_coord(:,:),km_gauss(:,:),mm_gimp(:,:),          &
    km(:,:),vel_change(:),eps_2(:),props(:)

  REAL(iwp),ALLOCATABLE:: mm(:,:),mm_s(:,:),kv_CSR(:),kv_CSR_aux(:),eps_s(:),   &
    km_mf(:,:),mm_mf(:,:),mm_mf_acum(:,:)
  
  INTEGER:: maxfct, mnum, mtype, phase, nrhs, error, msglvl
  INTEGER:: nels_bar,nip_1D=3,ndof_1D=6,nband,bandval
  REAL(iwp),ALLOCATABLE::kd(:),mm_acum(:,:),mvis(:),Cc(:,:)
  REAL(iwp) :: dparm(64)
  REAL(iwp) :: ddum,valone
  INTEGER   :: idum(1)
 
  !** Define the MPM objects to be used in the simulation
  TYPE(mpm_body):: mbod(3) ! data structure is defined in utilities

  ! Additional variables
  INTEGER::nyp,emps,nmps,tep,newel,slopeel,A,Bi
  INTEGER::nshr,nstatv=100,npt,layer,kspt,kstep,kinc,npropsum,void
  INTEGER,PARAMETER::ndi=3,nstumat=6

  INTEGER,ALLOCATABLE:: sum_vol(:),c_ele(:),m_num(:),a_ele(:),                  &
    k_ele(:),b_ele(:),d_ele(:),dp_ele(:),flag(:),g_s(:)

  INTEGER,ALLOCATABLE::penpos_count(:)
  REAL(iwp),ALLOCATABLE::penpos_displacement(:),penalized_stiffness(:)
  
  REAL(iwp),ALLOCATABLE::ini_volume(:),gm_coord(:,:),m_coord(:,:),              &
    mweights(:),m_mass(:),m_volume(:),ini_density(:),m_stress(:,:),             &
    m_velocity(:,:),ins_acum(:,:),nod_stress(:,:),d1x1(:,:),                    &
    diag(:,:),d2x1(:),m_mvp(:,:),m_mva(:),ddylds(:),                            &
    accp(:,:),acc(:),vcc(:),vccp(:,:),ins(:,:),stress(:),ecm_acum(:,:),ecm(:,:)

  REAL(iwp),ALLOCATABLE::ddsddt(:),drplde(:),drpldt(:),stran(:),                &
    time(:),predef(:),dpred(:),drot(:,:),dfgrd0(:,:),dfgrd1(:,:),               &
    stressumat(:),epsumat(:),deeumat(:,:),statev(:,:)

  
  !=============================================================================
  ! Input and Initialisation                              
  !=============================================================================
    
  nlen=7
  argv='Results'
  OPEN(800,FILE='Output/mpm_disp.dat')
  OPEN(810,FILE='Output/mpm_vel.dat')
  OPEN(820,FILE='Output/mpm_acc.dat')
  OPEN(830,FILE='Output/ff_disp.dat')
  OPEN(840,FILE='Output/ff_vel.dat')
  OPEN(850,FILE='Output/ff_acc.dat')
  
  OPEN(10,FILE='Input/Benchmark/Datafound.dat',status='old')
  OPEN(300,FILE='Input/Benchmark/Groundacc.dat',status='old')
  
  ! read material properties
  OPEN(400,FILE='Input/Benchmark/parumat.dat',status='old')
  READ(400,*)npropsum
  ALLOCATE(                                                 &
    ddsddt(nstumat),drplde(nstumat),stran(nstumat),         &
    props(npropsum),stressumat(nstumat),epsumat(nstumat),   &
    deeumat(nstumat,nstumat)                                &
  )
  READ(400,*)props

  ! Read parameter for the MPM body (bod=1) [AS]
  bod=1          
  mbod(1)%nprops=6 !For Von Mises with softening == 7
  mbod(1)%ndof=ndof
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
    mbod(bod)%ndof=mbod(1)%ndof

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
      mbod(bod)%slopeel=mbod(bod)%slopeel+i  !Slopeel is the number of elements in the slope / Part of the slope
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
    mbod(bod)%skylength=500*mbod(bod)%nels ! for freefield bodies
    
    ! Allocate array dimension dictated by nels and nn
    ALLOCATE(                                           &
      mbod(bod)%g_coord(ndim,mbod(bod)%nn),             &
      mbod(bod)%g_coord_aux(ndim,mbod(bod)%nn),         &  
      mbod(bod)%g_num(nod,mbod(bod)%nels),              &
      mbod(bod)%g_g(ndof,mbod(bod)%nels),               &
      mbod(bod)%g(mbod(bod)%ndof),                      &  
      mbod(bod)%nf(nodof,mbod(bod)%nn),                 &
      mbod(bod)%deeinv(nst,nst),                        &
      mbod(bod)%dee(nst,nst),                           &
      mbod(bod)%ini_volume(mbod(bod)%nels),             &
      mbod(bod)%c_ele(mbod(bod)%nels))                    
      
    ! Tied-FF Variables
    ALLOCATE(                                           &
      mbod(bod)%tied_nn(mbod(bod)%nn,2),                & 
      mbod(bod)%ia_aux(mbod(bod)%skylength),            &  
      mbod(bod)%ja_aux(mbod(bod)%skylength),            & 
      mbod(bod)%km(mbod(bod)%ndof,mbod(bod)%ndof),      & 
      mbod(bod)%KGC(mbod(bod)%ndof,mbod(bod)%ndof),     & 
      mbod(bod)%MMS(mbod(bod)%ndof,mbod(bod)%ndof),     & 
      mbod(bod)%CCS(mbod(bod)%ndof,mbod(bod)%ndof),     & 
      mbod(bod)%MOD_MTX(mbod(bod)%ndof,mbod(bod)%ndof)  & 
    )

    !- IF mbod(bod)%kconst=1, then the stiffness in the full domain will be constant
    !- any other value avoid using a constant stiffness
    ALLOCATE(mbod(bod)%kconst(1))
    mbod(bod)%kconst=0
  END DO count_nodes

  
  !===========================================================================AS
  ! GLOBAL VARIABLE ALLOCATIONS (Based on element properties)
  !===========================================================================AS
  ! nip = 4; nod = 4; nodof = 2; ndof = nod*nodof = 8; ndim = 2; nst = 4; 
  
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
    
  ALLOCATE(m_num(nip),vcc(ndof),m_coord(nip,ndim),gc(ndim),ecm(ndof,ndof),ecm_acum(ndof,ndof))
  ALLOCATE(gp(1,2))
  
  ALLOCATE(mm(ndof,ndof),mm_s(ndof,ndof)) ! mass matrices

  ALLOCATE(                           &
    acc(ndof),                        &
    delta(ndim,ndim),                 &
    delta_1(ndim),                    &
    g_matrix(ndim)) ! gravity field 
  
  !===========================================================================AS
  ! Geometry Generation
  !===========================================================================AS
  
  Mesh_create: DO bod=1,size(mbod)
    !-------------------------------------------------------------------------AS
    ! Determine Node Locations and Generate Mesh (nn, nels, nmps are known)
    !-------------------------------------------------------------------------AS

    !- Set up the global nf
    mbod(bod)%nf=0
    CALL emb_2d_bc2(mbod(bod)%nex,0,mbod(bod)%ney,mbod(bod)%newel,mbod(bod)%nf)

    !- Loop the elements to determine the global node and element numbering
    row=1; column=1
    A=0; Bi=0; dist_x=0
    DO iel=1,mbod(bod)%nels
      CALL emb_2d_geom2(iel,bod,mbod(bod)%nex,mbod(bod)%ney,mbod(bod)%s1,mbod(bod)%newel,       &
        mbod(bod)%slopeel,row,column,mbod(bod)%slope1,mbod(bod)%w1,mbod(bod)%h1,coord,num,A,Bi, &
        mbod(bod)%dist_x,mbod(bod)%dist_y,mbod(bod)%slopeopt)

      CALL num_to_g(num,mbod(bod)%nf,g)
      mbod(bod)%g_num(:,iel)=num
      mbod(bod)%g_coord(:,num)=TRANSPOSE(coord)
      mbod(bod)%g_g(:,iel)=g
    END DO

    
    !-------------------------------------------------------------------------AS
    ! Allocate Material Points Variables
    !-------------------------------------------------------------------------AS

    !- ALlocate material points tracers
    ALLOCATE(                                        &
      mbod(bod)%gm_coord(ndim,mbod(bod)%nmps),       &
      mbod(bod)%mpoints(mbod(bod)%nmps,ndim),        &
      mbod(bod)%a_ele(mbod(bod)%nmps),               &
      mbod(bod)%b(mbod(bod)%nmps)                    &
    )
      
    !- Allocate material points kinematics
    ALLOCATE(                                        &
      mbod(bod)%ins_acum(ndim,mbod(bod)%nmps),       &
      mbod(bod)%ins(ndim,mbod(bod)%nmps),            &
      mbod(bod)%accp(ndim,mbod(bod)%nmps),           &
      mbod(bod)%accb(ndim,mbod(bod)%nmps),           &
      mbod(bod)%vccp(ndim,mbod(bod)%nmps)            &
    )
        
    !- Allocate material points state variables
    ALLOCATE(                                        &
      mbod(bod)%ini_density(mbod(bod)%nmps),         &
      mbod(bod)%m_dens(mbod(bod)%nmps),              &
      mbod(bod)%mp_dens(mbod(bod)%nmps),             &
      mbod(bod)%m_volume(mbod(bod)%nmps),            &
      mbod(bod)%m_velocity(nodof,mbod(bod)%nmps),    &
      mbod(bod)%m_mass(mbod(bod)%nmps),              &
      mbod(bod)%m_stress(nst,mbod(bod)%nmps),        &
      mbod(bod)%m_acc(nodof,mbod(bod)%nmps),         &
      mbod(bod)%m_stress_ini(nst,mbod(bod)%nmps),    &
      mbod(bod)%eps_acum(nst,mbod(bod)%nmps),        &
      mbod(bod)%epsinv_acum(mbod(bod)%nmps)          &
    )
      
    !- Allocate additional MP state variables    
    ALLOCATE(                                        &
      mbod(bod)%eps_2(nst,mbod(bod)%nmps)            &
    )
        
    !- Allocate variables for yield functions
    ALLOCATE(                                        &
      mbod(bod)%flag(mbod(bod)%nmps),                &
      mbod(bod)%Devstress(mbod(bod)%nmps),           &
      mbod(bod)%mean_stress(mbod(bod)%nmps),         &
      mbod(bod)%mpyield(mbod(bod)%nmps)              &
    )
    
    ALLOCATE(                                        &
      mbod(bod)%m_pore(1,mbod(bod)%nmps),            &
      mbod(bod)%m_stress_efe(nst,mbod(bod)%nmps),    &
      mbod(bod)%statev(mbod(bod)%nmps,100),          &
      mbod(bod)%mpcp(mbod(bod)%nmps)                 &
    )
    
    !- Allocate GIMP variables
    ALLOCATE(                                        &
      mbod(bod)%lp_mp(ndim,mbod(bod)%nmps),          &
      mbod(bod)%valuesg(mbod(bod)%nmps),             &
      mbod(bod)%GIMP_nodes(9,mbod(bod)%nmps),        & ! maximum number of affected NODES under a particle is 9
      mbod(bod)%elemmpoints(mbod(bod)%nmps,4)        & ! maximum number of affected CELLS under a particle is 4
    )
    
      
    !-------------------------------------------------------------------------AS
    ! Insert Material Points / Stress Points
    !-------------------------------------------------------------------------AS
    
    ALLOCATE(mbod(bod)%MPPE(nip))
    mbod(bod)%a_ele=0
    mbod(bod)%c_ele=0
    mbod(bod)%MPPE=0
    DO i=1,nip
      mbod(bod)%MPPE(i)=i
    END DO

    IF(bod==1)THEN
      CALL sample2(element,points,weights)  
    ELSE
      CALL sample(element,points,weights)
    END IF
    
    InsertMP: DO iel=1,mbod(bod)%nels
      num=mbod(bod)%g_num(:,iel)
      coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
      DO i=1,nip
        CALL shape_fun(fun,points,i)
        m_coord(i,:)=MATMUL(fun,coord)
      END DO
      m_num=mbod(bod)%MPPE
      mbod(bod)%gm_coord(:,m_num)=TRANSPOSE(m_coord)
      mbod(bod)%a_ele(m_num)=iel
      mbod(bod)%c_ele(iel)=mbod(bod)%emps
      DO j=1,nip
        mbod(bod)%MPPE(j)=mbod(bod)%MPPE(j)+nip
      END DO
    END DO InsertMP

    
    !-------------------------------------------------------------------------AS
    ! Calculate initial volume, density and mass of each MPs
    !-------------------------------------------------------------------------AS
    
    !- Compute the initial volume of the CELLS
    Initialvol: DO iel=1,mbod(bod)%nels
      num=mbod(bod)%g_num(:,iel)
      coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
      mbod(bod)%ini_volume(iel)=cellsize*cellsize
    END DO Initialvol
    
    !- Compute the intiial volume of the PARTICLES & density & mass 
    k=1
    DO i=1,mbod(bod)%nmps
      iel=mbod(bod)%a_ele(i)

      ! Calculate m_volume using Gaussian Weights (FEM) or Cell Volume (MPM)
      IF(bod==1)THEN
        mbod(bod)%m_volume(i)=mbod(bod)%ini_volume(iel)/mbod(bod)%emps 
      ELSE
        CALL shape_der(der,points,k)
        jac=MATMUL(der,coord)
        det=determinant(jac)
        mbod(bod)%m_volume(i)=det*weights(k) 
        k=k+1; IF(k>nip) k=1
      END IF
      
      mbod(bod)%ini_density(i)=mbod(bod)%prop(1,1)
      mbod(bod)%mp_dens(i)=(mbod(bod)%prop(1,1)+0.75)/(1.75) !- rho=(Gs+e)/(1+e); Note that this equation assume a void ratio (e) of 0.75
      mbod(bod)%m_mass(i)=mbod(bod)%mp_dens(i)*mbod(bod)%m_volume(i)
    END DO
    
    mbod(bod)%m_dens=mbod(bod)%mp_dens !- For some reasons, the freefield use m_dens instead mp_dens
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
  offset_y = 1; offset_x=4
  nx1=dist_x + offset_x ! add 30 columns of elements to the right of the domain
  ny1=dist_y + offset_y ! add 1 row of elements at the bottom of the domain
  w1=cellsize*nx1 ! Size of domain in x-dir in units of distance
  h1=cellsize*ny1 ! Size of domain in y-dir in units of distance
  nx2=0; ny2=0; s1=0 ! no additional elements and no slope
  
  ! Calculate number of elements and number of nodes (nn & nels)
  nels=nx1*ny1
  nn=(ny1+1)*(nx1+1)
  nnxe=nx1
  nnye=ny1
  nnx=nx1+1
  nny=ny1+1

  !- Reallocate Global Variable with new nn & nels
  ALLOCATE(                           &
    g_coord(ndim,nn),                 &
    g_num(nod,nels),                  &
    g_g(ndof,nels),                   &
    nf(nodof,nn),                     &
    c_ele(nels),                      &
    d_ele(nels),                      &
    k_ele(0:nels),                    &
    etype(nels),                      &
    dp_ele(nels),                     &
    nod_stress(nst,nn),               &
    neighb(nels,8))

  ALLOCATE(m_mvp(0:nn,ndim), diag(0:nn,ndim), d1x1(0:nn,ndim))
    
  !- Allocate MPM body variable that tied with the background mesh
  DO bod=1,1
    ALLOCATE(mbod(bod)%g_g(ndof,nels))      
    ALLOCATE(mbod(bod)%nf(nodof,nn))        
    ALLOCATE(mbod(bod)%c_ele(nels))         
    ALLOCATE(mbod(bod)%d_ele(nels))         
    ALLOCATE(mbod(bod)%k_ele(0:nels))       
    ALLOCATE(mbod(bod)%GIMP_node_mp(mbod(bod)%nmps,nn))
    ALLOCATE(mbod(bod)%nodecont(nn))        
    ALLOCATE(mbod(bod)%etype(nels))         
    ALLOCATE(mbod(bod)%dp_ele(nels))        
    ALLOCATE(mbod(bod)%v_field(0:nn,ndim))  
    ALLOCATE(mbod(bod)%a_field(0:nn,ndim))  
    ALLOCATE(mbod(bod)%m_field(0:nn))       
    ALLOCATE(mbod(bod)%g_matrix(ndim))      
    ALLOCATE(mbod(bod)%tangent(nodof,nn))   
    ALLOCATE(mbod(bod)%boundnod(nn))        
    ALLOCATE(mbod(bod)%base_nn(nn))
  END DO
  
  ALLOCATE(iel_boundary(2, nels))

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
  right_boundary = 0
  left_boundary = MAXINT
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
      
      !*** AS's Thesis
      ! Check whether element is a boundary element
      IF(rowy - mbod(bod)%dist_y > 0)THEN
        IF(right_boundary(rowy-mbod(bod)%dist_y) < colx) right_boundary(rowy-mbod(bod)%dist_y) = colx
        IF(left_boundary(rowy-mbod(bod)%dist_y) > colx) left_boundary(rowy-mbod(bod)%dist_y) = colx
      END IF

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
  
  !*** AS's Thesis
  ! Check boundary element integrity and determine actual element number
  DO i=1,size(right_boundary)
    IF( right_boundary(i)<1 ) print *, "corresponding right boundary element", i, "is not found"
    right_boundary(i) = ((i+mbod(1)%dist_y) - 1.0)*nx1 + right_boundary(i)
  END DO
  DO i=1,size(left_boundary)
    IF( left_boundary(i)>MAXINT-1 ) print *, "corresponding left boundary element", i, "is not found"
    left_boundary(i) = ((i+mbod(1)%dist_y) - 1.0)*nx1 + left_boundary(i)
  END DO
  ! List all the boundary MPM cells and its corresponding freefield cells
  iel_boundary=0; k=1
  DO i=1,size(left_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = left_boundary(i)
    k = k+1
  END DO
  DO i=1,size(left_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = left_boundary(i) - 1
    k = k+1
  END DO
  DO i=1,size(right_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = right_boundary(i)
    k = k+1
  END DO
  DO i=1,size(right_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = right_boundary(i) + 1
    k = k+1
  END DO
  DO j=size(left_boundary),size(left_boundary)
    DO i=left_boundary(j), right_boundary(j)
      iel_boundary(1,k) = j
      iel_boundary(2,k) = i
      k = k+1
    END DO
  END DO

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
  ! Setup boundary conditions (determine nf and calculate neq)
  !===========================================================================AS
  
  nf=0; g_g=0
  MPM_fixities: DO bod=1,1
    mbod(bod)%g_g=0; mbod(bod)%nf=0
    IF(nip==4) mbod(bod)%lp_mp=cellsize/4.0   
    IF(nip==9) mbod(bod)%lp_mp=cellsize/6.0   
    
    ! mbod(bod)%nf is marked to 1 for active elements in Elem_suport()
    CALL GIMP_activenode(g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp, &
      mbod(bod)%gimptol,neighb,mbod(bod)%a_ele,mbod(bod)%nf,mbod(bod)%GIMP_node_mp)
    CALL Elem_suport(mbod(bod)%c_ele,mbod(bod)%nf,g_num,cellsize,nx1,nip,nels,  &
      g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%gimptol,smethod,mbod(bod)%elemmpoints)

    lowbound=minval(g_coord(2,:))
    upbound=maxval(g_coord(2,:))

    ! Set equation marker for zero-displacement boundary (y-dir on bottom elements)
    DO i=1,nn
      IF(g_coord(2,i)<lowbound+cellsize+0.01) mbod(bod)%nf(2,i)=0 
    END DO

    ! forming nf, and determine number of equations
    CALL formnf(mbod(bod)%nf)
    mbod(bod)%neq=MAXVAL(mbod(bod)%nf)
  END DO MPM_fixities

  Freefield_fixities: DO bod=2,size(mbod)
    ! Find the boundary of the geometry from g_coord
    leftbound = minval(mbod(bod)%g_coord(1,:))
    rightbound= maxval(mbod(bod)%g_coord(1,:))
    lowbound  = minval(mbod(bod)%g_coord(2,:))
    upbound   = maxval(mbod(bod)%g_coord(2,:))
    
    ! g_coord save the nodal coordinates in the initial setup
    ! g_coord_aux save the current step nodal coordinates
    mbod(bod)%g_coord_aux = mbod(bod)%g_coord

    ! Setup basic boundary conditions
    mbod(bod)%nf = 1 ! first mark every equations as active
    j=1; k=1
    DO i=1,mbod(bod)%nn
      ! Bottom boundaries
      IF(mbod(bod)%g_coord(2,i)<lowbound+0.01) THEN
        mbod(bod)%nf(2,i) = 0  
      END IF
      ! Right Boundaries (will be tied with the left boundary nodes)
      IF(mbod(bod)%g_coord(1,i)>rightbound-0.01) THEN 
        mbod(bod)%nf(:,i) = 0
      END IF
    END DO
    
    !-Repeat equations of left side of the domain to the right side to 
    ! simulate Tied-Degrees, essentially making a circular 1D column.
    mbod(bod)%tied_nn=0
    j=1;k=1
    DO i=1,mbod(bod)%nn
      IF(mbod(bod)%g_coord(1,i)<leftbound+0.01_iwp)THEN
        mbod(bod)%tied_nn(j,1)=i  ! Mark left nodes
        j=j+1
      ELSE IF(mbod(bod)%g_coord(1,i)>rightbound-0.01_iwp)THEN
        mbod(bod)%tied_nn(k,2)=i ! Mark right nodes
        k=k+1
      ENDIF
    END DO
    
    ! forming nf, and determine number of equations
    CALL formnf(mbod(bod)%nf)
    mbod(bod)%neq=MAXVAL(mbod(bod)%nf)
    
    !- Mark nf of right side nodes to be equal with the left side nodes
    DO i=1,size(mbod(bod)%tied_nn, 1)
      IF(mbod(bod)%tied_nn(i,1)>0)THEN
        mbod(bod)%nf(:,mbod(bod)%tied_nn(i,2))=mbod(bod)%nf(:,mbod(bod)%tied_nn(i,1))
      END IF
    END DO
  END DO Freefield_fixities

  !===========================================================================AS
  ! Variable initialization
  !===========================================================================AS

  ALLOCATE_VARIABLES: DO bod=1,size(mbod)
    !-------------------------------------------------------------------------AS
    ! Allocate single vector variables
    !-------------------------------------------------------------------------AS

    CALL allocate_body(mbod(bod)) ! utilities.f90

    !-------------------------------------------------------------------------AS
    ! Allocate sparse matrices
    !  Data Structures [Check Grifftihs, 2014, Fig 3.20]:
    !  * Fullband Nonsymmetric Matrix : c_matrix, mf_matrix
    !  * Skyline Matrix Lowertriangle : mv, mvis
    !  * CSR/Three Arrays (PARDISO)   : kv (kv_CSR_aux)
    !-------------------------------------------------------------------------AS
    
    ! Determine global streeing factor (g_g) and number of active elements per equation (kdiag)
    mbod(bod)%kdiag=zero
    mbod(bod)%g_g=0
    IF(bod==1)THEN
      DO iel=1,nels
        num=g_num(:,iel)
        CALL num_to_g(num,mbod(bod)%nf,g)
        CALL fkdiag(mbod(bod)%kdiag,g)
        mbod(bod)%g_g(:,iel)=g  
      END DO
    ELSE
      nband=0; bandval=0
      DO iel=1,mbod(bod)%nels
        mbod(bod)%g = 0
        num = mbod(bod)%g_num(:,iel)
        mbod(bod)%g(1:ndof:2) = mbod(bod)%nf(1,num(:)) ! this one in x-dir
        mbod(bod)%g(2:ndof:2) = mbod(bod)%nf(2,num(:)) ! this one in y-dir
        mbod(bod)%g_g(:,iel)  = mbod(bod)%g
        CALL fkdiag(mbod(bod)%kdiag,mbod(bod)%g)
        
        ! bandval is calculated by taking the difference of g vector minmax
        ! only if it is a positive number. (AS) 
        bandval = MAXVAL(mbod(bod)%g,1,mbod(bod)%g>0) -                     &
                  MINVAL(mbod(bod)%g,1,mbod(bod)%g>0)
        IF(nband<bandval) THEN
          nband = bandval
        END IF
      END DO
    END IF
     
    ! Accumulate number of active elements per equation to obtain diagonal index (kdiag)
    DO i=2,mbod(bod)%neq
      mbod(bod)%kdiag(i)=mbod(bod)%kdiag(i)+mbod(bod)%kdiag(i-1)
    END DO
  
    ! Allocate Fullband Nonsystemtric matrices
    ALLOCATE(                                           &
      mbod(bod)%c_matrix(mbod(bod)%neq,2*(nband+1)),    &
      mbod(bod)%mf_matrix(mbod(bod)%neq,2*(nband+1)))

    ! Allocate skyline matrices
    ALLOCATE(                                           &
      mbod(bod)%kv(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%mv(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%cv(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%kp(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%mvis(mbod(bod)%kdiag(mbod(bod)%neq)))
    
    ! Allocate CSR matrices
    ALLOCATE(mbod(bod)%kv_CSR_aux(mbod(bod)%skylength))
  END DO ALLOCATE_VARIABLES


  !*** AS's Thesis
  ! allocate penalty position marker
  ALLOCATE(                               &
    penpos_displacement(0:mbod(1)%neq),   &
    penpos_count(0:mbod(1)%neq),          &
    penalized_stiffness(0:mbod(1)%neq))

  !===========================================================================AS
  ! Material and Simulation Parameters
  !===========================================================================AS
  READ(10,*)dtim,k0,aiters,tol,limit,nsrf

  ! Material parameter valid also valid for Free Field elements [AS]
  DO bod=1,size(mbod)
    mbod(bod)%Young = mbod(bod)%prop(2,1)
    mbod(bod)%Poiss = mbod(bod)%prop(3,1)
    mbod(bod)%mpcp  = mbod(bod)%prop(4,1)
    cpeak = mbod(bod)%prop(4,1)
    mpcr  = mbod(bod)%prop(5,1) 
    phi   = mbod(bod)%prop(6,1) 
  END DO
    
  Maxplastic=0.3
  Maxfric=1.0_iwp
  Minfric=0.50
  Fricmod=(Minfric-Maxfric)/Maxplastic
  SModulus=(mpcr-cpeak)/Maxplastic
  
  fmax=-10.0
  fnew=-10
  stol=0.01_iwp
  evpt=zero
  elso=zero

  ! Create consitutive matrix (dee) for each body
  Constitutive: DO bod=1,size(mbod)
    CALL deemat(mbod(bod)%dee,mbod(bod)%Young,mbod(bod)%Poiss)
  END DO Constitutive

  ! Define gravity field
  Gravf=10.00_iwp
  k0=0.50_iwp
  g_matrix=(/0.0_iwp,-10.0_iwp/)  !--Gravity acting in the vertical direction
  DO bod=1,size(mbod)
    mbod(bod)%g_matrix=(/0.0_iwp,-10.0_iwp/)
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
  
  ! Adjust ground acceleration multiplier
  DO bod=1,size(mbod)
    mbod(bod)%ground_acc=mbod(bod)%ground_acc*3.0_iwp
  END DO

  !===========================================================================AS
  ! Initial Conditions (variables defined here must not be reset, only updated)
  !===========================================================================AS
  
  initial_conditions: DO bod=1,size(mbod)
    ! Stresses
    mbod(bod)%m_stress      = zero
    mbod(bod)%m_stress_ini  = zero
    
    ! Strains and displacements
    mbod(bod)%x1_acum       = zero
    mbod(bod)%ins_acum      = zero
    mbod(bod)%eps_acum      = zero
    mbod(bod)%epsinv_acum   = zero

    ! Acceleration and Velocity
    mbod(bod)%m_acc         = zero
    mbod(bod)%m_velocity    = zero
  END DO initial_conditions
  
  !---------------------------------------------------------------------------AS
  ! Initial stresses (k0 procedure)
  !---------------------------------------------------------------------------AS
  
  HSURF = maxval(mbod(2)%g_coord(2,:))
  DO bod=1,size(mbod)
    Ini_stress_loop:DO i=1,mbod(bod)%nmps
      mbod(bod)%m_stress(2,i) = mbod(bod)%m_dens(i)*Gravf * (mbod(bod)%gm_coord(2,i)-HSURF)
      mbod(bod)%m_stress(1,i) = mbod(bod)%m_stress(2,i)*k0
      mbod(bod)%m_stress(3,i) = zero
      mbod(bod)%m_stress(4,i) = mbod(bod)%m_stress(1,i)
    END DO Ini_stress_loop
    mbod(bod)%m_stress_ini = mbod(bod)%m_stress
  END DO

  !---------------------------------------------------------------------------AS
  ! Initial Configuration (print to paraview)
  !---------------------------------------------------------------------------AS

  DO bod=1,size(mbod)
    IF(bod==1)THEN 
      CALL IO_PARAVIEW(                     &
        input=0,                            &
        node_type=4,                        &
        coord=g_coord,                      &
        num=g_num,                          &
        nf=mbod(bod)%nf,                    &
        kdiag=mbod(bod)%kdiag,              &
        diag=mbod(bod)%diag,                &
        loads=mbod(bod)%loads,              &
        ddylds=mbod(bod)%ddylds,            &
        gravlo=mbod(bod)%gravlo,            &
        vcm=mbod(bod)%vcm,                  &
        cdamp=mbod(bod)%cdamp,              &
        f_earth=mbod(bod)%f_earth,          &
        f_ff=mbod(bod)%f_ff,                &
        x1=mbod(bod)%x1,                    &
        d1x1=mbod(bod)%d1x1,                &
        d2x1=mbod(bod)%d2x1,                &
        directory="Output\Paraview_Cells\", &
        argv="MPM_Body"                     &
      )
      CALL IO_POINT_VIZ(                    &
        input=0,                            &
        coord=mbod(bod)%gm_coord,           &
        a_ins=mbod(bod)%ins_acum,           &
        evpt=mbod(bod)%eps_acum,            &
        m_stress=mbod(bod)%m_stress,        &
        m_stress_inc=mbod(bod)%m_stress,    &
        acc=mbod(bod)%m_acc,                &
        velocity=mbod(bod)%m_velocity,      &
        cohesion=mbod(bod)%mpcp,            &
        devstress=mbod(bod)%Devstress,      &
        meanstress=mbod(bod)%mean_stress,   &
        mpyield=mbod(bod)%mpyield,          &
        directory="Output\Paraview_Point\", &
        argv="MPM_Body"                     &
      )                                     
    ELSE IF(bod==2)THEN
      CALL IO_PARAVIEW(                     &
        input=0,                            &
        node_type=4,                        &
        coord=mbod(bod)%g_coord,            &
        num=mbod(bod)%g_num,                &
        nf=mbod(bod)%nf,                    &
        kdiag=mbod(bod)%kdiag,              &
        diag=mbod(bod)%diag,                &
        loads=mbod(bod)%loads,              &
        ddylds=mbod(bod)%ddylds,            &
        gravlo=mbod(bod)%gravlo,            &
        vcm=mbod(bod)%vcm,                  &
        cdamp=mbod(bod)%cdamp,              &
        f_earth=mbod(bod)%f_earth,          &
        f_ff=mbod(bod)%f_ff,                &
        x1=mbod(bod)%x1,                    &
        d1x1=mbod(bod)%d1x1,                &
        d2x1=mbod(bod)%d2x1,                &
        directory="Output\Paraview_Cells\", &
        argv="FF_Left"                      &
      )
      CALL IO_POINT_VIZ(                    &
        input=0,                            &
        coord=mbod(bod)%gm_coord,           &
        a_ins=mbod(bod)%ins_acum,           &
        evpt=mbod(bod)%eps_acum,            &
        m_stress=mbod(bod)%m_stress,        &
        m_stress_inc=mbod(bod)%m_stress,    &
        acc=mbod(bod)%m_acc,                &
        velocity=mbod(bod)%m_velocity,      &
        cohesion=mbod(bod)%mpcp,            &
        devstress=mbod(bod)%Devstress,      &
        meanstress=mbod(bod)%mean_stress,   &
        mpyield=mbod(bod)%mpyield,          &
        directory="Output\Paraview_Point\", &
        argv="FF_Left"                      &
      )    
    ELSE
      CALL IO_PARAVIEW(                     &
        input=0,                            &
        node_type=4,                        &
        coord=mbod(bod)%g_coord,            &
        num=mbod(bod)%g_num,                &
        nf=mbod(bod)%nf,                    &
        kdiag=mbod(bod)%kdiag,              &
        diag=mbod(bod)%diag,                &
        loads=mbod(bod)%loads,              &
        ddylds=mbod(bod)%ddylds,            &
        gravlo=mbod(bod)%gravlo,            &
        vcm=mbod(bod)%vcm,                  &
        cdamp=mbod(bod)%cdamp,              &
        f_earth=mbod(bod)%f_earth,          &
        f_ff=mbod(bod)%f_ff,                &
        x1=mbod(bod)%x1,                    &
        d1x1=mbod(bod)%d1x1,                &
        d2x1=mbod(bod)%d2x1,                &
        directory="Output\Paraview_Cells\", &
        argv="FF_Right"                     &
      )
      CALL IO_POINT_VIZ(                    &
        input=0,                            &
        coord=mbod(bod)%gm_coord,           &
        a_ins=mbod(bod)%ins_acum,           &
        evpt=mbod(bod)%eps_acum,            &
        m_stress=mbod(bod)%m_stress,        &
        m_stress_inc=mbod(bod)%m_stress,    &
        acc=mbod(bod)%m_acc,                &
        velocity=mbod(bod)%m_velocity,      &
        cohesion=mbod(bod)%mpcp,            &
        devstress=mbod(bod)%Devstress,      &
        meanstress=mbod(bod)%mean_stress,   &
        mpyield=mbod(bod)%mpyield,          &
        directory="Output\Paraview_Point\", &
        argv="FF_Right"                     &
      )    
    END IF
  END DO

  
!===============================================================================
! ----------------------- BEGINNING OF THE TIME STEP ---------------------------
!===============================================================================

  stable=.true.
  step=0
  time_steps: DO w=1,accdata + 2000
  step=step+1 

  !===========================================================================AS
  ! Domain Setups
  !===========================================================================AS
  
  ! Reset variables
  Reset_Variables: DO bod=1,size(mbod)
    ! Nodal Kinematics and Forces
    mbod(bod)%x1         = zero    
    mbod(bod)%kinup_d1x1 = zero
    mbod(bod)%kinup_d2x1 = zero
    mbod(bod)%kinup_Ground_d2x1 = zero
    mbod(bod)%m_mvp      = zero
    mbod(bod)%m_mva      = zero
    mbod(bod)%diag       = zero ! nodal mass
    
    ! Body forces
    mbod(bod)%vcm        = zero
    mbod(bod)%gravlo     = zero
    mbod(bod)%cdamp      = zero
    mbod(bod)%ddylds     = zero
    mbod(bod)%f_ff       = zero
    mbod(bod)%f_earth    = zero
    mbod(bod)%loads      = zero

    ! Additional body loads for FF
    mbod(bod)%ground_loads = zero
    mbod(bod)%c_force    = zero
    mbod(bod)%mf_force   = zero
    
    ! freefield variables
    penpos_count         = 0
    penpos_displacement  = zero
    penalized_stiffness  = zero
    
    ! GIMP variables
    mbod(bod)%valuesg    = zero
    mbod(bod)%GIMP_nodes = zero
    
    ! additional variables (unused)
    mbod(bod)%eps_2      = zero
    mbod(bod)%v_field    = zero
    mbod(bod)%m_field    = zero
    mbod(bod)%a_field    = zero
    mbod(bod)%vel_change = zero
    mbod(bod)%nodecont   = zero
  END DO Reset_Variables
  

  !---------------------------------------------------------------------------AS
  ! Construct Mass Vector (diag) and Calculate Gravity loads (gravlo)
  !---------------------------------------------------------------------------AS
  Body_Solution: DO bod=1,size(mbod)
    IF(bod==1)THEN
      ! The following algorithm are using smethod = 2 [GIMP] or 3 [CMPM]
      DO i=1,mbod(bod)%nmps
        iel=mbod(bod)%a_ele(i)
          
        ! Determine number of support nodes and allocate variable for GIMP shape functions
        CALL GIMP_nodsup(i,g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%mpoints,mbod(bod)%valuesg,mbod(bod)%gimptol,neighb,nny,mbod(bod)%a_ele,mbod(bod)%GIMP_nodes)
        values=mbod(bod)%valuesg(i)

        ! Allocate Temporary Variable
        ALLOCATE(derextend(nodof,values),funextend2(values*2,2),eldddylds(values*2))
        ALLOCATE(beeextend(nst,values*nodof),jac_coord(values,nodof))
        ALLOCATE(mm_gimp(values*2,values*2))
          
        ! Calculate GIMP shape functions and derivatives and get equation index (nf) of the supporting nodes
        eldddylds=zero;funextend2=zero;beeextend=zero;derextend=zero;jac_coord=zero
        CALL GIMP_funder2(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)  
        CALL gimpfunform(i,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values)  
                      
        ! Calculate nodal gravity loads (gravlo)
        mbod(bod)%gravlo(eldddylds) = mbod(bod)%gravlo(eldddylds) + MATMUL(funextend2,mbod(bod)%m_mass(i)*MATMUL(delta,g_matrix))  
        mbod(bod)%gravlo(0) = zero

        ! Form nodal masses (diag; size(diag) = neq)
        mbod(bod)%diag(eldddylds) = mbod(bod)%diag(eldddylds) + MATMUL(funextend2,mbod(bod)%m_mass(i)*delta_1)
        mbod(bod)%diag(0) = zero

        DEALLOCATE(derextend,funextend2,eldddylds,jac_coord,beeextend,mm_gimp)
      END DO
    ELSE !- Free Fields bodies
      i=1
      CALL sample(element,points,weights)
      DO k=1,mbod(bod)%nmps
        mbod(bod)%diag=zero
        iel=mbod(bod)%a_ele(k) 
        
        ! build element consistent mass matrix (ecm) with ecmat2
        num=mbod(bod)%g_num(:,iel)
        coord=TRANSPOSE(mbod(bod)%g_coord(:,num)) 
        CALL shape_der(der,mbod(bod)%mpoints,k)
        jac=MATMUL(der,coord)
        jac(2,1)=zero; jac(1,2)=zero ! ensure the off diagonal is zero
        det=determinant(jac)
        CALL shape_fun(fun,mbod(bod)%mpoints,k)
        CALL ecmat2(ecm,fun,ndof,nodof)

        ! element mass lumping
        DO m=1,size(ecm, 1)
          ecm_acum=zero
          DO j=1,size(ecm, 2)
            ecm_acum(m,m)=ecm_acum(m,m)+ecm(m,j)
          END DO
          ecm(m,:)=zero
          ecm(m,m)=ecm_acum(m,m)
        END DO
        mm_s=mm_s+ecm*mbod(bod)%m_dens(k)*det*weights(i)
      
        ! increment i to cycle weigths for next stress points
        ! Note: this will break if a_ele is not ordered correctly
        i=i+1; IF(i>nip)THEN
          ! The element steering vector is set such that in x-direction, g is 0
          ! because the gravity load and body mass only go downwards
          mbod(bod)%g=mbod(bod)%g_g(:,iel) 
          mbod(bod)%g(1:size(mbod(bod)%g):2)=0
          
          ! form mass vector (diag)
          CALL formlump(mbod(bod)%diag,mm_s,mbod(bod)%g)
          mbod(bod)%gravlo(mbod(bod)%g) = mbod(bod)%gravlo(mbod(bod)%g) + mbod(bod)%diag(mbod(bod)%g)*(-Gravf)  
  
          mbod(bod)%diag(0)=zero
          mm_s=zero
          i=1
        END IF
      END DO
    END IF
  END DO Body_Solution
   

  !---------------------------------------------------------------------------AS
  ! Construct Decomposed Modified Stiffness Matrix (kp)
  !---------------------------------------------------------------------------AS

  ! Compute domain stiffnes matrices (kv), domain mass matrices (mv), and their modified stiffness (kp)
  ! This matrices are computed using only the double mapping technique and GIMP, regular MPM is not used
  Stiffness:DO bod=1,size(mbod)
    FEM_MOD_STIFF: IF(bod>1)THEN
      mm_s         = zero
      mm_mf_acum   = zero
      mm_acum      = zero      
      mbod(bod)%km = zero 
      
      mbod(bod)%kp = zero 
      mbod(bod)%kv = zero
      mbod(bod)%cv = zero
      mbod(bod)%mv = zero

      mbod(bod)%KGC     = zero
      mbod(bod)%CCS     = zero
      mbod(bod)%MMS     = zero
      mbod(bod)%MOD_MTX = zero

      nnzacum              = 0
      mbod(bod)%ia_aux     = 0
      mbod(bod)%ja_aux     = 0
      mbod(bod)%skylength  = 0
      mbod(bod)%kv_CSR_aux = zero
      
      i=1
      CALL sample(element,points,weights)
        KM_MV_2D:DO k=1,mbod(bod)%nmps
          iel=mbod(bod)%a_ele(k)
          num=mbod(bod)%g_num(:,iel)
          coord=TRANSPOSE(mbod(bod)%g_coord(:,num))
          
          mbod(bod)%g=0
          mbod(bod)%g=mbod(bod)%g_g(:,iel)

          !-------------------------------------------------------------------AS
          ! Build stiffness matrix (kv)
          !-------------------------------------------------------------------AS
        
          CALL shape_der(der,mbod(bod)%mpoints,k)
          jac=MATMUL(der,coord)
          jac(2,1)=zero; jac(1,2)=zero ! ensure the off diagonal is zero
          det=determinant(jac)
          CALL invert(jac)
          deriv=MATMUL(jac,der)
          CALL beemat(bee,deriv)
              
          mbod(bod)%km = mbod(bod)%km + MATMUL(MATMUL(TRANSPOSE(bee),mbod(bod)%dee),bee)*det*weights(i)

          !-------------------------------------------------------------------AS
          ! Build diagonal mass matrix (mv)
          !-------------------------------------------------------------------AS

          CALL shape_fun(fun,mbod(bod)%mpoints,k)
          CALL ecmat2(ecm,fun,ndof,nodof)
          DO m=1,size(ecm,1)
            ecm_acum=zero
            DO j=1,size(ecm,2)
              ecm_acum(m,m)=ecm_acum(m,m)+(ecm(m,j))
            END DO
            ecm(m,:)=zero
            ecm(m,m)=ecm_acum(m,m)
          END DO
          mm_s=mm_s+ecm*mbod(bod)%m_dens(k)*det*weights(i)

          !-------------------------------------------------------------------AS
          ! Build modified stiffness matrix (kp / MOD_MTX)
          !-------------------------------------------------------------------AS
          ! at the end of each element form the modified stiffness matrix (kp or MOD_MTK)
          ! according to the incremental FEM iterations (MOD_MTK = 4MMS/dtim**2 + 2C/dtim + KGC) ; C = fk*KGC + fm*MMS 
          IF(i==nip)THEN
            IF(use_damping)THEN
              fm=0.4936788455641104_iwp; fk=0.0013158595278225688_iwp
            ELSE
              fm=0.0_iwp; fk=0.0_iwp
            END IF

            mbod(bod)%KGC = zero
            mbod(bod)%MMS = zero
            mbod(bod)%MMS = mm_s
            mbod(bod)%KGC = mbod(bod)%km
            mbod(bod)%CCS = mbod(bod)%KGC*fk+mbod(bod)%MMS*fm
            mbod(bod)%MOD_MTX = 4.0_iwp*mbod(bod)%MMS/dtim**2.0_iwp +  &
                                2.0_iwp/dtim*mbod(bod)%CCS +           &
                                mbod(bod)%KGC
                                
            CALL fsparv(mbod(bod)%cv,mbod(bod)%CCS,mbod(bod)%g,mbod(bod)%kdiag)
            CALL fsparv(mbod(bod)%mv,mbod(bod)%MMS,mbod(bod)%g,mbod(bod)%kdiag)
            CALL formspars_unsym(   &
              mbod(bod)%ndof,       &
              mbod(bod)%g,          &
              mbod(bod)%MOD_MTX,    &
              mbod(bod)%ia_aux,     &
              mbod(bod)%ja_aux,     &
              mbod(bod)%kv_CSR_aux, &
              mbod(bod)%skylength )

            mm_s              = zero
            mm_mf_acum        = zero
            mm_acum           = zero
            mbod(bod)%km      = zero 
            mbod(bod)%MOD_MTX = zero
          END IF    
        
        i=i+1; IF(i>nip)i=1
      END DO KM_MV_2D
    END IF FEM_MOD_STIFF


    MPM_MOD_STIFF: IF(bod==1)THEN
      mm_s         = zero
      mm_mf_acum   = zero
      mm_acum      = zero      
      mbod(bod)%km = zero 
      
      mbod(bod)%kp = zero 
      mbod(bod)%kv = zero
      mbod(bod)%cv = zero
      mbod(bod)%mv = zero

      KM_MV:DO k=1,mbod(bod)%nmps
      
        !---------------------------------------------------------------------AS
        !--- Build stiffness matrix (kv)
        !---------------------------------------------------------------------AS
        ! Stiffness matrix is constructed using non-constant stiffness (kconst /= 0)
        ! This loop goes until 4 since 4 is the maximum number of elements affected by a material point
        El_stiff:DO i=1,4 
          iel=mbod(bod)%elemmpoints(k,i)
          Full_el:IF(mbod(bod)%elemmpoints(k,i)>0)THEN !--If is true, an element is affected by a material point (material point can be outside the element)
            num=g_num(:,iel)
            nod_num(:,1)=num
            coord=TRANSPOSE(g_coord(:,num))
            g=mbod(bod)%g_g(:,iel)
            km_gauss=zero   
            waverage=4.0/nip

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
                
              km_gauss = km_gauss + MATMUL(MATMUL(TRANSPOSE(bee_gauss),dee_scal),bee_gauss)*det*weights(s)
            END DO
              
            CALL fsparv(mbod(bod)%kv,km_gauss,g,mbod(bod)%kdiag) ! stiffness matrix (KV)
          END IF Full_el
        END DO El_stiff
        
        !---------------------------------------------------------------------AS
        !--- Build diagonal mass matrix (mv)
        !---------------------------------------------------------------------AS

        values=mbod(bod)%valuesg(k) ! values have been calculated in the previous section
        ALLOCATE(derextend(nodof,values),funextend2(values*2,2),eldddylds(values*2),mm_gimp(values*2,values*2))
          
        CALL GIMP_funder2(k,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)
        CALL gimpfunform(k,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values)  
        
        mm_gimp=zero
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

      !-----------------------------------------------------------------------AS
      ! Construct modified stiffness matrix (kp) and Rayleigh damping (cv)
      !-----------------------------------------------------------------------AS
      IF(use_damping)THEN
        fm=0.4936788455641104_iwp; fk=0.0013158595278225688_iwp
      ELSE
        fm=0.0_iwp; fk=0.0_iwp
      END IF
      mbod(bod)%cv = fm*mbod(bod)%mv + fk*mbod(bod)%kv  ! cv is already in skyline form
      mbod(bod)%kp = 4.0_iwp*mbod(bod)%mv/dtim**2.0_iwp +  &
                     2.0_iwp/dtim*(mbod(bod)%cv) + mbod(bod)%kv
      
      !-----------------------------------------------------------------------AS
      ! Apply penalty stiffness to modified stiffness matrix
      !-----------------------------------------------------------------------AS
  
      !*** AS's Thesis
      !--- Determine freefield boundaries displacement
      CALL get_ff_displacement(         &
        mpm_disp=penpos_displacement,   &
        mpm_counter=penpos_count,       &
        iel_boundary=iel_boundary,      &
        mpm_g_num=g_num,                &
        mpm_nf=mbod(1)%nf,              &
        ff_disp=mbod(2)%x1,             &
        ff_g_num=mbod(2)%g_num,         &
        ff_nf=mbod(2)%nf                &
      )
        
      !--- Apply penalty to modified stiffness matrix (kp) and form boundary force
      DO i=1,mbod(1)%neq
        IF(penpos_count(i)>0)THEN
          mbod(1)%kp(mbod(1)%kdiag(i)) = mbod(1)%kp(mbod(1)%kdiag(i)) + penalty
          penalized_stiffness(i) = penalized_stiffness(i) + mbod(1)%kp(mbod(1)%kdiag(i))
        END IF
      END DO

      !-----------------------------------------------------------------------AS
      ! Cholesky Decomposition
      !-----------------------------------------------------------------------AS
      CALL sparin(mbod(bod)%kp,mbod(bod)%kdiag)
    END IF MPM_MOD_STIFF

  END DO Stiffness
 
  
  !---------------------------------------------------------------------------AS
  ! Pardiso CSR format and LU Decomposition
  !---------------------------------------------------------------------------AS
  
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
    CALL pardiso(           &
      mbod(bod)%pt,         &
      maxfct,               &
      mnum,                 &
      mtype,                &
      phase,                &
      mbod(bod)%neq,        &
      mbod(bod)%kv_CSR,     & 
      mbod(bod)%ia,         &
      mbod(bod)%ja,         &
      idum,                 &
      nrhs,                 &
      mbod(bod)%iparm,      &
      msglvl,               &
      ddum,                 &
      ddum,                 &
      error,                &
      dparm                 &
    )

    ! A*x = b
    ! pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
    phase = 22  !// LU decompose
    CALL pardiso(           &
      mbod(bod)%pt,         &
      maxfct,               &
      mnum,                 &
      mtype,                &
      phase,                &
      mbod(bod)%neq,        &
      mbod(bod)%kv_CSR,     &
      mbod(bod)%ia,         &
      mbod(bod)%ja,         &
      idum,                 &
      nrhs,                 &
      mbod(bod)%iparm,      &
      msglvl,               &
      ddum,                 &
      ddum,                 &
      error,                &
      dparm                 &
    )
  END DO

  !===========================================================================AS
  ! Mapping from Material Points to Nodes (MPM only)
  !===========================================================================AS

  ! calculate initial nodal acceleration (d2x1) and velocity (d1x1)
  DO bod=1,1
    mbod(bod)%d1x1     = zero
    mbod(bod)%d2x1     = zero
    MP_VelAccP: DO i=1,mbod(bod)%nmps
      iel=mbod(bod)%a_ele(i)    
      CALL GIMP_nodsup(i,g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%mpoints, &
        mbod(bod)%valuesg,mbod(bod)%gimptol,neighb,nny,mbod(bod)%a_ele,mbod(bod)%GIMP_nodes)
        
      values=mbod(bod)%valuesg(i)
      ALLOCATE(derextend(nodof,values),jac_coord(values,nodof),funextend2(values*2,2),	&
        eldddylds(values*2),beeextend(nst,values*nodof))

      eldddylds=zero;funextend2=zero;beeextend=zero;derextend=zero;jac_coord=zero
  
      CALL GIMP_funder2(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2)  !--creates the matrix with the values of shape functions and derivatives
      CALL gimpfunform(i,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values)  !--collect in eldddylds the nuber of the equations in the support domain
                              
      mbod(bod)%m_mvp(eldddylds)=mbod(bod)%m_mvp(eldddylds)+MATMUL(MATMUL(funextend2,mbod(bod)%m_mass(i)*delta),mbod(bod)%m_velocity(:,i))
      mbod(bod)%m_mva(eldddylds)=mbod(bod)%m_mva(eldddylds)+MATMUL(MATMUL(funextend2,mbod(bod)%m_mass(i)*delta),mbod(bod)%accb(:,i))

      mbod(bod)%m_mvp(0)=zero;mbod(bod)%m_mva(0)=zero
      DEALLOCATE(derextend,funextend2,eldddylds,jac_coord,beeextend)
    END DO MP_VelAccP

    DO i=1,mbod(bod)%neq
      mbod(bod)%d1x1(i)=mbod(bod)%m_mvp(i)/mbod(bod)%diag(i) ! m_mvp & diag is the contribution from all bodies and at the contact d1x1 could be also zero
      mbod(bod)%d2x1(i)=mbod(bod)%m_mva(i)/mbod(bod)%diag(i)
    END DO
  END DO 
  

  !============================================================================AS
  ! Iteration loops 
  !============================================================================AS
  
  !--- Solve freefield bodies stresses and displacement levels
  iters=0; newrapconv1=.false.

  COMBINED_NR: DO WHILE(iters<limit)
    iters = iters + 1
    !-------------------------------------------------------------------------AS
    ! Determine Nodal Kinematics
    !-------------------------------------------------------------------------AS

    ! calculate nodal acceleration (kinup_d2x1) and velocity (kinup_d1x1) increment
    DO bod=1,size(mbod)
      mbod(bod)%kinup_d2x1 = (4.0_iwp*mbod(bod)%x1/dtim**2.0_iwp) - &
                             (4.0_iwp*mbod(bod)%d1x1/dtim) -        &
                             mbod(bod)%d2x1         
      mbod(bod)%kinup_d1x1 = 2.0_iwp*mbod(bod)%x1/dtim - mbod(bod)%d1x1
      mbod(bod)%kinup_d2x1(0) = zero
      mbod(bod)%kinup_d1x1(0) = zero
    END DO 

    ! calculate nodal acceleration considering ground motion (kinup_Ground_d2x1)
    DO bod=2,size(mbod)
      lowbound=minval(mbod(bod)%g_coord(2,:)) 
      IF(w<=accdata.and.w>=1)THEN
        DO i=1,(mbod(bod)%ney+1)*(mbod(bod)%nex)
          IF(mbod(bod)%g_coord(2,i)<lowbound+0.01_iwp)THEN
            mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i)) = mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i)) - &
                                                             mbod(bod)%kinup_d2x1(mbod(bod)%nf(1,i)) +        &    
                                                             mbod(bod)%ground_acc(w)
            mbod(bod)%kinup_Ground_d2x1(0)=zero
          END IF    
        END DO
      END IF
    END DO

    !-------------------------------------------------------------------------AS
    ! Determine Nodal Forces
    !-------------------------------------------------------------------------AS

    ! Calculate internal force (ddylds) from particle stresses (m_stress)
    DO bod=1,1 ! MPM Body
      mbod(bod)%ddylds=zero
      DO i=1,mbod(bod)%nmps
        iel=mbod(bod)%a_ele(i)   
        num=g_num(:,iel)
        coord=TRANSPOSE(g_coord(:,num))
        g=mbod(bod)%g_g(:,iel)
        
        values=mbod(bod)%valuesg(i)
        ALLOCATE(derextend(nodof,values),funextend2(values*2,2),              &
          jac_coord(values,nodof),eldddylds(values*2),beeextend(nst,values*nodof))
        
        eldddylds=zero;funextend2=zero;beeextend=zero;derextend=zero;jac_coord=zero
        CALL GIMP_funder2(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp, &
          mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2) 
        CALL gimpfunform(i,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values)
        CALL Sup_coord(i,mbod(bod)%GIMP_nodes,g_coord,jac_coord)
        jac=MATMUL(derextend,jac_coord)
        jac(2,1)=0.0;jac(1,2)=0.0
        CALL invert(jac)
        derextend=MATMUL(jac,derextend)
        CALL beematgimp(i,beeextend,derextend,values)

        ! Calculate Body Loads (ddylds)
        sigma=mbod(bod)%m_stress(:,i) 
        mbod(bod)%ddylds(eldddylds)=mbod(bod)%ddylds(eldddylds) + MATMUL(sigma,beeextend)*(4.0*mbod(bod)%lp_mp(1,i)*mbod(bod)%lp_mp(2,i))
        mbod(bod)%ddylds(0)=zero 

        DEALLOCATE(derextend,funextend2,eldddylds,jac_coord,beeextend)
      END DO
    END DO 

    DO bod=2,size(mbod) ! Freefield Bodies
      mbod(bod)%ddylds=zero
      a=1
      CALL sample2(element,points,weights)
      DO i=1,mbod(bod)%nmps
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
        
        ! Calculate Body Loads (ddylds)
        sigma=mbod(bod)%m_stress(:,i) 
        mbod(bod)%ddylds(mbod(bod)%g) = mbod(bod)%ddylds(mbod(bod)%g) + MATMUL(TRANSPOSE(bee),sigma)*det*weights(a)
        mbod(bod)%ddylds(0)=zero

        a=a+1; IF(a>nip)a=1
      END DO
    END DO 


    ! calculate damping force (cdamp; F=cv)
    Damping_force: DO bod=1,size(mbod)
      mbod(bod)%cdamp=zero
      IF(mbod(bod)%nmps>1)THEN
        CALL linmul_sky(mbod(bod)%cv,mbod(bod)%kinup_d1x1,mbod(bod)%cdamp,mbod(bod)%kdiag)
        mbod(bod)%cdamp(0)=zero
      END IF
    END DO Damping_force


    ! calculate inertial force (vcm; F=ma)
    Force_MA: DO bod=1,size(mbod)
      mbod(bod)%vcm=zero
      IF(mbod(bod)%nmps>1)THEN 
        CALL linmul_sky(mbod(bod)%mv,mbod(bod)%kinup_d2x1,mbod(bod)%vcm,mbod(bod)%kdiag)  !--Multiplication of the mass matrix per the acceleration vector vcm=Ma 
        mbod(bod)%vcm(0)=zero   
      END IF
    END DO Force_MA 


    ! calculate ground force (f_earth; F=ma)
    Ground_F: DO bod=2,size(mbod)
      mbod(bod)%f_earth = zero
      IF(mbod(bod)%nmps>1)THEN 
        CALL linmul_sky(mbod(bod)%mv,mbod(bod)%kinup_Ground_d2x1,mbod(bod)%f_earth,mbod(bod)%kdiag)  
        mbod(bod)%f_earth(0)=zero
      END IF
    END DO Ground_F 


    !-------------------------------------------------------------------------AS
    ! Solve system of equations
    !-------------------------------------------------------------------------AS

    ! Solve Freefields system of equations, obtain displacements
    FF_Displacements: DO bod=2,size(mbod)
      mbod(bod)%loads = mbod(bod)%gravlo - mbod(bod)%ddylds -    &
                        mbod(bod)%vcm - mbod(bod)%cdamp +        &
                        mbod(bod)%f_earth                 
      mbod(bod)%loads(0)=zero

      ! Check whether the calculated loads has converge (saved in newrapconv1)
      IF(iters==1) mbod(bod)%loads_ini = mbod(bod)%loads

      ! Solve the system of equations (Ma+Cv+Ku = F)
      ! pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
      ! A*X = B ; X = mbod(bod)%loads(1:mbod(bod)%neq) ; B = mbod(bod)%residual(1:mbod(bod)%neq)
      mbod(bod)%residual = zero
      mbod(bod)%iparm(8) = 10 ! max numbers of iterative refinement steps
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
      mbod(bod)%loads(0)=zero; mbod(bod)%residual(0)=zero 

      mbod(bod)%loads       = mbod(bod)%residual  ! assign displacement to loads (weird indeed)
      mbod(bod)%x1          = mbod(bod)%x1+mbod(bod)%loads ! add displacement increments
      mbod(bod)%x1(0)       = zero
    END DO FF_Displacements

    
    !*** AS's Thesis
    ! Determine freefield boundaries displacement in the MPM domain
    CALL get_ff_displacement(         &
      mpm_disp=penpos_displacement,   &
      mpm_counter=penpos_count,       &
      iel_boundary=iel_boundary,      &
      mpm_g_num=g_num,                &
      mpm_nf=mbod(1)%nf,              &
      ff_disp=mbod(2)%residual,       &
      ff_g_num=mbod(2)%g_num,         &
      ff_nf=mbod(2)%nf                &
    )
      

    ! Solve MPM system of equations, obtain displacements
    MPM_DISPLACEMENTS: DO bod=1,1
      mbod(bod)%loads=zero
      mbod(bod)%loads=mbod(bod)%gravlo - mbod(bod)%ddylds - &
                      mbod(bod)%vcm - mbod(bod)%cdamp
      
      !*** AS's Thesis
      ! Calculate boundary penalized force
      DO i=1,mbod(bod)%neq
        IF(penpos_count(i) > 0)THEN
          mbod(bod)%loads(i) = penpos_displacement(i) * penalized_stiffness(i)
        END IF
      END DO
      
      mbod(bod)%loads(0)=zero 
      
      CALL spabac(mbod(bod)%kp,mbod(bod)%loads,mbod(bod)%kdiag)  
      mbod(bod)%loads(0)=zero 
      mbod(bod)%x1=mbod(bod)%x1+mbod(bod)%loads
      mbod(bod)%x1(0)=zero
    END DO MPM_DISPLACEMENTS  


    !-------------------------------------------------------------------------AS
    ! Update Nodal Kinematics
    !-------------------------------------------------------------------------AS

    ! Calculate new acceleration and velocity increments (kinup_d2x1 & kinup_d1x1)
    ! based on the new displacement increment (x1) with Newmark's Equations
    DO bod=1,size(mbod)
      mbod(bod)%kinup_d2x1=(4.0_iwp*mbod(bod)%x1/dtim**2)-(4.0_iwp*mbod(bod)%d1x1/dtim)-mbod(bod)%d2x1 
      mbod(bod)%kinup_d1x1=2.0_iwp*mbod(bod)%x1/dtim-mbod(bod)%d1x1  
      mbod(bod)%kinup_d1x1(0)=zero;mbod(bod)%kinup_d2x1(0)=zero
    END DO 

    !-------------------------------------------------------------------------AS
    ! Stress recovery
    !-------------------------------------------------------------------------AS
    DO bod=2,size(mbod)
      k=1
      DO i=1,mbod(bod)%nmps
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
        
        k=k+1; IF(k>nip) k=1
        
        ! calculate updated strains distribution, remember loads is displacement
        eld   = mbod(bod)%loads(g_s)
        eps_s = MATMUL(bee,eld)

        ! calculate the stress increase based on the increase strain
        mbod(bod)%eps_acum(:,i) = mbod(bod)%eps_acum(:,i) + eps_s
        sigma = MATMUL(mbod(bod)%dee,eps_s)

        ! Accumulate stresses
        mbod(bod)%m_stress(:,i) = mbod(bod)%m_stress(:,i) + sigma

        ! calculate the new stresses invariants
        CALL invar(mbod(bod)%m_stress(:,i),sigm,dsbar,lode_theta)
        ps1 = sigm+(2.0/3.0)*dsbar*sin(lode_theta-(2.0*3.1415926/3.0))
        ps2 = sigm+(2.0/3.0)*dsbar*sin(lode_theta)
        ps3 = sigm+(2.0/3.0)*dsbar*sin(lode_theta+(2.0*3.1415926/3.0))
        mbod(bod)%Devstress(i)   = (1.0/sqrt(two))*sqrt((ps1-ps2)**2+&
                                   (ps2-ps3)**2+(ps3-ps1)**2)
        mbod(bod)%mean_stress(i) = (mbod(bod)%m_stress(1,i) +  &
                                    mbod(bod)%m_stress(2,i) +  &
                                    mbod(bod)%m_stress(4,i))/3.0_iwp
      END DO
    END DO

    DO bod=1,1
      DO i=1,mbod(bod)%nmps
        iel=mbod(bod)%a_ele(i)
        num=g_num(:,iel)
        coord=TRANSPOSE(g_coord(:,num))
        g=mbod(bod)%g_g(:,iel)
    
        !- Compute stresses using neighbouring element strains (eps) with CMPM
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

            IF(m==1) sigma_1=MATMUL(mbod(bod)%dee,eps)
          END DO average
          IF(aver==1) eps=eps
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
        
        ! Stress updated considering the new "eps_acum" level
        mbod(bod)%eps_acum(:,i) = mbod(bod)%eps_acum(:,i) + eps
        sigma = MATMUL(mbod(bod)%dee,eps)  

        ! Accumulate stresses
        mbod(bod)%m_stress(:,i) = mbod(bod)%m_stress(:,i) + sigma
      
        ! calculate the new stresses invariants
        CALL invar(mbod(bod)%m_stress(:,i),sigm,dsbar,lode_theta)
        ps1 = sigm+(2.0/3.0)*dsbar*sin(lode_theta-(2.0*3.1415926/3.0))
        ps2 = sigm+(2.0/3.0)*dsbar*sin(lode_theta)
        ps3 = sigm+(2.0/3.0)*dsbar*sin(lode_theta+(2.0*3.1415926/3.0))
        mbod(bod)%Devstress(i)   = (1.0/sqrt(two))*sqrt((ps1-ps2)**2+&
                                   (ps2-ps3)**2+(ps3-ps1)**2)
        mbod(bod)%mean_stress(i) = (mbod(bod)%m_stress(1,i) +  &
                                    mbod(bod)%m_stress(2,i) +  &
                                    mbod(bod)%m_stress(4,i))/3.0_iwp
      END DO
    END DO   
  END DO COMBINED_NR
 
  
  !=============================================================================AS
  ! Calculate new nodal acceleration (d2x1) and velocity (d1x1)
  !=============================================================================AS
  
  ! x1_acum is here just for post-simulation purpose. Not for actual calculation
  ! With the FEM, nodal kinematics is kept on the nodes.
  DO bod=2,size(mbod)
    mbod(bod)%x1_acum = mbod(bod)%x1_acum+mbod(bod)%x1
    mbod(bod)%d1x1    = mbod(bod)%kinup_d1x1
    mbod(bod)%d2x1    = mbod(bod)%kinup_d2x1
  END DO
  
  ! Update the node displacement for visualization
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
  ! For FEM, this step is solely for visualization purposes
  ! Eitherway, the Stress/Material Points (gm_coord) coordinates are updated here
  
  REVERSE_MAPPING: DO bod=1,size(mbod)
    IF(bod==1)THEN
      mbod(bod)%vccp=zero
      DO i=1,mbod(bod)%nmps
        values=mbod(bod)%valuesg(i)
        ALLOCATE(derextend(nodof,values),beeextend(nst,values*nodof),funextend2(values*2,2),	&
          eldddylds(values*2),funextend(values),jac_coord(values,nodof),eldCMPM(values*2))
        CALL GIMP_funder(i,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend) !--creates the matrix with the values of shape functions and derivatives
        
        CALL eldformgimp(i,eldCMPM,mbod(bod)%kinup_d2x1,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g)
        mbod(bod)%accp(1,i)=DOT_PRODUCT(eldCMPM(1:ndof:2),funextend)
        mbod(bod)%accp(2,i)=DOT_PRODUCT(eldCMPM(2:ndof:2),funextend)
        mbod(bod)%m_acc(:,i)=mbod(bod)%accp(:,i)
        mbod(bod)%m_velocity(:,i)=mbod(bod)%m_velocity(:,i) + 0.5_iwp*(mbod(bod)%accp(:,i) + mbod(bod)%accb(:,i))*dtim

        CALL eldformgimp(i,eldCMPM,mbod(bod)%x1,mbod(bod)%nf,mbod(bod)%GIMP_nodes,values,mbod(bod)%g_g)
        mbod(bod)%ins(1,i)=DOT_PRODUCT(eldCMPM(1:ndof:2),funextend)
        mbod(bod)%ins_acum(:,i)=mbod(bod)%ins_acum(:,i)+mbod(bod)%ins(:,i)
        mbod(bod)%gm_coord(:,i)=mbod(bod)%gm_coord(:,i)+mbod(bod)%ins(:,i)
        
        DEALLOCATE(derextend,beeextend,funextend2,eldddylds,funextend,jac_coord,eldCMPM)
      END DO
      mbod(bod)%accb=mbod(bod)%accp
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
        mbod(bod)%m_velocity(:,i) = mbod(bod)%m_velocity(:,i) + 0.5_iwp*(mbod(bod)%accp(:,i) + mbod(bod)%accb(:,i))*dtim
        
        !** Convect material points
        eld=mbod(bod)%x1(mbod(bod)%g)
        mbod(bod)%ins(1,i)   = DOT_PRODUCT(eld(1:ndof:2),fun)  
        mbod(bod)%ins(2,i)   = DOT_PRODUCT(eld(2:ndof:2),fun) 
        mbod(bod)%ins_acum(:,i) = mbod(bod)%ins_acum(:,i) + mbod(bod)%ins(:,i)
        mbod(bod)%gm_coord(:,i) = mbod(bod)%gm_coord(:,i) + mbod(bod)%ins(:,i)
      END DO
      mbod(bod)%accb=mbod(bod)%accp
    END IF
  END DO REVERSE_MAPPING
  


  !=============================================================================AS
  ! Numerical Cleanup
  !=============================================================================AS
  ! Clean up memory for Pardiso solver
  phase = -1
  DO bod=2,size(mbod)
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
 
  write(800, '(*(E15.5 ","))') mbod(1)%ins_acum
  write(810, '(*(E15.5 ","))') mbod(1)%m_velocity
  write(820, '(*(E15.5 ","))') mbod(1)%m_acc
  write(830, '(*(E15.5 ","))') mbod(2)%ins_acum
  write(840, '(*(E15.5 ","))') mbod(2)%m_velocity
  write(850, '(*(E15.5 ","))') mbod(2)%m_acc
  
  ! -- Loop to save data from both bodies and print it in point_vis
  IF(step/printval*printval==step)THEN
    PRINT '("Steps :" (I10) "/" (I10))', step, accdata
    DO bod=1,size(mbod)
      IF(bod==1)THEN 
        CALL IO_PARAVIEW(                     &
          input=step,                         &
          node_type=4,                        &
          coord=g_coord,                      &
          num=g_num,                          &
          nf=mbod(bod)%nf,                    &
          kdiag=mbod(bod)%kdiag,              &
          diag=mbod(bod)%diag,                &
          loads=mbod(bod)%loads,              &
          ddylds=mbod(bod)%ddylds,            &
          gravlo=mbod(bod)%gravlo,            &
          vcm=mbod(bod)%vcm,                  &
          cdamp=mbod(bod)%cdamp,              &
          f_earth=mbod(bod)%f_earth,          &
          f_ff=mbod(bod)%f_ff,                &
          x1=mbod(bod)%x1,                    &
          d1x1=mbod(bod)%d1x1,                &
          d2x1=mbod(bod)%d2x1,                &
          directory="Output\Paraview_Cells\", &
          argv="MPM_Body"                     &
        )
        CALL IO_POINT_VIZ(                    &
          input=step,                         &
          coord=mbod(bod)%gm_coord,           &
          a_ins=mbod(bod)%ins_acum,           &
          evpt=mbod(bod)%eps_acum,            &
          m_stress=mbod(bod)%m_stress,        &
          m_stress_inc=mbod(bod)%m_stress,    &
          acc=mbod(bod)%m_acc,                &
          velocity=mbod(bod)%m_velocity,      &
          cohesion=mbod(bod)%mpcp,            &
          devstress=mbod(bod)%Devstress,      &
          meanstress=mbod(bod)%mean_stress,   &
          mpyield=mbod(bod)%mpyield,          &
          directory="Output\Paraview_Point\", &
          argv="MPM_Body"                     &
        )                                     
      ELSE IF(bod==2)THEN
        CALL IO_PARAVIEW(                     &
          input=step,                         &
          node_type=4,                        &
          coord=mbod(bod)%g_coord_aux,        &
          num=mbod(bod)%g_num,                &
          nf=mbod(bod)%nf,                    &
          kdiag=mbod(bod)%kdiag,              &
          diag=mbod(bod)%diag,                &
          loads=mbod(bod)%loads,              &
          ddylds=mbod(bod)%ddylds,            &
          gravlo=mbod(bod)%gravlo,            &
          vcm=mbod(bod)%vcm,                  &
          cdamp=mbod(bod)%cdamp,              &
          f_earth=mbod(bod)%f_earth,          &
          f_ff=mbod(bod)%f_ff,                &
          x1=mbod(bod)%x1,                    &
          d1x1=mbod(bod)%d1x1,                &
          d2x1=mbod(bod)%d2x1,                &
          directory="Output\Paraview_Cells\", &
          argv="FF_Left"                      &
        )
        CALL IO_POINT_VIZ(                    &
          input=step,                         &
          coord=mbod(bod)%gm_coord,           &
          a_ins=mbod(bod)%ins_acum,           &
          evpt=mbod(bod)%eps_acum,            &
          m_stress=mbod(bod)%m_stress,        &
          m_stress_inc=mbod(bod)%m_stress,    &
          acc=mbod(bod)%m_acc,                &
          velocity=mbod(bod)%m_velocity,      &
          cohesion=mbod(bod)%mpcp,            &
          devstress=mbod(bod)%Devstress,      &
          meanstress=mbod(bod)%mean_stress,   &
          mpyield=mbod(bod)%mpyield,          &
          directory="Output\Paraview_Point\",  &
          argv="FF_Left"                      &
        )    
      ELSE
        CALL IO_PARAVIEW(                     &
          input=step,                         &
          node_type=4,                        &
          coord=mbod(bod)%g_coord_aux,        &
          num=mbod(bod)%g_num,                &
          nf=mbod(bod)%nf,                    &
          kdiag=mbod(bod)%kdiag,              &
          diag=mbod(bod)%diag,                &
          loads=mbod(bod)%loads,              &
          ddylds=mbod(bod)%ddylds,            &
          gravlo=mbod(bod)%gravlo,            &
          vcm=mbod(bod)%vcm,                  &
          cdamp=mbod(bod)%cdamp,              &
          f_earth=mbod(bod)%f_earth,          &
          f_ff=mbod(bod)%f_ff,                &
          x1=mbod(bod)%x1,                    &
          d1x1=mbod(bod)%d1x1,                &
          d2x1=mbod(bod)%d2x1,                &
          directory="Output\Paraview_Cells\", &
          argv="FF_Right"                     &
        )
        CALL IO_POINT_VIZ(                    &
          input=step,                         &
          coord=mbod(bod)%gm_coord,           &
          a_ins=mbod(bod)%ins_acum,           &
          evpt=mbod(bod)%eps_acum,            &
          m_stress=mbod(bod)%m_stress,        &
          m_stress_inc=mbod(bod)%m_stress,    &
          acc=mbod(bod)%m_acc,                &
          velocity=mbod(bod)%m_velocity,      &
          cohesion=mbod(bod)%mpcp,            &
          devstress=mbod(bod)%Devstress,      &
          meanstress=mbod(bod)%mean_stress,   &
          mpyield=mbod(bod)%mpyield,          &
          directory="Output\Paraview_Point\", &
          argv="FF_Right"                     &
        )    
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
      IF(rowy - mbod(bod)%dist_y > 0)THEN
        IF(right_boundary(rowy-mbod(bod)%dist_y) < colx) right_boundary(rowy-mbod(bod)%dist_y) = colx
        IF(left_boundary(rowy-mbod(bod)%dist_y) > colx) left_boundary(rowy-mbod(bod)%dist_y) = colx
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
    right_boundary(i) = ((i+mbod(bod)%dist_y) - 1.0)*nx1 + right_boundary(i)
  END DO
  DO i=1,size(left_boundary)
    IF( left_boundary(i)>MAXINT-1 ) print *, "corresponding left boundary element", i, "is not found"
    left_boundary(i) = ((i+mbod(bod)%dist_y) - 1.0)*nx1 + left_boundary(i)
  END DO
  ! List all the boundary MPM cells and its corresponding freefield cells
  iel_boundary=0; k=1
  DO i=1,size(left_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = left_boundary(i)
    k = k+1
  END DO
  DO i=1,size(left_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = right_boundary(i)
    k = k+1
  END DO
  DO i=1,size(right_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = right_boundary(i)
    k = k+1
  END DO
  DO i=1,size(right_boundary) 
    iel_boundary(1,k) = i
    iel_boundary(2,k) = right_boundary(i) + 1
    k = k+1
  END DO
  DO j=size(left_boundary),size(left_boundary)
    DO i=left_boundary(j), right_boundary(j)
      iel_boundary(1,k) = j
      iel_boundary(2,k) = i
      k = k+1
    END DO
  END DO
  
  !=============================================================================AS
  ! Recreate boundary conditions (determine nf and calculate neq)
  !=============================================================================AS
  g_g=0; nf=0
  ChangeEl: DO bod=1,1
    !- Free all nodes associated with filled element,freeze the freedom of the empty element
    mbod(bod)%g_g=0; mbod(bod)%nf=0

    !-------------------------------------------------------------------------AS
    ! Find GIMP active nodes/elements
    !-------------------------------------------------------------------------AS
    mbod(bod)%GIMP_node_mp=0
    CALL GIMP_activenode(g_num,nip,g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,  &
      mbod(bod)%gimptol,neighb,mbod(bod)%a_ele,mbod(bod)%nf,mbod(bod)%GIMP_node_mp)
    CALL Elem_suport(mbod(bod)%c_ele,mbod(bod)%nf,g_num,cellsize,nx1,nip,nels,  &
      g_coord,mbod(bod)%gm_coord,mbod(bod)%lp_mp,mbod(bod)%gimptol,smethod,mbod(bod)%elemmpoints)

    !-------------------------------------------------------------------------AS
    ! Locate boundary nodes
    !-------------------------------------------------------------------------AS
    lowbound=minval(g_coord(2,:))
    upbound=maxval(g_coord(2,:))
    DO i=1,nn
      IF(g_coord(2,i)<lowbound+cellsize+0.01) mbod(bod)%nf(2,i)=0 
    END DO

    !-------------------------------------------------------------------------AS
    ! Determine equation numbering nf and reallocate memory
    !-------------------------------------------------------------------------AS
    CALL formnf(mbod(bod)%nf)
    mbod(bod)%neq=MAXVAL(mbod(bod)%nf)

    !-------------------------------------------------------------------------AS
    ! Reallocate variable sizes according to the new neq
    !-------------------------------------------------------------------------AS
    CALL deallocate_body(mbod(bod))
    CALL allocate_body(mbod(bod))
    
    !*** AS's Thesis
    ! allocate penalty position marker
    DEALLOCATE(penpos_displacement,penpos_count,penalized_stiffness)
    ALLOCATE(                                 &
      penpos_displacement(0:mbod(bod)%neq),   &
      penpos_count(0:mbod(bod)%neq),          &
      penalized_stiffness(0:mbod(bod)%neq))
  
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

    DEALLOCATE(mbod(bod)%kv,mbod(bod)%mv,mbod(bod)%cv,mbod(bod)%kp,mbod(bod)%mvis)
    ! Allocate skyline matrices
    ALLOCATE(                                           &
      mbod(bod)%kv(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%mv(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%cv(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%kp(mbod(bod)%kdiag(mbod(bod)%neq)),     &
      mbod(bod)%mvis(mbod(bod)%kdiag(mbod(bod)%neq)))
  END DO ChangeEl
  
  END DO time_steps
STOP
END PROGRAM Implicit_MPM_eartquake

