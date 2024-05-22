INCLUDE 'mkl_pardiso.f90'
PROGRAM Tied_Free_Field_Basic
    !---------------------------------------------------------------------------
    ! This code simulates 1D-seismic propagation using the newly depeloped 
    ! boundary conditions Tied-Free-Fields. 
    !---------------------------------------------------------------------------
    USE main
    USE geom
    USE mpm
    USE gimp
    USE sparse_lib
    USE omp_LIB
    USE Logging
    IMPLICIT NONE
    
    INTEGER, PARAMETER::iwp=SELECTED_REAL_KIND(15)

    !** Counters and placeholders
    INTEGER:: i,j,k,m,n,q,s,iel,iters,nsrf,ale,limit,bod,plc,printval,          &
    colx,rowy,ielloc,nnx,nny,tcel,step,nnz,nnzacum,accval,eqno,substep,tstep,   &
    sstep,v_presc,h_presc
      
    REAL(iwp)                 :: dparm(64)
    INTEGER*8                 :: pt(64)   
    INTEGER                   :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
    INTEGER,ALLOCATABLE       :: iparm(:)
    INTEGER, EXTERNAL         :: mkl_get_max_threads
    REAL(iwp)                 :: ddum
    INTEGER  idum(1)

    !** Element/mesh properties
    INTEGER:: ndim=2,ndof=8,nels,neq,nip=4,nn,nod=4,nodf=4,                     & ! nodf is not used in the simulation
    nodof=2,nprops=7,np_types,nst=4,nx1,nx2,nye,ny1,ny2,nlen,                   &
    row,column,newnodes,nxe,cont2,nnxe,nnye,dist_y,aiters,                      &
    count,dist_x,skylength,meshmeth,nels_2D,elastic
        
    REAL(iwp)::h1,h2,s1,w1,w2,lowbound,maxbound,cellsize,Gravf,HSURF,upbound,   &
    prom,mmass,three=3.0_iwp,wff,loadval,valone,penalty=10.0e20_iwp,            &
    leftbound,rightbound

    LOGICAL:: shape=.false.,slope=.false.,newrapconv1

    !** Simulation setup (steps, intervals, filename lengths etc)
    INTEGER:: gravbod1,gravbod2,slopeopt,ploop,accdata
      
    REAL(iwp):: det,tol,dtim,pi,fm,fk,Vpress,Vshear,mu,theta=1
      
    CHARACTER(LEN=15)::element='quadrilateral',argv
      
    LOGICAL::converged,stable=.false.
      
    REAL(iwp),ALLOCATABLE:: nt(:,:),g_matrix(:),ecm(:,:),ecm_acum(:,:),eld_p(:)
      
    REAL(iwp),ALLOCATABLE::kd(:),mm_acum(:,:),mvis(:),Cc(:,:)
                            
    INTEGER::nels_bar,nip_1D=3
    INTEGER::ndof_1D=6,nband,bandval

    !** CMPM/GIMP variables
    INTEGER:: mid,values,bound,aver,boundel

    !** Material parameters/constitutive model
    INTEGER:: fail !Softening parameters
    REAL(iwp)::dsbar,f,lode_theta,sigm,k0,fac,ps1,ps2,ps3,                      &
    dens,Tret,dampsign,pt5=0.5_iwp,Poiss_h,Young_h

    !** Constants/multipliers
    REAL(iwp):: zero=0.0_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp

    !** SYSTEM ARRAYS (FEM arrays for the entire mesh)
    INTEGER,ALLOCATABLE::g(:),g_g(:,:),g_num(:,:),nf(:,:),num(:),MPPE(:),b(:),  &
    g_aux(:),g_s(:),numf(:)

    !MPPE=material points per element
    REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),dee(:,:),eld(:),eps(:),fun(:),   &
    gravlo(:),g_coord(:,:),loads(:),points(:,:),prop(:,:),sigma(:),weights(:),  &
    jac(:,:),der(:,:),deriv(:,:),mpoints(:,:),gc(:),Devstress(:),sp_coord(:,:), &
    lp_coord(:,:),ground_acc(:),mm(:,:),mm_s(:,:),kv_CSR(:),kv_CSR_aux(:),      &
    eps_s(:),km_mf(:,:),mm_mf(:,:),mm_mf_acum(:,:)
      
      
    INTEGER,PARAMETER::ndi=3,nstumat=6
    INTEGER::nshr,nstatv=100,npt,layer,kspt,kstep,kinc,npropsum,void
    REAL(iwp)::sse,spd,scd,rpl,dtime,temp,dtemp,pnewdt,celent
    CHARACTER(LEN=80)::cmname
    REAL(iwp),ALLOCATABLE::ddsddt(:),drplde(:),drpldt(:),stran(:),              &
    time(:),predef(:),dpred(:),props(:),drot(:,:),dfgrd0(:,:),dfgrd1(:,:),      &
    stressumat(:),epsumat(:),deeumat(:,:),statev(:,:)

    !** Data structures
    TYPE::mpm_body
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
    
    !** Define the MPM objects to be used in the simulation
    TYPE(mpm_body):: mbod(3)

    !-----------------------------mpm integers----------------------------------
    ! These are only used for mesh generation. Each body will be 
    ! generated separately, so need to be in body structure.

    INTEGER::mpart,npart,nyp,emps,nmps,tep,newel,slopeel,A,Bi

    !-----------------------------mpm arrays------------------------------------
    INTEGER,ALLOCATABLE::m_num(:),a_ele(:),k_ele(:),b_ele(:),d_ele(:),flag(:)

    REAL(iwp),ALLOCATABLE::gm_coord(:,:),m_coord(:,:),m_mass(:),m_volume(:),    &
    m_stress(:,:),m_velocity(:,:),a_ins(:,:),d1x1(:,:),diag(:,:),d2x1(:),       &
    accp(:,:),acc(:),vcc(:),vccp(:,:),ins(:,:),stress(:),x1(:)
      
    !-----------------------input and initialisation--------------------------
    nlen=7
    argv='Results'
    OPEN(10,FILE='Input/Datafound.dat',status='old')
    OPEN(300,FILE='Input/Groundacc.dat',status='old')
    OPEN(800,FILE='Output/FEM/mpcoord'//'.txt')
    OPEN(11,FILE='Output/FEM/'//argv(1:nlen)//'.txt')

    !--First the variables of the centre body (bod=1) is read
    ! These variables are related to the geometry of the body
    READ(10,*)mbod(1)%name,mbod(1)%slopeopt,mbod(1)%w1,mbod(1)%h1,mbod(1)%s1,   &
    mbod(1)%w2,mbod(1)%h2,mbod(1)%nex,mbod(1)%ney,mbod(1)%nx2,mbod(1)%ny2,      &
    mbod(1)%np_types

    mbod(1)%nx1=mbod(1)%nex
    mbod(1)%ny1=mbod(1)%ney
    cellsize=mbod(1)%w1/mbod(1)%nx1
    mbod(1)%nprops=5

    ALLOCATE(mbod(1)%prop(mbod(1)%nprops,mbod(1)%np_types))
    READ(10,*)mbod(1)%prop
    mbod(1)%ntot=ndof
      
    gravbod1=1
    gravbod2=1
    ploop=2
    meshmeth=1
    v_presc=1
    h_presc=1
    elastic=0
      
    !-The variables of both Free-Fields (bod=2 and bod=3) are read
    !-Number of nodes and elements of each Free-Field body (element)

    ! Define Body Geometry for Free Field Bodies (AS)
    ! Notice that the width is set to 1 and the number of element in x dir is 2
    DO bod=2,size(mbod)
        mbod(bod)%ntot=ndof
        mbod(bod)%nx1=2
        mbod(bod)%w1=1.0_iwp
        mbod(bod)%ny1=mbod(1)%ney
        mbod(bod)%h1=mbod(1)%h1
        mbod(bod)%s1=mbod(1)%s1
        mbod(bod)%w2=mbod(1)%w2
        mbod(bod)%h2=mbod(1)%h2
        ALLOCATE(mbod(bod)%prop(mbod(1)%nprops,mbod(1)%np_types))
        mbod(bod)%prop=mbod(1)%prop
        mbod(bod)%ney=mbod(1)%ny1
        mbod(bod)%nex=2
    END DO
         
    printval=20
      
    ALLOCATE(mpoints(1,2))
      
    READ(300,*)accdata

    ! Create body geometry using mesh_size_slope (AS)
    DO bod=1,size(mbod)
        mbod(bod)%emps=0
        mbod(bod)%nn=0
        mbod(bod)%nels=0
        mbod(bod)%nmps=0
        IF (meshmeth==1) THEN
            CALL mesh_size_slope4(                                              &
                mbod(bod)%nx1, 0,                                               &
                mbod(bod)%ny1, 0,                                               &
                mbod(bod)%nn, mbod(bod)%nels                                    &
            )
        ENDIF
        mbod(bod)%emps=nip
        !--Adding free-field nodes
        mbod(bod)%skylength=500*mbod(bod)%nels
        mbod(bod)%nmps=mbod(bod)%nels*nip
    END DO

    ! Now that nn, nels, and nmps has been defined, allocate all necessary memory (AS)
    count_nodes: DO bod=1,size(mbod)
        ALLOCATE(                                                               &
            mbod(bod)%nf(3,mbod(bod)%nn),                                       &
            mbod(bod)%g_coord(ndim,mbod(bod)%nn),                               &
            mbod(bod)%dee(nst,nst),mbod(bod)%deeinv(nst,nst),                   &
            mbod(bod)%tied_nn(mbod(bod)%nn,2),                                  &
            mbod(bod)%ia_aux(mbod(bod)%skylength),                              &
            mbod(bod)%ja_aux(mbod(bod)%skylength),                              &
            mbod(bod)%g_g(mbod(bod)%ntot,mbod(bod)%nels),                       &
            mbod(bod)%g(mbod(bod)%ntot),                                        &
            mbod(bod)%g_coord_aux(ndim,mbod(bod)%nn),                           &
            mbod(bod)%km(mbod(bod)%ntot,mbod(bod)%ntot),                        &
            mbod(bod)%KGC(mbod(bod)%ntot,mbod(bod)%ntot),                       &
            mbod(bod)%MMS(mbod(bod)%ntot,mbod(bod)%ntot),                       &
            mbod(bod)%CCS(mbod(bod)%ntot,mbod(bod)%ntot),                       &
            mbod(bod)%MOD_MTX(mbod(bod)%ntot,mbod(bod)%ntot)                    &
        )
        IF(bod>=1) ALLOCATE(mbod(bod)%g_num(nod,mbod(bod)%nels))
    END DO count_nodes
       
    ! These allocations are for variables that are not unique to the body but
    ! necessary for doing calculations
    ALLOCATE(                                                                   &
        weights(nip),points(nip,ndim),num(nod),coord(nod,ndim),fun(nod),        &
        bee(nst,ndof),eps(nst),sigma(nst),jac(ndim,ndim),der(ndim,nod),         &
        eld(ndof),deriv(ndim,nod),stress(nst),sp_coord(ndim,1),                 &
        lp_coord(1,ndim),g_s(ndof),eld_p(nodf),numf(nodf),eps_s(nst),           &
        g_aux(16),mm_acum(ndof_1D,ndof_1D),km_mf(ndof,ndof_1D),Cc(2,2),         &
        mm_mf(ndof,ndof_1D),mm_mf_acum(ndof,ndof_1D)                            &
    )   

    ALLOCATE(                                                                   &
        m_num(nip),vcc(ndof),m_coord(nip,ndim),                                 &
        gc(ndim),ecm(ndof,ndof),ecm_acum(ndof,ndof)                             &
    )

    ! read ground motion data
    ALLOCATE(ground_acc(accdata))
    READ(300,*)ground_acc

    ! read material properties
    OPEN(400,FILE='Input/parumat.dat',status='old')
    READ(400,*)npropsum
    ALLOCATE(                                                                   &
        ddsddt(nstumat),drplde(nstumat),stran(nstumat),                         &
        props(npropsum),stressumat(nstumat),epsumat(nstumat),                   &
        deeumat(nstumat,nstumat)                                                &
    )
    READ(400,*)props

    ! Now we want to get the actual coordinates of the mesh (AS)
    Mesh_create: DO bod=1,size(mbod)
        IF(bod>=1)THEN     
            mbod(bod)%nf=0
            mbod(bod)%g_num=0
            mbod(bod)%g_coord=0
            row=1
            column=1
            A=0
            Bi=0
            dist_x=0

            !----- loop the elements to find global arrays sizes -----
            colx=1
            rowy=1
            coord=zero
            num=0
            
            ! For every element, get the node coordinates and node number
            ! Then, saved it to the structure variable 
            elements_2D: DO iel=1,mbod(bod)%nels
                    !- coordenates and numbering main body   
                    Call emb_2d_geom_4n(                                        &
                        iel,mbod(bod)%nx1,                                      &
                        mbod(bod)%nx2,                                          &
                        mbod(bod)%ny1,                                          &
                        mbod(bod)%ny2,                                          &
                        mbod(bod)%w1,                                           &
                        mbod(bod)%s1,                                           &
                        mbod(bod)%w2,                                           &
                        mbod(bod)%h1,                                           &
                        mbod(bod)%h2,                                           &
                        coord,num                                              &
                    )
                mbod(bod)%g_num(:,iel)=num
                mbod(bod)%g_coord(:,num)=TRANSPOSE(coord)

                ! The generated geometry centers at 0. In order to not overlap,
                ! the main body and the right free-field body needs to be displaced
                ! for the main body, displace the whole coordinate to the right
                ! with a value of left-ff-width + one, for the right boundary 
                ! displace the whole body with a value of 
                ! left-ff-width + one + main-body-width + one
                IF(bod==1)THEN ! main body
                    valone=mbod(2)%w1+one 
                ELSE IF(bod==2)THEN ! left ff body
                    valone=zero
                ELSE ! right ff body
                    valone=mbod(2)%w1+mbod(1)%nex+2.0_iwp
                END IF
                mbod(bod)%g_coord(1,num)=mbod(bod)%g_coord(1,num)+valone
                
                ! increment rowy to get the next element in the y-dir
                rowy=rowy+1

                ! if its over the boundary, reset rowy, and start with new column
                IF(rowy>mbod(bod)%ney.and.colx<=mbod(bod)%nx1)THEN
                    rowy=1
                    colx=colx+1
                END IF

                IF(colx>mbod(bod)%nx1.and.rowy==1)THEN
                    rowy=mbod(bod)%ny1+1
                END IF  
                
                IF(colx>mbod(bod)%nx1.and.rowy>mbod(bod)%ney)THEN
                    rowy=mbod(bod)%ny1+1
                    colx=colx+1
                END IF
                upbound=maxval(mbod(bod)%g_coord)
            END DO elements_2D
        END IF
    END DO Mesh_create
  
    ! ------------------------- insert material points -------------------------

    ! In this section, the integration point of each element is defined
    Gauss_vectors: DO bod=1,size(mbod)

        ALLOCATE(mbod(bod)%MPPE(nip))
        
        ALLOCATE(                                                               &
            mbod(bod)%gm_coord(ndim,mbod(bod)%nmps),                            &
            mbod(bod)%m_volume(mbod(bod)%nmps),                                 &
            mbod(bod)%m_mass(mbod(bod)%nmps),                                   &
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
            mbod(bod)%fk(mbod(bod)%nmps)                                        &
        )

        mbod(bod)%a_ele=0
        mbod(bod)%gm_coord=0.0_iwp

        !---------------- Gauss points weights and coordenates -----------------
        
        
        ! for each element, define the integration locations and assign
        ! a numbering sequence based on MPPE
        mbod(bod)%MPPE=0
        DO i=1,nip
            mbod(bod)%MPPE(i)=i
        END DO

        CALL sample(element,points,weights) !should be 2 equaled spaced points
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
            ! get the next gaussian point indexes for the next element
            DO j=1,nip
                mbod(bod)%MPPE(j)=mbod(bod)%MPPE(j)+nip
            END DO
        END DO InsertMP
    
        !--------- get initial volume and mass for each material point ---------
        A=0
        k=1   
        
        DO i=1,mbod(bod)%nmps
            iel                   = mbod(bod)%a_ele(i)
            num                   = mbod(bod)%g_num(:,iel)
            coord                 = TRANSPOSE(mbod(bod)%g_coord(:,num)) 

            CALL shape_der(der,points,k)
            jac                   = MATMUL(der,coord)
            det                   = determinant(jac)

            ! Notice here the mp volume, density, and mass are calculated
            mbod(bod)%m_volume(i) = det*weights(k)
            mbod(bod)%m_dens(i)   = ((mbod(bod)%prop(1,1))/(one+props(16)))            
            mbod(bod)%m_mass(i)   = mbod(bod)%m_dens(i)*mbod(bod)%m_volume(i)

            k=k+1
            IF(k>nip)k=1
        END DO
    END DO Gauss_vectors
      
    !--------------------------------------   
    !--CREATE FREE-FIELD steering vectors
    !---------------------------------------       

    ALLOCATE(                                                                   &
        g_matrix(ndim),                                                         &
        acc(ndof),                                                              &
        nt(ndof,nodof),                                                         &
        mm(ndof,ndof),                                                          &
        mm_s(ndof,ndof)                                                         &
    )

    ! ------- find out in which element the material point is located in -------

    Flags: DO bod=1,size(mbod)
        IF(mbod(bod)%nmps>1)THEN 
            ielloc        = 1
            m             = 0
            mbod(bod)%ale = 1
            DO i=1,mbod(bod)%nmps
                m = m+1
                IF(m>nip)THEN 
                    ielloc = ielloc+1
                    m      = 1
                END IF
                inner: DO iel=1,mbod(bod)%nels
                    converged=.true.
                    IF(converged) THEN
                        mbod(bod)%a_ele(i) = ielloc
                        num                = mbod(bod)%g_num(:,ielloc)
                        coord              = TRANSPOSE(mbod(bod)%g_coord(:,num))
                        sp_coord(:,1)      = mbod(bod)%gm_coord(:,i)
                        lp_coord(1,:)      = mbod(bod)%mpoints(i,:)

                        ! get the local coordinates of each integration points
                        CALL floc(coord,sp_coord,lp_coord,i)
                        mbod(bod)%mpoints(i,:)=lp_coord(1,:)

                        ! This part is for early exit and counting the ammount of
                        ! active elements in mbod%ale
                        IF(i>1)THEN
                            DO j=1,i-1
                                IF(mbod(bod)%a_ele(j)==mbod(bod)%a_ele(i)) THEN
                                    EXIT inner
                                END IF
                            END DO
                            mbod(bod)%ale=mbod(bod)%ale+1
                        END IF
                        
                        EXIT inner
                    END IF
                END DO inner
            END DO
        END IF
    END DO Flags

    !----------------------- Global boundary conditions ------------------------

    Body_fixities: DO bod=1,size(mbod)
        ! mark nf to 0 for non displacing boundary
        mbod(bod)%nf = 1
        mbod(bod)%nf(3,:) = 0 !No pore pressure
        ! Find the boundary of the geometry from g_coord
        lowbound = minval(mbod(bod)%g_coord)
        upbound  = maxval(mbod(bod)%g_coord)    ! Computed earlier to consider 
                                                ! only main domain coordinates
        
        ! g_coord_aux save the mesh coordinates in the initial condition
        mbod(bod)%g_coord_aux = mbod(bod)%g_coord

        ! Setup vertical boundary conditions
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
    
        j=1; k=1; m=0
        
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
        ALLOCATE(                                                               &
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
        
        ALLOCATE(                                                               &
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
            mbod(bod)%g(1:ndof:2)  = mbod(bod)%nf(1,num(:)) ! this one in x-dir
            mbod(bod)%g(2:ndof:2)  = mbod(bod)%nf(2,num(:)) ! this one in y-dir
            mbod(bod)%g_g(:,iel) = mbod(bod)%g
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
        ALLOCATE(                                                               &
            mbod(bod)%c_matrix(mbod(bod)%neq,2*(nband+1)),                      &
            mbod(bod)%mf_matrix(mbod(bod)%neq,2*(nband+1))                      &
        )
        
        j=1
        m=1
        mbod(bod)%kdiag=zero
        
        ! for each element, determine the diagonal index using fkdiag
        DO iel=1,mbod(bod)%nels
            mbod(bod)%g=mbod(bod)%g_g(:,iel) 
            CALL fkdiag(mbod(bod)%kdiag,mbod(bod)%g)
        END DO
        
        ! the step before only valid for each element, to build the whole matrix
        ! acummulate sum the diagonal terms
        DO i=2,mbod(bod)%neq
            mbod(bod)%kdiag(i)=mbod(bod)%kdiag(i)+mbod(bod)%kdiag(i-1)
        END DO
    
        ALLOCATE(                                                               &
            mbod(bod)%kv_CSR_aux(mbod(bod)%skylength),                          &
            mbod(bod)%mv(mbod(bod)%kdiag(mbod(bod)%neq)),                       &
            mbod(bod)%mvis(mbod(bod)%kdiag(mbod(bod)%neq))                      &
        )
   
    END DO Body_fixities

    !--------------------variables initialisation ------------------------------
    READ(10,*)dtim,k0,aiters,tol,limit,nsrf

    !--- Global initial conditions ---
    d1x1  = zero
    pi    = ACOS(-one)     
    cont2 = 0
    
    !--Body initial conditions--
    initial_conditions: DO bod=1,size(mbod)
        mbod(bod)%m_acc         = zero
        mbod(bod)%m_stress      = zero
        mbod(bod)%m_stress_efe  = zero
        mbod(bod)%m_pore        = zero
        mbod(bod)%m_stress_ini  = zero
        mbod(bod)%accp          = zero
        mbod(bod)%a_ins         = zero
        mbod(bod)%x1            = zero
        mbod(bod)%x1_acum       = zero
        mbod(bod)%d1x1          = zero
        mbod(bod)%d2x1          = zero
        mbod(bod)%accb          = zero
        mbod(bod)%vccp          = zero
        mbod(bod)%flag          = 0
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
        accval                  = 0
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
    END DO initial_conditions
           
    Gravf = 10.0_iwp
    HSURF = 0.0
    k0    = 0.5_iwp

    !------------------------------------
    !-Initial stresses (k0 procedure)----
    !------------------------------------
    DO bod=1,size(mbod)
        mbod(bod)%g_matrix=(/zero,-10.0_iwp/)
        Ini_stress_loop:DO i=1,mbod(bod)%nmps
            mbod(bod)%m_stress_efe(2,i) = mbod(bod)%m_dens(i)*10.0_iwp *        &
                                          (mbod(bod)%gm_coord(2,i)-HSURF)       
            mbod(bod)%m_stress_efe(1,i) = mbod(bod)%m_dens(i)*10.0_iwp*         &
                                          (mbod(bod)%gm_coord(2,i)-HSURF)*k0
            mbod(bod)%m_stress_efe(3,i) = zero
            mbod(bod)%m_stress_efe(4,i) = mbod(bod)%m_stress_efe(1,i)
                
            mbod(bod)%m_stress(2,i) = mbod(bod)%m_dens(i)*10.0_iwp *            &
                                      (mbod(bod)%gm_coord(2,i)-HSURF)
            mbod(bod)%m_stress(1,i) = mbod(bod)%m_stress(2,i)*k0
            mbod(bod)%m_stress(3,i) = zero
            mbod(bod)%m_stress(4,i) = mbod(bod)%m_stress(1,i)
        END DO Ini_stress_loop
        mbod(bod)%m_stress_ini = mbod(bod)%m_stress
        CONTINUE
    END DO
     
    DO i=1,mbod(1)%nmps
        WRITE(800,'(3f15.6)')mbod(1)%gm_coord(1,i),mbod(1)%gm_coord(2,i)
    END DO

    g_matrix=(/zero,-10.0_iwp/)  !--Gravity acting in the vertical direction
    fm=zero!0.0616_iwp*2.0_iwp   ! Reighleigh Mass Coef
    fk=zero!0.003121_iwp*2.0_iwp ! Reighleigh Stiffness Coef
    substep=1
      
    dtim=dtim/substep
     
    IF(gravbod1==2) mbod(1)%g_matrix = zero

    stable=.true.

    !!==========================================================================
    !! -------------------- BEGINNING OF THE TIME STEP -------------------------
    !!==========================================================================

    step=0
    time_steps: DO tstep=1,2000
        PRINT '("Steps :" (I5) "/" (I4))', tstep, 2000
        sub_step: DO sstep =1,substep
            step=step+1 
            wff=1.0_iwp

            CALL sample(element,points,weights)
            Body_Solution: DO bod=1,size(mbod)
                
                mbod(bod)%gravlo=zero
                mbod(bod)%diag=zero
                mbod(bod)%f_fint=zero
                mbod(bod)%vcm=zero
                mbod(bod)%c_force=zero
                mbod(bod)%mf_force=zero
                mbod(bod)%loads=zero
                mbod(bod)%x1=zero
                mbod(bod)%kinup_d1x1=zero
                mbod(bod)%kinup_d2x1=zero
                newrapconv1=.false.
                mbod(bod)%loads_end=zero
                mbod(bod)%loads_base=zero
                mbod(bod)%kinup_Ground_d2x1=zero
                mbod(bod)%ground_loads=zero
                mbod(bod)%penpos=0
                mbod(bod)%penpos_v=0
                mbod(bod)%eq_pen=0
                mbod(bod)%eq_pen_v=0
                        
                !-------------------------------
                !-External forces Loop (gravity)
                !-------------------------------
                i=1
                j=1
                n=1
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
                    mbod(bod)%gravlo(mbod(bod)%g(1:ndof)) =                       &
                        mbod(bod)%gravlo(mbod(bod)%g(1:ndof)) +                   &
                        mbod(bod)%diag(mbod(bod)%g(1:ndof))*(-Gravf)  
                
                    mbod(bod)%gravlo(0)=zero

                    ! reset mm_s for next integration point
                    ! the index is to cycle through the weigths
                    ! Note: this will break if a_ele is not ordered correctly
                    mm_s=zero
                    i=i+1
                    IF(i>nip)i=1
            
                END DO Fext_domain_ff
            END DO Body_Solution

            !-------------------------------------------------------------------
            ! This Loop is to create the system M+C+K (Mass+Damping+Stiffness)
            !-------------------------------------------------------------------
            
            CALL sample(element,points,weights)
            Stiffness_2D:DO bod=1,size(mbod)
        
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
                        
                    ! Create the element stiffness matrix
                    CALL deemat(mbod(bod)%dee,mbod(bod)%Young,mbod(bod)%Poiss)
                    mbod(bod)%km = mbod(bod)%km +                               &
                        MATMUL(MATMUL(TRANSPOSE(bee),mbod(bod)%dee),bee)*det*weights(i)

                    ! Recreate the mass matrix
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

                    ! at the end of each element form the modified stiffness matrix
                    ! according to the Newmark Time Integration approach
                    IF(i>nip-1)THEN
                        mbod(bod)%KGC = zero
                        mbod(bod)%MMS = zero
                        mbod(bod)%KGC = mbod(bod)%km
                        mbod(bod)%MMS = mm_s
                        mbod(bod)%MOD_MTX = mbod(bod)%KGC +                     &
                            4.0_iwp*mbod(bod)%MMS/dtim**2.0_iwp +               &
                            2.0_iwp*(mbod(bod)%KGC*fk +                         &
                            mbod(bod)%MMS*fm)/dtim
                            
                        CALL formspars_unsym(                                   &
                            mbod(bod)%ntot,                                     &
                            mbod(bod)%g,                                        &
                            mbod(bod)%MOD_MTX,                                  &
                            mbod(bod)%ia_aux,                                   &
                            mbod(bod)%ja_aux,                                   &
                            mbod(bod)%kv_CSR_aux,                               &
                            mbod(bod)%skylength                                 &
                        )
                        mm_acum    = zero
                        mm_mf_acum = zero
                        mm_s       = zero
                    END IF    
                    i=i+1
                    IF(i>nip)i=1
                END DO KM_MV_2D
            END DO Stiffness_2D 

            !-------------------------------------------------------------------
            ! This Loop is to preparing the solvers by constructing the CSR Form
            !-------------------------------------------------------------------
            
            DO bod=1,size(mbod)
                nnzacum=0
                CALL sortadd(                                                   &
                    mbod(bod)%skylength,                                        &
                    mbod(bod)%ia_aux,                                           &
                    mbod(bod)%ja_aux,                                           &
                    mbod(bod)%kv_CSR_aux,                                       &
                    mbod(bod)%neq+1,                                            &
                    nnzacum,                                                    &
                    mbod(bod)%penpos                                            &
                )
                ALLOCATE(                                                       &
                    mbod(bod)%ia(mbod(bod)%neq+1),                              &
                    mbod(bod)%ja(nnzacum),                                      &
                    mbod(bod)%kv_CSR(nnzacum)                                   &
                )
                mbod(bod)%kv_CSR = zero;mbod(bod)%ja=0
                mbod(bod)%ia     = 0
                mbod(bod)%kv_CSR(1:nnzacum)     = mbod(bod)%kv_CSR_aux
                mbod(bod)%ja(1:nnzacum)         = mbod(bod)%ja_aux
                mbod(bod)%ia(1:mbod(bod)%neq+1) = mbod(bod)%ia_aux(1:mbod(bod)%neq+1)
            END DO
        
            !-------------------------------------------------------------------
            !-Find equation number to apply penalty method (presc. displacement)
            ! <h_presc = 1> and <v_presc = 1> indicate that the vertical and/or 
            ! horizontal boundary nodes will have prescribed displacements
            !-------------------------------------------------------------------

            DO bod=1,1
                m = 0
                j = 0
                DO i=1,mbod(bod)%nn-1
                    ! check whether the node is a boundary node
                    IF (                                                        &
                        (i<=mbod(bod)%ny1) .or.                               &
                        (i>(mbod(bod)%ny1+1)*(mbod(bod)%nex))               &
                    ) THEN
                    ! apply prescribed displacements on respective directions
                        IF(h_presc==1)THEN
                            mbod(bod)%eq_pen(m+1)=mbod(bod)%nf(1,i)
                            m=m+1
                        END IF
                        IF(v_presc==1)THEN
                            mbod(bod)%eq_pen(m+1)=mbod(bod)%nf(2,i)
                            m=m+1
                        END IF
                    END IF     
                END DO
        
                !---------------------------------------------------------------
                !- This Do-Loop is to prescribe a displacement to 
                !  the lower boundary nodes
                !---------------------------------------------------------------
        
                DO i=1,mbod(bod)%nn
                    IF(mbod(bod)%g_coord(2,i)<lowbound+0.01)THEN
                        mbod(bod)%eq_pen_v(j+1)=mbod(bod)%nf(1,i)
                        j=j+1
                    END IF    
                END DO

                !---------------------------------------------------------------
                !- The following Do-Loop is to find the equations associated to 
                !  the stiffnes matrix diagonal in the horizontal 
                !  degree of freedom (DOF) mbod(bod)%penpos and 
                !  vertical DOF mbod(bod)%penpos_v 
                !---------------------------------------------------------------
            
                mbod(bod)%penpos=0
                
                DO i=1,m
                    n=mbod(bod)%eq_pen(i)
                    IF(n>0)THEN
                        k=mbod(bod)%ia(n)
                        DO q=1,500
                            IF(mbod(bod)%ja(k).ne.n)THEN
                                k=k+1
                            ELSE
                                mbod(bod)%penpos(i)=k
                                EXIT
                            END IF
                        END DO
                    END IF
                END DO
                
                mbod(bod)%penpos_v=0
                DO i=1,j
                    n=mbod(bod)%eq_pen_v(i)
                    IF(n>0)THEN
                        k=mbod(bod)%ia(n)
                        DO q=1,500
                            IF(mbod(bod)%ja(k).ne.n)THEN
                                k=k+1
                            ELSE
                                mbod(bod)%penpos_v(i)=k
                                EXIT
                            END IF
                        END DO
                    END IF
                END DO
                
                !----------------------------------------------  
                !- Apply penalty values to the stiffness matrix
                !----------------------------------------------  
                
                DO i=1,m
                    mbod(bod)%kv_CSR(mbod(bod)%penpos(i)) =                     &
                        mbod(bod)%kv_CSR(mbod(bod)%penpos(i)) + penalty
                END DO
                
                DO i=1,j
                    mbod(bod)%kv_CSR(mbod(bod)%penpos_v(i)) =                   &
                        mbod(bod)%kv_CSR(mbod(bod)%penpos_v(i)) + penalty
                END DO
            END DO
        
            !-------------------------------------------------------------------
            !- End penalty fixities
            !-------------------------------------------------------------------
                
                
            !-------------------------------------------------------------------
            !-The following Do-Loop is to create the Mass and 
            ! Damping force vectors F=Ma and F=Cv
            !-------------------------------------------------------------------
            
            DO bod=1,size(mbod)
                mbod(bod)%mv        = zero
                mbod(bod)%c_matrix  = zero
                mbod(bod)%mvis      = zero
                mbod(bod)%mf_matrix = zero
                mbod(bod)%km        = zero
                mbod(bod)%CCS       = zero
                ecm                 = zero
                i = 1
                n = 1
                j = 1
                q = 1
                mm_s = zero

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

                    mbod(bod)%km = mbod(bod)%km +                               &
                        MATMUL(MATMUL(TRANSPOSE(bee),mbod(bod)%dee),bee)*det*weights(i)
                    
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
        

            !-------------------------------------------------------------------
            !-After creating the Ma+Cv+Ku matrix system, the following "Do loop" 
            !-decompose the system into a LU system to perform matrix operations 
            !-------------------------------------------------------------------
        
            DO bod=1,size(mbod) 
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
                !iparm(3) = 4 
                
                phase = 11  !// check matrix consistency
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

    !---------------------------------------------------------------------------
    !---------------------------- ITERATION LOOP -------------------------------
    !---------------------------------------------------------------------------
            
            iters=0 

            Newton_Rhapson: DO WHILE(iters<=limit)        
                iters=iters+1

                !---------------------
                !-Internal forces Loop
                !---------------------
                
                ! first get the nodal force by integrating the stresses
                DO bod=1,size(mbod)
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
                        mbod(bod)%ddylds(mbod(bod)%g(1:ndof)) =                   &
                            mbod(bod)%ddylds(mbod(bod)%g(1:ndof)) +               &
                            MATMUL(TRANSPOSE(bee),sigma)*det*weights(a)
                        a=a+1
                        IF(a>nip)a=1
                        mbod(bod)%ddylds(0)=zero 
                    END DO Fint_load
                END DO 


                DO bod=1,size(mbod)
                    IF(tstep<=accdata.and.tstep>=1)THEN
                        IF(bod==1)THEN
                            ! This part is used to calculate the main body
                            DO i=1,mbod(bod)%nn
                                IF(mbod(bod)%g_coord(2,i)<lowbound+0.01_iwp)THEN
                                    valone = ground_acc(tstep)
                                    ! add ground force : f = 2*a - 
                                    mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i))=                     &
                                        mbod(bod)%kinup_Ground_d2x1(mbod(bod)%nf(1,i))+                 &
                                        (2.0_iwp*ground_acc(tstep)-mbod(bod)%kinup_d2x1(mbod(bod)%nf(1,i)))
                                    mbod(bod)%kinup_Ground_d2x1(0)=zero
                                END IF    
                            END DO  
                        ELSE
                            ! This part is used for calculate the free field bodies
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
                    END IF
                END DO

                !-Kinematic update of acceleration and velocity without 
                ! considering any boundary condition (contact)

                DO bod=1,size(mbod)
                    mbod(bod)%kinup_d2x1   = (4.0_iwp*mbod(bod)%x1/dtim**2.0_iwp)- &
                                             (4.0_iwp*mbod(bod)%d1x1/dtim)-     &
                                             mbod(bod)%d2x1         
                    mbod(bod)%kinup_d1x1   = 2.0_iwp*mbod(bod)%x1/dtim-         &
                                             mbod(bod)%d1x1
                    mbod(bod)%kinup_d2x1(0) = zero
                    mbod(bod)%kinup_d1x1(0) = zero
                    mbod(bod)%d1p1(0)       = zero
                END DO 

                Force_MA: DO bod=1,size(mbod)
                    ! This is force from the damping F=Cv
                    CALL bantmul(                                               &
                        mbod(bod)%c_matrix,                                     &
                        mbod(bod)%kinup_d1x1,                                   &
                        mbod(bod)%c_force                                       &
                    )
                    ! This is force from the general body motion F=Ma
                    CALL bantmul(                                               &
                        mbod(bod)%mf_matrix,                                    &
                        mbod(bod)%kinup_d2x1,                                   &
                        mbod(bod)%mf_force                                      &
                    )
                    mbod(bod)%vcm(0)      = zero
                    mbod(bod)%c_force(0)  = zero
                    mbod(bod)%mf_force(0) = zero
                END DO Force_MA 

                Ground_F: DO bod=1,size(mbod)
                    mbod(bod)%kinup_Ground_d2x1(0)=zero
                    ! This is force from the ground motion also F=Ma
                    CALL bantmul(                                               &
                        mbod(bod)%mf_matrix,                                    &
                        mbod(bod)%kinup_Ground_d2x1,                            &
                        mbod(bod)%ground_loads                                  &
                    )
                    mbod(bod)%ground_loads(0)=zero
                END DO Ground_F

                CONTINUE
                
                DISPLACEMENTS: DO bod=size(mbod),1,-1
                    mbod(bod)%loads    = zero
                    mbod(bod)%residual = zero
                    ! Combine all the forces: gravity loads, internal loads, ground loads, 
                    ! displacement loads, and damping load
                    mbod(bod)%loads    = mbod(bod)%gravlo - mbod(bod)%ddylds +  &
                                         mbod(bod)%ground_loads -               &
                                         mbod(bod)%mf_force - mbod(bod)%c_force
                
                    mbod(bod)%loads(0)=zero

                    ! Check whether the calculated loads has converge (saved in newrapconv1)
                    IF(iters==1) mbod(bod)%loads_ini = mbod(bod)%loads
                    IF(iters>1.and.bod==1) THEN
                        CALL checon_1(                                          &
                            mbod(bod)%loads_ini(1:mbod(bod)%neq),               &
                            mbod(bod)%loads(1:mbod(bod)%neq),                   &
                            newrapconv1,                                        &
                            tol,                                                &
                            mbod(bod)%neq                                       &
                        ) 
                    END IF

                    !-----------------------------------------------------------
                    !- Apply (prescribed) free-fields displacements (bod 2 & 3) 
                    !  to the boundaries of the main domain (bod 1) 
                    !-----------------------------------------------------------

                    IF(bod==1)THEN
                        j=1

                        ! This one corresponds to the left boundary
                        m=1
                        DO i=1,mbod(bod)%nn
                            IF(i<=mbod(bod)%ny1)THEN
                                IF(h_presc==1)THEN
                                    mbod(bod)%loads(mbod(bod)%nf(1,i)) =        &
                                        mbod(bod)%kv_CSR(mbod(bod)%penpos(j)) * &
                                        mbod(2)%loads(mbod(2)%nf(1,m))
                                    j=j+1
                                END IF
                                IF(v_presc==1)THEN
                                    mbod(bod)%loads(mbod(bod)%nf(2,i)) =        &
                                        mbod(bod)%kv_CSR(mbod(bod)%penpos(j))*  &
                                        mbod(2)%loads(mbod(2)%nf(2,m))
                                    j=j+1
                                END IF
                            m=m+1
                            END IF  
                        END DO
                        
                        ! this one corresponds to the right boundary
                        m=1
                        DO i=1,mbod(bod)%nn-1
                            IF(i>(mbod(bod)%ny1+1)*(mbod(bod)%nex))THEN
                                IF(h_presc==1)THEN
                                    mbod(bod)%loads(mbod(bod)%nf(1,i)) =        &
                                        mbod(bod)%kv_CSR(mbod(bod)%penpos(j)) * &
                                        mbod(3)%loads(mbod(3)%nf(1,m))
                                    j = j+1
                                END IF
                                IF(v_presc==1)THEN
                                    mbod(bod)%loads(mbod(bod)%nf(2,i)) =        &
                                        mbod(bod)%kv_CSR(mbod(bod)%penpos(j)) * &
                                        mbod(3)%loads(mbod(3)%nf(2,m))
                                    j = j+1
                                END IF
                                m = m+1
                            END IF
                        END DO
                        
                        ! this one corresponds to the bottom boundary
                        m=1
                        DO i=1,mbod(bod)%nn
                            IF(mbod(bod)%g_coord(2,i)<lowbound+0.01_iwp)THEN
                                mbod(bod)%loads(mbod(bod)%nf(1,i)) =            &
                                    mbod(bod)%kv_CSR(mbod(bod)%penpos_v(m)) *   &
                                    mbod(2)%loads(mbod(2)%nf(1,mbod(1)%ny1+1))
                                m = m+1
                            END IF    
                        END DO
                    END IF
            
                    !-----------------------------------------------------------
                    !------- Solve the system of equations Ma+Cv+Ku = F --------
                    !-----------------------------------------------------------

                    mbod(bod)%iparm(8)  = 10 ! max numbers of iterative refinement steps
                    
                    ! pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
                    ! A*X = B
                    ! X = mbod(bod)%loads(1:mbod(bod)%neq)
                    ! B = mbod(bod)%residual(1:mbod(bod)%neq)
                    phase = 33  ! solve equations
                    call pardiso(                                               &
                        mbod(bod)%pt,                                           &
                        maxfct,                                                 &
                        mnum,                                                   &
                        mtype,                                                  &
                        phase,                                                  &
                        mbod(bod)%neq,                                          &
                        mbod(bod)%kv_CSR,                                       &
                        mbod(bod)%ia,                                           &
                        mbod(bod)%ja,                                           &
                        idum,                                                   &
                        nrhs,                                                   &
                        mbod(bod)%iparm,                                        &
                        msglvl,                                                 &
                        mbod(bod)%loads(1:mbod(bod)%neq),                       &
                        mbod(bod)%residual(1:mbod(bod)%neq),                    &
                        error                                                   &
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

                END DO DISPLACEMENTS

                ! calculate the acceleration and velocity increments
                DO bod=1,size(mbod)
                    ! acceleration
                    mbod(bod)%kinup_d2x1    = (4.0_iwp*mbod(bod)%x1/dtim**2.0_iwp) - &
                                              (4.0_iwp*mbod(bod)%d1x1/dtim) -        &
                                              mbod(bod)%d2x1
                    ! velocity
                    mbod(bod)%kinup_d1x1    = 2.0_iwp*mbod(bod)%x1/dtim-mbod(bod)%d1x1
                    mbod(bod)%kinup_d2x1(0) = zero
                    mbod(bod)%kinup_d1x1(0) = zero
                    mbod(bod)%d1p1(0)       = zero
                END DO 

                ! update stresses and check yield criterion
                DO bod=1,size(mbod)
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
              END DO Newton_Rhapson

            ! reset the convergence state for next time step
            newrapconv1=.false.
            IF(sstep==substep) cont2 = cont2+1

            ! x1_acum is here just for post-simulation purpose. Not for actual calculation
            DO bod=1,size(mbod)
                mbod(bod)%x1_acum = mbod(bod)%x1_acum+mbod(bod)%x1
                mbod(bod)%d2x1    = mbod(bod)%kinup_d2x1
                mbod(bod)%d1x1    = mbod(bod)%kinup_d1x1
            END DO
        
            ! Update the node displacement
            !- NOTE: Changing nodes coordenates can damage simulations 
            DO bod=1,size(mbod)
                DO i=1,mbod(bod)%nn
                    mbod(bod)%g_coord_aux(1,i) = mbod(bod)%g_coord_aux(1,i) +   &
                                                 mbod(bod)%x1(mbod(bod)%nf(1,i))
                    mbod(bod)%g_coord_aux(2,i) = mbod(bod)%g_coord_aux(2,i) +   &
                                                 mbod(bod)%x1(mbod(bod)%nf(2,i))
                END DO
            END DO
        
            ! These are all variables for presentation and update the location
            ! of the integration point
            REVERSE_MAPPING: DO bod=1,size(mbod) 
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
            END DO REVERSE_MAPPING
        
            ! The following Do-Loop is necessary to erase the accumulated memory 
            ! associated to the matrices used in the pardiso solver
        
            phase = -1
            DO bod=1,size(mbod)
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

            IF(cont2/printval*printval==cont2.or.tstep==1)THEN
                DO bod=1,size(mbod)
                    CALL paraview2(                                             &
                        cont2,                                                  &
                        1,                                                      &
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
                        nels_2D,                                                &
                        mbod(bod)%penpos,bod                                    &
                    )     
                    CALL point_viz2(                                            &
                        cont2,                                                  &
                        1,                                                      &
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
                        bod                                                     &
                    )
                END DO
            END IF
        END DO sub_step
    END DO time_steps
       
    STOP

END PROGRAM Tied_Free_Field_Basic

 