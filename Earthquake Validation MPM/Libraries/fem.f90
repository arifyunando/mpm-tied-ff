MODULE fem
  SAVE
CONTAINS
  SUBROUTINE ecmat2(ecm,fun,ndof,nodof)
    !
    ! This subroutine forms the element consistent mass matrix.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::fun(:)
    REAL(iwp),INTENT(OUT)::ecm(:,:)
    INTEGER,INTENT(IN)::nodof,ndof
    INTEGER::nod,i,j
    REAL::nt(ndof,nodof),tn(nodof,ndof),zero=0.0_iwp
    ecm=zero
    nod=ndof/nodof
    nt=zero
    tn=zero
    DO i=1,nod
      DO j=1,nodof
        nt((i-1)*nodof+j,j)=fun(i)
        tn(j,(i-1)*nodof+j)=fun(i)
      END DO
    END DO
    
    ecm=MATMUL(nt,tn)
    RETURN
  END SUBROUTINE ecmat2

  SUBROUTINE bantmul(kb,loads,ans)
    !
    ! This subroutine multiplies an unsymmetrical band kb by the vector loads.
    ! Could be much improved for vector processors.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::kb(:,:),loads(0:)
     REAL(iwp),INTENT(OUT)::ans(0:)
     INTEGER::i,j,k,l,m,n,iw
     REAL(iwp)::x,zero=0.0_iwp
     n=SIZE(kb,1)
     l=SIZE(kb,2)
     iw=(l-1)/2
     DO i=1,n
       x=zero
       k=iw+2
       DO j=1,l
         k=k-1
         m=i-k+1
         IF(m<=n.AND.m>=1)x=x+kb(i,j)*loads(m)
       END DO
       ans(i)=x
     END DO
    RETURN
  END SUBROUTINE bantmul
  
  SUBROUTINE formlump(diag,emm,g)
    !
    ! This subroutine forms the lumped global mass matrix as a vector diag.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::emm(:,:)
     REAL(iwp),INTENT(OUT)::diag(0:)
     INTEGER,INTENT(IN)::g(:)
     INTEGER::i,ndof
     ndof=UBOUND(emm,1)
     DO i=1,ndof
       diag(g(i))=diag(g(i))+emm(i,i)
     END DO
    RETURN
  END SUBROUTINE formlump
  
  SUBROUTINE formtb(pb,km,g)
    !
    ! This subroutine assembles an unsymmetrical band matrix pb from
    ! element constituent matrices km.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::km(:,:)
     INTEGER,INTENT(IN)::g(:)
     REAL(iwp),INTENT(OUT)::pb(:,:)
     INTEGER::i,j,idof,icd,iw
     idof=SIZE(km,1)
     iw=(SIZE(pb,2)-1)/2
     DO i=1,idof
       IF(g(i)/=0)THEN
         DO j=1,idof
           IF(g(j)/=0)THEN
             icd=g(j)-g(i)+iw+1
             pb(g(i),icd)=pb(g(i),icd)+km(i,j)
           END IF
         END DO
       END IF
     END DO
    RETURN
  END SUBROUTINE formtb
  
  SUBROUTINE paraview(input,realisation,argv,nff,nelsy,g_coord,g_num,nf,nels,nod,nn,nlen,  &
    diag,ddylds,d1x1,d2x1,gravlo,loads,fcont,mv,kdiag,vcm,residual,damp,KVC,penpos,bod)

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::g_coord(:,:),mv(:)
    REAL(iwp),INTENT(IN)::diag(:),loads(:),d1x1(:),d2x1(:),fcont(:),vcm(:),residual(:),damp(:),KVC(:)
    REAL(iwp),INTENT(IN)::ddylds(:),gravlo(:)
    INTEGER,INTENT(IN)::input,nels,nod,nn,nlen,g_num(:,:),nf(:,:),realisation,kdiag(:),nff,nelsy,penpos(:),bod
    CHARACTER(*),INTENT(IN)::argv
    INTEGER::i,iel,ss,m,n
    REAL(iwp):: zero=0.0_iwp,dis_vec(2,nn),x
    character (len=8) ::  cnumber,cnumber1
    ss=input!+100
    write(cnumber,'(i8.6)') ss
    write(cnumber1,'(i8.6)') realisation
    ss=15
    !open(ss,FILE = argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
    IF(bod==1)THEN
    OPEN(ss,FILE="Output/Paraview2_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
    ELSE IF(bod==2)THEN
    OPEN(ss,FILE="Output/Paraview_2DM_1/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
    ELSE
    OPEN(ss,FILE="Output/Paraview_2DM_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
    END IF

    WRITE(ss,'(a)')'# vtk DataFile Version 3.0'
    WRITE(ss,'(a)')"vtk output"
    WRITE(ss,'(a)')"ASCII"
    WRITE(ss,'(a)')""
    WRITE(ss,'(a)')"DATASET UNSTRUCTURED_GRID"

    WRITE(ss,'(1A6,1I5,1A7)')"POINTS", nn , "float"
    DO i=1, nn
    WRITE(ss,'(3f9.4)') g_coord(:,i), zero
    END DO

    WRITE(ss,'(a)')""
    SELECT CASE(nod)
    CASE(4)
    WRITE(ss,'(1A6,2I6)')"CELLS", nels, nels*(1+4)
    DO iel=1, nels
    WRITE(ss,'(5I5)')nod,  &
    g_num(1,iel)-1, g_num(4,iel)-1, g_num(3,iel)-1, g_num(2,iel)-1
    END DO

    WRITE(ss,'(a)')""
    WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
    DO iel = 1 , nels
    WRITE(ss,*)"9"
    END DO

    CASE(8)
    WRITE(ss,'(1A5,2I5)')"CELLS", nels, nels*(1+8)
    DO iel=1, nels
    WRITE(ss,'(9I5)')nod,  &
    g_num(3,iel)-1, g_num(5,iel)-1, g_num(7,iel)-1, g_num(1,iel)-1, &
                g_num(4,iel)-1, g_num(6,iel)-1, g_num(8,iel)-1, g_num(2,iel)-1
    END DO

    WRITE(ss,'(a)')""
    WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
    DO iel = 1 , nels
    WRITE(ss,*)"23"
    END DO

    CASE(9)
    WRITE(ss,'(1A5,2I5)')"CELLS", (nels), (nels)*(1+9)
    DO iel=1, nels
    WRITE(ss,'(10I5)')nod,  &
    g_num(3,iel)-1, g_num(5,iel)-1, g_num(7,iel)-1, g_num(1,iel)-1, &
                g_num(4,iel)-1, g_num(6,iel)-1, g_num(8,iel)-1, g_num(2,iel)-1, &
    g_num(9,iel)-1 
    END DO

    WRITE(ss,'(a)')""
    WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
    DO iel = 1 , nels
    WRITE(ss,*)"28"
    END DO

    CASE DEFAULT
    WRITE(*,*)"wrong number of nodes input in paraview"
    END SELECT

    WRITE(ss,'(a)')""
    !
    WRITE(ss,'(1A10,1I5)')"POINT_DATA", nn

    WRITE(ss,'(a)')"vectors Ftotal float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')loads(nf(1:2,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors Fint_s float "
    DO  i=1, nn
    WRITE(ss,'(3f15.6)')ddylds(nf(1:2,i)+1), zero
    END DO


    WRITE(ss,'(a)')"vectors Fext_s float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')gravlo(nf(1:2,i)+1), zero
    END DO

    !  
    WRITE(ss,'(a)')"vectors fkin_s float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')vcm(nf(1:2,i)+1), zero
    END DO

        WRITE(ss,'(a)')"vectors f_damp float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')damp(nf(1:2,i)+1), zero
    END DO


    WRITE(ss,'(a)')"vectors R_s float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')residual(nf(1:2,i)+1), zero
    END DO


    !WRITE(ss,'(a)')"vectors normals float "
    ! DO i=1, nn
    !   WRITE(ss,'(3f15.6)')normals(1,i),normals(2,i), zero
    !  !WRITE(ss,'(3f9.4)')normals(nf(:,i)+1), zero
    ! END DO 

    WRITE(ss,'(a)')"vectors Fcont float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')fcont(nf(1:2,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors Vel_s float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')d1x1(nf(1:2,i)+1), zero
    END DO


    WRITE(ss,'(a)')"vectors Acc_s float "
    DO i=1, nn
    WRITE(ss,'(3f15.6)')d2x1(nf(1:2,i)+1), zero
    END DO 


    WRITE(ss,'(a)')"vectors Mv_s float "
    DO i=1, nn
    IF(nf(1,i)==0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')zero,zero, zero
    IF(nf(1,i)>0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')mv(kdiag(nf(1,i))),zero, zero
    IF(nf(1,i)==0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')zero,mv(kdiag(nf(2,i))), zero
    IF(nf(1,i)>0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')mv(kdiag(nf(1,i))),mv(kdiag(nf(2,i))), zero
    END DO

    !WRITE(ss,'(a)')"vectors K_coupled float "
    !DO i=1, nn
    !   WRITE(ss,'(3f15.4)')KVC(penpos(i)),KVC(penpos(i)), zero
    !END DO



    close(ss)
    RETURN
    END SUBROUTINE paraview
    
  SUBROUTINE point_viz(input,realisation,argv,gm_coord,m_stress,m_eps_acum,evpt,a_ins,Devstress,   &
    mstress,mpyield,cohesion,m_velocity,acc,nmps,nlen,bod)
    !---- SUBROUTINE used to save visualise outputs to Paraview format ---------------------

    IMPLICIT NONE  ! sr,pore,estress, removed from subroutine since not used now
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::gm_coord(:,:),m_stress(:,:),a_ins(:,:),evpt(:),cohesion(:),m_velocity(:,:),   &
    Devstress(:),acc(:,:),m_eps_acum(:,:),mstress(:),mpyield(:)
    INTEGER,INTENT(IN)::input,nmps,nlen,bod
    integer,intent(in)::realisation
    CHARACTER(*),INTENT(IN)::argv
    INTEGER::i,j,ss
    REAL(iwp):: zero=0.0_iwp
    character (len=8) ::  cnumber,cnumber1
    ! fflag : flag of a point yielded or nor
    ! realisation: used for monte-carlo simulation
    ! adgama: incremental shear strain
    ! acstrain: accumulated shear strain
    ! sr: saturation ratio
    ! pore: pore pressure
    ! estress: elastic stress?!
    ss=input
    !ss=input/20
    !/100000
    write(cnumber,'(i8.6)') ss
    write(cnumber1,'(i8.6)') realisation
    !open(ss,FILE=argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp')
    ss=15
    IF(bod==1)THEN
    OPEN(ss,FILE="Output/Paraview2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber
    ELSE IF(bod==2)THEN
    OPEN(ss,FILE="Output/Paraview_2DP_1/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber   
    ELSE
    OPEN(ss,FILE="Output/Paraview_2DP_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber   
    END IF

    write(ss,*) '#'
    write(ss,*) '# Simple AVS UCD File'
    write(ss,*) '#'
    write(ss,'(2i7.1,a)') nmps, nmps, ' 3 0 0' ! Number of nodes, number of cells, number of data items per node, number of data items per cell, model number (always 0)

    do j=1,nmps ! coordinates
    write(ss,'(i5.1,3g13.5)') j, gm_coord(:,j), zero
    ! MP number, x-coordinate y coordinate zero
    end do

    do j=1,nmps ! points
    write(ss,'(2i5.1,a,i5.1)') j, 0, '1 pt ', j
    end do

    write(ss,*) '10 2 1 4 4 2 2 1 1 1 1'              !1 1 1  for parameters not included (see below)
    ! 1st nr = number of variables you want to display, afterwards - how many components for each
    write(ss,*) 'displacement, meters'
    write(ss,*) 'Void, meters'
    write(ss,*) 'Effective-Stress, kpa'
    write(ss,*) 'Total-strains, no-unit'
    write(ss,*) 'Acc, kpa'
    write(ss,*) 'velocity, m/s'
    write(ss,*) 'Young, kpa'
    !write(ss,*) 'pore, kpa'
    !write(ss,*) 'strain invariant, no-unit'
    !write(ss,*) 'sr, no-unit'
    write(ss,*) 'Dev-stress, no-unit'
    write(ss,*) 'Mean-stress, no-unit'
    write(ss,*) 'Pore-pressure, no-unit'

    ! 1 yield; 0 nor yield
    do j=1,nmps ! properties
    write(ss,'(i5.1,20f15.4)') j, a_ins(:,j),evpt(j),m_stress(:,j),m_eps_acum(:,j),acc(:,j), &
            m_velocity(:,j),cohesion(j),Devstress(j),mstress(j),mpyield(j)
    end do

    CLOSE(ss)
    RETURN
    END SUBROUTINE point_viz  
  
  SUBROUTINE get_ff_displacement(mpm_disp, mpm_counter, mpm_bc_elements, mpm_g_num, mpm_nf, ff_disp, ff_g_num, ff_nf)
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(OUT)::mpm_disp(0:)
    INTEGER,INTENT(OUT)::mpm_counter(0:)
    REAL(iwp),INTENT(IN)::ff_disp(0:)
    INTEGER,INTENT(IN)::mpm_bc_elements(:),mpm_g_num(:,:),ff_g_num(:,:),mpm_nf(:,:),ff_nf(:,:)
    
    INTEGER::i, iel_bc_mpm, iel_bc_ff
    INTEGER,ALLOCATABLE::g_mp(:), g_ff(:)
    mpm_disp = 0.0_iwp
    mpm_counter = 0
    DO i=1,size(mpm_bc_elements)
      ! get boundary element index in ff & mpm
      iel_bc_mpm = mpm_bc_elements(i)
      iel_bc_ff  = i*2             ! 2 column for each row
    
      ! get g in mpm & ff from element number
      g_mp = mpm_nf(1,mpm_g_num(:,iel_bc_mpm))
      g_ff = ff_nf(1,ff_g_num(:,iel_bc_ff))
    
      ! take displacement with from ff to mpm and add counter
      mpm_disp(g_mp) = mpm_disp(g_mp) + ff_disp(g_ff)
      mpm_counter(g_mp) = mpm_counter(g_mp) + 1
    END DO
    ! cleanup memory
    mpm_disp(0)=0.0_iwp; mpm_counter(0)=0
    
    ! normalize displacement
    DO i=0, size(mpm_disp)-1
      IF(mpm_counter(i) > 0) mpm_disp(i) = mpm_disp(i) / mpm_counter(i)
    END DO
    mpm_disp(0)=0.0_iwp; mpm_counter(0)=0
  END SUBROUTINE
END MODULE fem