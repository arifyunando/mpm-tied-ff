module mpm
    save

contains

SUBROUTINE point_viz2(input,realisation,argv,gm_coord,m_stress,m_eps_acum,evpt,a_ins,Devstress,   &
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
END SUBROUTINE point_viz2  


SUBROUTINE floc(coord,gm_coord,mpoints,i)

    ! local coordinate of the material points

    USE MAIN
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:)
    INTEGER,INTENT(IN)::i
    REAL(iwp),INTENT(IN OUT)::mpoints(:,:)
    REAL(iwp),ALLOCATABLE::guess_local(:,:),fun(:),guess_global(:,:), &
    der(:,:),jac(:,:)
    INTEGER::count,nod,ndim,nmps
    REAL(iwp)::tolerance1,tolerance2,beta,zero=0.0_iwp,pt9=0.9_iwp,one=1.0_iwp, &
            tolerance
    LOGICAL:: find_local
    nod=UBOUND(coord,1)
    ndim=UBOUND(coord,2)
    nmps=UBOUND(gm_coord,2)
    ! --------------------
    ALLOCATE(fun(nod),guess_local(1,ndim),guess_global(1,ndim), &
    der(ndim,nod),jac(ndim,ndim))
    find_local=.FALSE.
    guess_local=zero
    CALL shape_fun(fun,guess_local,1)
    guess_global(1,:)=MATMUL(fun,coord)
    count=0
    DO
    count=count+1
    beta=one
    IF(count>=1) beta=pt9
    guess_global(1,:)=gm_coord(:,1)-guess_global(1,:)
    CALL shape_der(der,guess_local,1)
    jac=MATMUL(der,coord)
    CALL invert(jac)
!guess_local(i,:)=guess_local(i,:)+beta*MATMUL(guess_global(i,:),jac)
    guess_local(1,:)=guess_local(1,:)+beta*MATMUL(TRANSPOSE(jac),guess_global(1,:))
    CALL shape_fun(fun,guess_local,1)
    guess_global(1,:)=MATMUL(fun,coord)
    tolerance1=ABS(MAXVAL(guess_global(1,:)-gm_coord(:,1)))
    tolerance2=ABS(MINVAL(guess_global(1,:)-gm_coord(:,1)))
    !find_local=(tolerance1<=1.0D-6)
    !find_local=(tolerance2<=1.0D-6)
    !IF(find_local) THEN
    IF(tolerance1<=1.0D-10.and.tolerance2<=1.0D-10)THEN
        mpoints(1,:)=guess_local(1,:)
        EXIT
    END IF
    IF(count>10) THEN
        write(*,*)" wrong numerical method"
        PAUSE
        EXIT
    END IF
    END DO
    DEALLOCATE(fun,guess_local,guess_global,der,jac)
    RETURN
END SUBROUTINE floc
    
    
SUBROUTINE paraview2(input,realisation,argv,nff,nelsy,g_coord,g_num,nf,nels,nod,nn,nlen,  &
    diag,ddylds,d1x1,d2x1,gravlo,loads,fcont,mv,kdiag,vcm,residual,damp,KVC,nels_2D,penpos,bod)

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::g_coord(:,:),mv(:)
    REAL(iwp),INTENT(IN)::diag(:),loads(:),d1x1(:),d2x1(:),fcont(:),vcm(:),residual(:),damp(:),KVC(:)
    REAL(iwp),INTENT(IN)::ddylds(:),gravlo(:)
    INTEGER,INTENT(IN)::input,nels,nod,nn,nlen,g_num(:,:),nf(:,:),realisation,kdiag(:),nff,nelsy,nels_2D,penpos(:),bod
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
END SUBROUTINE paraview2
    
    
    
end module mpm