module mpm_mpm
    save

contains


SUBROUTINE alphavalue2(cohesion,stress_int,stress_ext,alpha_1)
    USE MPM_MAIN
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::cohesion
    REAL(iwp),INTENT(OUT)::alpha_1
    REAL(iwp),INTENT(IN)::stress_int(:),stress_ext(:)
    REAL(iwp)::sigm,dsbar,lode_theta,fa,fb,d3=3.0
 
    CALL invar(stress_int,sigm,dsbar,lode_theta)
    fa=dsbar-SQRT(d3)*cohesion
    CALL invar(stress_ext,sigm,dsbar,lode_theta)
    fb=dsbar-SQRT(d3)*cohesion
    IF(fa>0.0)THEN
        alpha_1=0.0
    ELSE
        alpha_1=((-1.0)*fa)/(fb-fa)
    END IF
END SUBROUTINE alphavalue2


SUBROUTINE beemat3(bee,deriv,values,SFUNCT)
    !
    ! This subroutine forms the bee matrix for 16 and 9 nodes.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     INTEGER,INTENT(IN)::values
     REAL(iwp),INTENT(IN)::SFUNCT(:),deriv(:,:)
     REAL(iwp),INTENT(OUT)::bee(:,:)
     INTEGER::l,m,n,i,nod,k

     bee=0.0_iwp

    IF(values==16)THEN
      bee(1,1)=deriv(1,1)
      bee(3,1)=deriv(2,1)
      bee(2,2)=deriv(2,1)
      bee(3,2)=deriv(1,1)
      bee(1,3)=deriv(1,2)
      bee(3,3)=deriv(2,2)
      bee(2,4)=deriv(2,2)
      bee(3,4)=deriv(1,2)
      bee(1,5)=deriv(1,3)
      bee(3,5)=deriv(2,3)
      bee(2,6)=deriv(2,3)
      bee(3,6)=deriv(1,3)
      bee(1,7)=deriv(1,4)
      bee(3,7)=deriv(2,4)
      bee(2,8)=deriv(2,4)
      bee(3,8)=deriv(1,4)
      bee(1,9)=deriv(1,5)
      bee(3,9)=deriv(2,5)
      bee(2,10)=deriv(2,5)
      bee(3,10)=deriv(1,5)
      bee(1,11)=deriv(1,6)
      bee(3,11)=deriv(2,6)
      bee(2,12)=deriv(2,6)
      bee(3,12)=deriv(1,6)
      bee(1,13)=deriv(1,7)
      bee(3,13)=deriv(2,7)
      bee(2,14)=deriv(2,7)
      bee(3,14)=deriv(1,7)
      bee(1,15)=deriv(1,8)
      bee(3,15)=deriv(2,8)
      bee(2,16)=deriv(2,8)
      bee(3,16)=deriv(1,8)
      bee(1,17)=deriv(1,9)
      bee(3,17)=deriv(2,9)
      bee(2,18)=deriv(2,9)
      bee(3,18)=deriv(1,9)
      bee(1,19)=deriv(1,10)
      bee(3,19)=deriv(2,10)
      bee(2,20)=deriv(2,10)
      bee(3,20)=deriv(1,10)
      bee(1,21)=deriv(1,11)
      bee(3,21)=deriv(2,11)
      bee(2,22)=deriv(2,11)
      bee(3,22)=deriv(1,11)
      bee(1,23)=deriv(1,12)
      bee(3,23)=deriv(2,12)
      bee(2,24)=deriv(2,12)
      bee(3,24)=deriv(1,12)
      bee(1,25)=deriv(1,13)
      bee(3,25)=deriv(2,13)
      bee(2,26)=deriv(2,13)
      bee(3,26)=deriv(1,13)
      bee(1,27)=deriv(1,14)
      bee(3,27)=deriv(2,14)
      bee(2,28)=deriv(2,14)
      bee(3,28)=deriv(1,14)
      bee(1,29)=deriv(1,15)
      bee(3,29)=deriv(2,15)
      bee(2,30)=deriv(2,15)
      bee(3,30)=deriv(1,15)
      bee(1,31)=deriv(1,16)
      bee(3,31)=deriv(2,16)
      bee(2,32)=deriv(2,16)
      bee(3,32)=deriv(1,16)
      !bee(4,1:values*2-1:2)=SFUNCT(:)/gc(1)
      bee(4,1:values*2-1:2)=0.0

     ELSE
      bee(1,1)=deriv(1,1)
      bee(3,1)=deriv(2,1)
      bee(2,2)=deriv(2,1)
      bee(3,2)=deriv(1,1)
      bee(1,3)=deriv(1,2)
      bee(3,3)=deriv(2,2)
      bee(2,4)=deriv(2,2)
      bee(3,4)=deriv(1,2)
      bee(1,5)=deriv(1,3)
      bee(3,5)=deriv(2,3)
      bee(2,6)=deriv(2,3)
      bee(3,6)=deriv(1,3)
      bee(1,7)=deriv(1,4)
      bee(3,7)=deriv(2,4)
      bee(2,8)=deriv(2,4)
      bee(3,8)=deriv(1,4)
      bee(1,9)=deriv(1,5)
      bee(3,9)=deriv(2,5)
      bee(2,10)=deriv(2,5)
      bee(3,10)=deriv(1,5)
      bee(1,11)=deriv(1,6)
      bee(3,11)=deriv(2,6)
      bee(2,12)=deriv(2,6)
      bee(3,12)=deriv(1,6)
      bee(1,13)=deriv(1,7)
      bee(3,13)=deriv(2,7)
      bee(2,14)=deriv(2,7)
      bee(3,14)=deriv(1,7)
      bee(1,15)=deriv(1,8)
      bee(3,15)=deriv(2,8)
      bee(2,16)=deriv(2,8)
      bee(3,16)=deriv(1,8)
      bee(1,17)=deriv(1,9)
      bee(3,17)=deriv(2,9)
      bee(2,18)=deriv(2,9)
      bee(3,18)=deriv(1,9)
      !bee(4,1:values*2-1:2)=SFUNCT(:)/gc(1)
      bee(4,1:values*2-1:2)=0.0
      !pause
     END IF

    RETURN
END SUBROUTINE beemat3


SUBROUTINE bound_detect(iel,neighb,bound,values,c_ele,nf,g_num,aver)
    !
    !   This subroutine produces derivatives of shape functions withe respect
    !   to local coordinates for 16 nodes.
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     INTEGER,INTENT(OUT)::bound,values,aver
     INTEGER,INTENT(IN)::neighb(:,:),iel,c_ele(:),g_num(:,:),nf(:,:)
     INTEGER::i,CONT
     REAL,ALLOCATABLE::NEIG(:)
     ALLOCATE(NEIG(8))
     CONT=0
     NEIG=0
     bound=0

    !Next loop is to know if there are neighbour elements and if those elements have MPs inside
    DO i=1,8
        IF(neighb(iel,i)==0)THEN
            NEIG(i)=0
        ELSE
            IF(c_ele(neighb(iel,i))>=1)NEIG(i)=1
            IF(c_ele(neighb(iel,i))==0)NEIG(i)=0
        END IF
    END DO

    IF(NEIG(1)==0.and.NEIG(2)==0.and.NEIG(3)==0.and.NEIG(4)>=1                  &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)==0.and.NEIG(8)==0)THEN
        bound=1;aver=1 !upper left corner
    END IF
    IF(NEIG(1)==0.and.NEIG(2)==0.and.NEIG(3)==0.and.NEIG(4)>=1                  &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=2;aver=2 !upper boundarie
    END IF
    IF(NEIG(1)==0.and.NEIG(2)==0.and.NEIG(3)==0.and.NEIG(4)==0                  &
        .and.NEIG(5)==0.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=3;aver=1 !upper right
    END IF
    IF(NEIG(1)==0.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1                  &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)==0.and.NEIG(8)==0)THEN
        bound=4;aver=2 !left boundarie
    END IF
    IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1                  &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=5;aver=1 !center element
    END IF
    IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)==0.and.NEIG(4)==0                  &
        .and.NEIG(5)==0.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=6;aver=2 !right boundarie
    END IF
    IF(NEIG(1)==0.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1                  &
        .and.NEIG(5)==0.and.NEIG(6)==0.and.NEIG(7)==0.and.NEIG(8)==0)THEN
        bound=7;aver=1 !lower left
    END IF
    IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1                  &
        .and.NEIG(5)==0.and.NEIG(6)==0.and.NEIG(7)==0.and.NEIG(8)>=1)THEN
        bound=8;aver=2 !lower boundarie
    END IF
    IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)==0.and.NEIG(4)==0                  &
        .and.NEIG(5)==0.and.NEIG(6)==0.and.NEIG(7)==0.and.NEIG(8)>=1)THEN
        bound=9;aver=1 !lower right
    END IF
    IF(bound>=1.and.bound/=5)THEN
        values=9
    ELSE IF(bound==5)THEN
        values=16
    ELSE
        values=4
    END IF
    RETURN
END SUBROUTINE bound_detect


SUBROUTINE couma(nels,nmps,a_ele,c_ele,k_ele,etype)

    ! count the number of material points inside an element.

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::a_ele(:)
    INTEGER,INTENT(IN)::nels,nmps
    INTEGER,INTENT(OUT)::c_ele(:),k_ele(0:),etype(:)
    INTEGER::i,iel,count
    k_ele=0
    DO iel=1,nels
    count=0
    DO i=1,nmps
        IF(iel==a_ele(i)) THEN
            count=count+1
        END IF
    END DO
    c_ele(iel)=count
    k_ele(iel)=k_ele(iel-1)+c_ele(iel)
    IF(count==0) THEN
        etype(iel)=1
    ELSE
        etype(iel)=2
    END IF
    END DO
    RETURN
END SUBROUTINE couma
    

subroutine formdlambda(flowf,flowg,dee,deps,a,dlambda)
    USE MPM_MAIN
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::dee(:,:),deps(:),a,flowf(:),flowg(:)
    REAL(iwp),INTENT(OUT)::dlambda
    REAL(iwp)::top,bot
    INTEGER::nst
    REAL(iwp),ALLOCATABLE::caflow(:),dsigmae(:)
     nst=UBOUND(deps,1)
    ALLOCATE(caflow(nst),dsigmae(nst))
    caflow=MATMUL(dee,flowg)
    bot=dot_product(flowf,caflow)
    bot=bot+a
    dsigmae=matmul(dee,deps)
    top=dot_product(flowf,dsigmae)
    dlambda=top/bot
    deallocate(caflow,dsigmae)
    RETURN
end subroutine formdlambda


SUBROUTINE inv(matrix,values)
    !
    ! This subroutine inverts the values x values matrix
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN OUT)::matrix(:,:)
    INTEGER,INTENT(IN)::values
    INTEGER,ALLOCATABLE::ipvt(:)
    REAL(iwp)::c,d,zero=0.0_iwp
    REAL(iwp),ALLOCATABLE::INVERSE(:,:),temp(:)
    INTEGER::l,k,N,M,j,imax(1)

    ALLOCATE(INVERSE(values,values),ipvt(values),temp(values))
    N=UBOUND(matrix,1)
    INVERSE=zero
    temp=zero
    ipvt=zero
    INVERSE = matrix
    ipvt = (/ (l, l = 1, N) /)
    DO k = 1,N
        imax = MAXLOC(ABS(INVERSE(k:n,k)))
        m = k-1+imax(1)
        IF (M /= k) THEN
            ipvt( (/M,k/) ) = ipvt( (/k,M/) )
            INVERSE((/M,k/),:) = INVERSE((/k,M/),:)
        END IF
        d = 1/INVERSE(k,k)
        temp = INVERSE(:,k)
        DO j = 1, N
            c = INVERSE(k,j)*d
            INVERSE(:,j) = INVERSE(:,j)-temp*c
            INVERSE(k,j) = c
        END DO
        INVERSE(:,k) = temp*(-d)
        INVERSE(k,k) = d
    END DO
    matrix(:,ipvt) = INVERSE
    RETURN
END SUBROUTINE inv


SUBROUTINE floc(coord,gm_coord,mpoints,i)
    !
    ! local coordinate of the material points
    !
    USE MPM_MAIN
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:)
    INTEGER,INTENT(IN)::i
    REAL(iwp),INTENT(IN OUT)::mpoints(:,:)
    REAL(iwp),ALLOCATABLE::guess_local(:,:),fun(:),guess_global(:,:),           &
        der(:,:),jac(:,:)
    INTEGER::count,nod,ndim,nmps
    REAL(iwp)::tolerance1,tolerance2,beta,zero=0.0_iwp,pt9=0.9_iwp,one=1.0_iwp, &
        tolerance
    LOGICAL:: find_local
    nod=UBOUND(coord,1)
    ndim=UBOUND(coord,2)
    nmps=UBOUND(gm_coord,2)
    ! --------------------
    ALLOCATE(fun(nod),guess_local(1,ndim),guess_global(1,ndim),                 &
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
        guess_local(1,:)=guess_local(1,:)+beta*MATMUL(TRANSPOSE(jac),guess_global(1,:))
        CALL shape_fun(fun,guess_local,1)
        guess_global(1,:)=MATMUL(fun,coord)
        tolerance1=ABS(MAXVAL(guess_global(1,:)-gm_coord(:,1)))
        tolerance2=ABS(MINVAL(guess_global(1,:)-gm_coord(:,1)))

        IF(tolerance1<=1.0D-10.and.tolerance2<=1.0D-10)THEN
            mpoints(1,:)=guess_local(1,:)
            EXIT
        END IF
        IF(count>10) THEN
            write(*,*)" wrong numerical method"
            EXIT
        END IF
    END DO
    DEALLOCATE(fun,guess_local,guess_global,der,jac)
    RETURN
END SUBROUTINE floc


SUBROUTINE point_viz2(input,realisation,argv,gm_coord,m_stress,m_stress_ini,evpt,a_ins,Devstress,   &
    mstress,mpyield,cohesion,m_velocity,acc,nmps,nlen)
    !---- SUBROUTINE used to save visualise outputs to Paraview format ---------------------

    IMPLICIT NONE  ! sr,pore,estress, removed from subroutine since not used now
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::gm_coord(:,:),m_stress(:,:),a_ins(:,:),evpt(:),cohesion(:),m_velocity(:,:),   &
    Devstress(:),acc(:,:),m_stress_ini(:,:),mstress(:),mpyield(:)
    INTEGER,INTENT(IN)::input,nmps,nlen
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
    write(cnumber,'(i8.6)') ss
    write(cnumber1,'(i8.6)') realisation
    ss=15
    OPEN(ss,FILE="Output/Paraview2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber
    write(ss,*) '#'
    write(ss,*) '# Simple AVS UCD File'
    write(ss,*) '#'
    write(ss,'(2i7.1,a)') nmps, nmps, ' 3 0 0' ! Number of nodes, number of cells, number of data items per node, number of data items per cell, model number (always 0)

    do j=1,nmps ! coordinates
        write(ss,'(i5.1,3g13.5)') j, gm_coord(:,j), zero
    end do

    do j=1,nmps ! points
        write(ss,'(2i5.1,a,i5.1)') j, 0, '1 pt ', j
    end do

    write(ss,*) '10 2 1 4 4 2 2 1 1 1 1' !1 1 1  for parameters not included (see below)
    write(ss,*) 'displacement, meters'
    write(ss,*) 'plastic strains, meters'
    write(ss,*) 'Stress, kpa'
    write(ss,*) 'Stress increase, kpa'
    write(ss,*) 'Acc, kpa'
    write(ss,*) 'velocity, m/s'
    write(ss,*) 'cohesion, kpa'
    write(ss,*) 'Dev-stress, no-unit'
    write(ss,*) 'Mean-stress, no-unit'
    write(ss,*) 'Yield-value, no-unit'
    do j=1,nmps ! properties
        write(ss,'(i5.1,19f13.3)') j, a_ins(:,j),evpt(j),m_stress(:,j),m_stress_ini(:,j),acc(:,j), &
                m_velocity(:,j),cohesion(j),Devstress(j),mstress(j),mpyield(j)
    end do

    CLOSE(ss)
    RETURN
END SUBROUTINE point_viz2


SUBROUTINE neigh_b(nnxe,nnye,nels,neighb,dir)

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::nnxe,nnye,nels
    INTEGER,INTENT(OUT)::neighb(:,:)
    CHARACTER(LEN=1),INTENT(IN)::dir
    INTEGER::row,col,i,neighbor

    row=1
    col=1
    IF(dir=='y')THEN
        DO i=1, nels
            IF(row==1.and.col==1)THEN
                neighb(i,1)=0
                neighb(i,2)=0
                neighb(i,3)=0
                neighb(i,7)=0
                neighb(i,8)=0
                neighb(i,4)=nnye*col+row
                neighb(i,5)=nnye*col+row+1
                neighb(i,6)=nnye*(col-1)+row+1
            ELSE IF(col==1.and.row>1.and.row<nnye)THEN
                neighb(i,1)=0
                neighb(i,7)=0
                neighb(i,8)=0
                neighb(i,2)=nnye*(col-1)+row-1
                neighb(i,3)=nnye*col+row-1
                neighb(i,4)=nnye*col+row
                neighb(i,5)=nnye*col+row+1
                neighb(i,6)=nnye*(col-1)+row+1
            ELSE IF(row==nnye.and.col==1)THEN
                neighb(i,1)=0
                neighb(i,5)=0
                neighb(i,6)=0
                neighb(i,7)=0
                neighb(i,8)=0
                neighb(i,2)=nnye*(col-1)+row-1
                neighb(i,3)=nnye*col+row-1
                neighb(i,4)=nnye*col+row
            ELSE IF(row==nnye.and.col>1.and.col<nnxe)THEN
                neighb(i,5)=0
                neighb(i,6)=0
                neighb(i,7)=0
                neighb(i,1)=nnye*(col-2)+row-1
                neighb(i,2)=nnye*(col-1)+row-1
                neighb(i,3)=nnye*col+row-1
                neighb(i,4)=nnye*col+row
                neighb(i,8)=nnye*(col-2)+row
            ELSE IF(row==nnye.and.col==nnxe)THEN
                neighb(i,3)=0
                neighb(i,4)=0
                neighb(i,5)=0
                neighb(i,6)=0
                neighb(i,7)=0
                neighb(i,1)=nnye*(col-2)+row-1
                neighb(i,2)=nnye*(col-1)+row-1
                neighb(i,8)=nnye*(col-2)+row
            ELSE IF(col==nnxe.and.row>1.and.row<nnye)THEN
                neighb(i,3)=0
                neighb(i,4)=0
                neighb(i,5)=0
                neighb(i,1)=nnye*(col-2)+row-1
                neighb(i,2)=nnye*(col-1)+row-1
                neighb(i,6)=nnye*(col-1)+row+1
                neighb(i,7)=nnye*(col-2)+row+1
                neighb(i,8)=nnye*(col-2)+row
            ELSE IF(row==1.and.col==nnxe)THEN
                neighb(i,1)=0
                neighb(i,2)=0
                neighb(i,3)=0
                neighb(i,4)=0
                neighb(i,5)=0
                neighb(i,6)=nnye*(col-1)+row+1
                neighb(i,7)=nnye*(col-2)+row+1
                neighb(i,8)=nnye*(col-2)+row
            ELSE IF(row==1.and.col>1.and.col<nnxe)THEN
                neighb(i,1)=0
                neighb(i,2)=0
                neighb(i,3)=0
                neighb(i,4)=nnye*col+row
                neighb(i,5)=nnye*col+row+1
                neighb(i,6)=nnye*(col-1)+row+1
                neighb(i,7)=nnye*(col-2)+row+1
                neighb(i,8)=nnye*(col-2)+row
            ELSE IF(row>1.and.row<nnye.and.col>1.and.col<nnxe)THEN
                neighb(i,1)=nnye*(col-2)+row-1
                neighb(i,2)=nnye*(col-1)+row-1
                neighb(i,3)=nnye*col+row-1
                neighb(i,4)=nnye*col+row
                neighb(i,5)=nnye*col+row+1
                neighb(i,6)=nnye*(col-1)+row+1
                neighb(i,7)=nnye*(col-2)+row+1
                neighb(i,8)=nnye*(col-2)+row
            END IF
            IF(row==nnye)THEN
                row=0
                col=col+1
            END IF
            row=row+1
        END DO
    ELSE
        DO i=1, nels
            IF(row==1.and.col==1)THEN
                neighb(i,1)=0
                neighb(i,2)=0
                neighb(i,3)=0
                neighb(i,7)=0
                neighb(i,8)=0
                neighb(i,4)=col+1
                neighb(i,5)=nnxe*row+col+1
                neighb(i,6)=nnxe*row+1
            ELSE IF(col==1.and.row>1.and.row<nnye)THEN
                neighb(i,1)=0
                neighb(i,7)=0
                neighb(i,8)=0
                neighb(i,2)=nnxe*(row-2)+1
                neighb(i,3)=nnxe*(row-2)+2
                neighb(i,4)=nnxe*(row-1)+2
                neighb(i,5)=nnxe*row+2
                neighb(i,6)=nnxe*row+1
            ELSE IF(row==nnye.and.col==1)THEN
                neighb(i,1)=0
                neighb(i,5)=0
                neighb(i,6)=0
                neighb(i,7)=0
                neighb(i,8)=0
                neighb(i,2)=nnxe*(row-2)+1
                neighb(i,3)=nnxe*(row-2)+2
                neighb(i,4)=nnxe*(row-1)+2
            ELSE IF(row==nnye.and.col>1.and.col<nnxe)THEN
                neighb(i,5)=0
                neighb(i,6)=0
                neighb(i,7)=0
                neighb(i,1)=nnxe*(row-2)+col-1
                neighb(i,2)=nnxe*(row-2)+col
                neighb(i,3)=nnxe*(row-2)+col+1
                neighb(i,4)=nnxe*(row-1)+col+1
                neighb(i,8)=nnxe*(row-1)+col-1
            ELSE IF(row==nnye.and.col==nnxe)THEN
                neighb(i,3)=0
                neighb(i,4)=0
                neighb(i,5)=0
                neighb(i,6)=0
                neighb(i,7)=0
                neighb(i,1)=nnxe*(row-2)+col-1
                neighb(i,2)=nnxe*(row-2)+col
                neighb(i,8)=nnxe*(row-1)+col-1
            ELSE IF(col==nnxe.and.row>1.and.row<nnye)THEN
                neighb(i,3)=0
                neighb(i,4)=0
                neighb(i,5)=0
                neighb(i,1)=nnxe*(row-2)+col-1
                neighb(i,2)=nnxe*(row-2)+col
                neighb(i,6)=nnxe*row+col
                neighb(i,7)=nnxe*row+col-1
                neighb(i,8)=nnxe*(row-1)+col-1
            ELSE IF(row==1.and.col==nnxe)THEN
                neighb(i,1)=0
                neighb(i,2)=0
                neighb(i,3)=0
                neighb(i,4)=0
                neighb(i,5)=0
                neighb(i,6)=nnxe*row+col
                neighb(i,7)=nnxe*row+col-1
                neighb(i,8)=nnxe*(row-1)+col-1
            ELSE IF(row==1.and.col>1.and.col<nnxe)THEN
                neighb(i,1)=0
                neighb(i,2)=0
                neighb(i,3)=0
                neighb(i,4)=col+1
                neighb(i,5)=nnxe*row+col+1
                neighb(i,6)=nnxe*row+col
                neighb(i,7)=nnxe*row+col-1
                neighb(i,8)=col-1
            ELSE IF(row>1.and.row<nnye.and.col>1.and.col<nnxe)THEN
                neighb(i,1)=nnxe*(row-2)+col-1
                neighb(i,2)=nnxe*(row-2)+col
                neighb(i,3)=nnxe*(row-2)+col+1
                neighb(i,4)=nnxe*(row-1)+col+1
                neighb(i,5)=nnxe*(row)+col+1
                neighb(i,6)=nnxe*(row)+col
                neighb(i,7)=nnxe*(row)+col-1
                neighb(i,8)=nnxe*(row-1)+col-1
            END IF

            IF(col==nnxe)THEN
                col=0
                row=row+1
            END IF
            col=col+1
        END DO
    END IF
    RETURN
END SUBROUTINE neigh_b


SUBROUTINE emb_2d_bc2(nx1,nx2,ny1,newel,nf)
    !
    ! This subroutine generates the nf array for a 2-d slope geometry.
    !
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nx1,nx2,ny1,newel
    INTEGER,INTENT(OUT)::nf(:,:)
    INTEGER::nm,ic,i,j,k,nye,val,nxe

    nxe=nx1+newel
    nye=ny1
    nm=1
    k=1
    val=1
    j=1
    nf=0
    DO i=1,nxe+1
        j=1
        IF(i>nx1+2)nye=nye-1
        DO WHILE (j<=nye+1)
            IF(i==1.and.j<nye+1)THEN
                nf(1,k)=0
                nf(2,k)=val
                val=val+1
            ELSE IF(i==1.and.j==nye+1)THEN
                nf(1,k)=0
                nf(2,k)=0
            ELSE IF(i>1.and.j<nye+1)THEN
                nf(1,k)=val
                val=val+1
                nf(2,k)=val
                val=val+1
            END IF
            k=k+1
            j=j+1
        END DO
    END DO
    RETURN
END SUBROUTINE emb_2d_bc2


SUBROUTINE shapefunct(points,i,CONSTANTS,SFUNCT)
    !
    !   This subroutine produces derivatives of shape functions withe respect
    !   to local coordinates for 16 and 9 nodes.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::i
    REAL(iwp)::eta,xi
    REAL(iwp),INTENT(IN)::points(:,:)
    REAL(iwp),INTENT(IN)::CONSTANTS(:,:)
    REAL(iwp),INTENT(OUT)::SFUNCT(:)
    REAL,PARAMETER::zero=0.0_iwp
    INTEGER::ndim,n
    xi=points(i,1)
    eta=points(i,2)
    ndim=UBOUND(CONSTANTS,1)

    IF(ndim==16)THEN
        DO n=1,ndim
        SFUNCT(n)=CONSTANTS(1,n)+CONSTANTS(2,n)*xi+eta*CONSTANTS(3,n)+xi*eta*CONSTANTS(4,n)+ &
                    CONSTANTS(5,n)*xi**2+CONSTANTS(6,n)*eta**2+CONSTANTS(7,n)*xi**2*eta+ &
                    CONSTANTS(8,n)*xi*eta**2+CONSTANTS(9,n)*xi**2*eta**2+CONSTANTS(10,n)*xi**3+ &
                    CONSTANTS(11,n)*eta**3+CONSTANTS(12,n)*eta*xi**3+CONSTANTS(13,n)*eta**3*xi+ &
                    CONSTANTS(14,n)*xi**3*eta**2+CONSTANTS(15,n)*xi**2*eta**3+CONSTANTS(16,n)*xi**3*eta**3
        END DO
    ELSE
        DO n=1,ndim
        SFUNCT(n)=CONSTANTS(1,n)+CONSTANTS(2,n)*xi+eta*CONSTANTS(3,n)+xi*eta*CONSTANTS(4,n)+ &
                    CONSTANTS(5,n)*xi**2+CONSTANTS(6,n)*eta**2+CONSTANTS(7,n)*xi**2*eta+ &
                    CONSTANTS(8,n)*xi*eta**2+CONSTANTS(9,n)*xi**2*eta**2
        END DO
    END IF
    RETURN
END SUBROUTINE shapefunct
    

SUBROUTINE yieldcorrection(stress,k,ci,gamap1,cp,gamap2,cr,hp,sp,dee,ft)
    USE MPM_MAIN
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::hp,sp,ci,gamap1,cp,gamap2,cr,dee(:,:)
    REAL(iwp),INTENT(IN OUT)::stress(:),k,ft
    REAL(iwp)::bot,dnamda,sigm,dsbar,lode_theta,dq1,dq2,dq3,a,f_trial,f0,ftol=1.0d-06, &
    zero=0.0d0,d3=3.0d0,two=2.0d0,k0,k_trial,h,cf,theta
    INTEGER::nst,iters,maxits=10
    REAL(iwp),ALLOCATABLE::caflow(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),flowf(:),flowg(:), &
    stress0(:),stress_trial(:),vmfl(:)
    nst=UBOUND(stress,1)
    ALLOCATE(caflow(nst),m1(nst,nst),m2(nst,nst),m3(nst,nst),flowf(nst),flowg(nst),flow(nst,nst), &
    stress0(nst),stress_trial(nst),vmfl(nst))
    stress0=stress
    k0=k
    iters=0
    !CALL formh(k0,gamap1,gamap2,hp,sp,h)
    CALL plasticmod(k0,gamap2,sp,h)
    a=sqrt(d3)*h
    iters_2: DO
        iters=iters+1
        CALL invar(stress0,sigm,dsbar,lode_theta)
        cf=k0*sp+cp
        IF(cf<cr) cf=cr
        f0=dsbar-sqrt(d3)*cf
        CALL vmflow(stress0,dsbar,vmfl)
        flowf=vmfl
        flowg= flowf
        CALL plasticmod(k0,gamap2,sp,h)
        a=sqrt(d3)*h
        caflow=MATMUL(dee,flowg)
        bot=dot_product(flowf,caflow)
        dnamda=f0/(a+bot)
        stress_trial=stress0-dnamda*caflow
        k_trial=k0+dnamda
        CALL invar(stress_trial,sigm,dsbar,lode_theta)
        cf=k_trial*sp+cp
        IF(cf<cr)cf=cr
        f_trial=dsbar-sqrt(d3)*cf   ! cf=k_trial ;;;;;;  my assumption not sure...
        IF(abs(f_trial)>abs(f0))THEN
            bot=dot_product(flowf,flowf)
            dnamda=f0/bot
            stress_trial=stress0-dnamda*flowf
            k_trial=k0
            CALL invar(stress_trial,sigm,dsbar,lode_theta)
            cf=k_trial*sp+cp
            IF(cf<cr) cf=cr
            f_trial=dsbar-sqrt(d3)*cf
        END IF
        IF(abs(f_trial)<=ftol) EXIT
        stress0=stress_trial
        k0=k_trial
        IF(iters==maxits) write(*,*) "correction stress not convergent"
    END DO iters_2

    stress=stress_trial
    k=k_trial
    CALL invar(stress,sigm,dsbar,lode_theta)
    cf=k*sp+cp
    IF(cf<cr)cf=cr
    ft=dsbar-sqrt(d3)*cf
    DEALLOCATE(caflow,m1,m2,m3,flowf,flowg,flow,stress0,stress_trial)
    RETURN
END SUBROUTINE yieldcorrection

end module mpm_mpm