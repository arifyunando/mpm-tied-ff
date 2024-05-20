MODULE FUNCTIONS
  IMPLICIT NONE
  ! Note: Naming Convention for dynamic 2D array saved in fixed array structure
  !   variable_v = actual value of variable
  !   variable_n = actual length of each index of the variable
  !   variable_i = If variable_v is saved as 1D-Array, variable_i indicates
  !                the index where each member array starts
  !   Example:
  !
  !     [[a, b, c, d, e],
  !      [f, g, h],
  !      [i, j, k, l]]
  !     
  !      variable_v:
  !        2D fixed array:
  !          [[a, b, c, d, e]
  !           [f, g, h, 0, 0]
  !           [i, j, k, l, 0]]
  !
  !        1D fixed array:
  !          [a, b, c, d, e, f, g, h, i, j, k, l, 0, 0, 0, 0, 0, 0]
  !
  !      variable_n:
  !          [5, 3, 4]
  !     
  !      variable_i
  !          [1, 6, 9]
  CONTAINS
  
!************************* COMMON FEM & MPM FUNCTIONS *************************!
!                                                                              !

  FUNCTION DETERMINANT(jac)RESULT(det)
    !
    ! This function returns the determinant of a 1x1, 2x2 or 3x3
    ! Jacobian matrix.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::jac(:,:)
    REAL(iwp)::det
    INTEGER::it
    it=UBOUND(jac,1)
    SELECT CASE(it) !-Matrix Size
    CASE(1)
      det=1.0_iwp
    CASE(2)
      det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
    CASE(3)
      det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
      det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
      det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
    CASE DEFAULT
      WRITE(*,*)' wrong dimension for Jacobian matrix'
    END SELECT
    RETURN
  END FUNCTION DETERMINANT


  SUBROUTINE BEEMAT(bee,deriv)
    !
    ! This subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6).
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::deriv(:,:)
    REAL(iwp),INTENT(OUT)::bee(:,:)
    INTEGER::k,l,m,n,ih,nod
    REAL::x,y,z
    bee=0.0_iwp
    ih=UBOUND(bee,1)
    nod=UBOUND(deriv,2)
    SELECT CASE (ih) ! number of stress/strain terms
    CASE(3,4) ! plane strain/stress and axisymmetric
      DO m=1,nod
        k=2*m
        l=k-1
        x=deriv(1,m)
        y=deriv(2,m)
        bee(1,l)=x
        bee(3,k)=x
        bee(2,k)=y
        bee(3,l)=y
      END DO
    CASE(6) ! 3-Dimensional Conditions
      DO m=1,nod
        n=3*m
        k=n-1
        l=k-1
        x=deriv(1,m)
        y=deriv(2,m)
        z=deriv(3,m)
        bee(1,l)=x
        bee(4,k)=x
        bee(6,n)=x
        bee(2,k)=y
        bee(4,l)=y
        bee(5,n)=y
        bee(3,n)=z
        bee(5,k)=z
        bee(6,l)=z
      END DO
    CASE DEFAULT
      WRITE(*,*) 'wrong dimension for nst in bee matrix'
    END SELECT
    RETURN
  END SUBROUTINE BEEMAT


  SUBROUTINE CHECK_CONVERGENCE(converged,DuC,DuP,tol)
    !
    ! This subroutine sets converged to .FALSE. if relative change in loads
    ! and oldlds is greater than tol and updates oldlds.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::DuC(0:),DuP(0:)
    REAL(iwp),INTENT(IN), OPTIONAL::tol
    LOGICAL,INTENT(OUT)::converged
    INTEGER::i
    REAL(iwp)::VAL,upvalB1,upvalB2,downval,energy_acum,energy_step,work
    ! Default Arguments
    REAL(iwp)::d_tol = 1.0e-8_iwp
    IF(present(tol)) d_tol = tol
    
    converged=.FALSE.
    upvalB1=0.0;upvalB2=0.0;downval=0.0;energy_acum=0.0;energy_step=0.0

    DO i=1,SIZE(DuC)-1
      upvalB1=upvalB1+DuC(i)**2
    END DO    
    
    upvalB2=MAXVAL(abs(DuP))
    upvalB2=upvalB2**2
    
    upvalB1=SQRT(upvalB1)
    upvalB2=SQRT(upvalB2)
    
    VAL=(ABS(upvalB1))/(ABS(upvalB2))
    CONVERGED=(VAL<=d_tol)
    RETURN
  END SUBROUTINE CHECK_CONVERGENCE


  SUBROUTINE DEEMAT(dee,e,v)
    !
    ! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
    ! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
    ! (three dimensions).
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::e,v
    REAL(iwp),INTENT(OUT)::dee(:,:)
    REAL(iwp)::v1,v2,c,vv,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp
    INTEGER::i,ih
    dee=zero
    ih=UBOUND(dee,1)
    v1=one-v
    c=e/((one+v)*(one-two*v))
    SELECT CASE(ih)
    CASE(3)
      dee(1,1)=v1*c
      dee(2,2)=v1*c
      dee(1,2)=v*c
      dee(2,1)=v*c
      dee(3,3)=pt5*c*(one-two*v)
    CASE(4)
      dee(1,1)=v1*c
      dee(2,2)=v1*c
      dee(4,4)=v1*c
      dee(3,3)=pt5*c*(one-two*v)
      dee(1,2)=v*c
      dee(2,1)=v*c
      dee(1,4)=v*c
      dee(4,1)=v*c
      dee(2,4)=v*c
      dee(4,2)=v*c
    CASE(6)
      v2=v/(one-v)
      vv=(one-two*v)/(one-v)*pt5
      DO i=1,3
        dee(i,i)=one
      END DO
      DO i=4,6
        dee(i,i)=vv
      END DO
      dee(1,2)=v2
      dee(2,1)=v2
      dee(1,3)=v2
      dee(3,1)=v2
      dee(2,3)=v2
      dee(3,2)=v2
      dee=dee*e/(two*(one+v)*vv)
    CASE DEFAULT
      WRITE(*,*)'wrong size for dee matrix'
    END SELECT
    RETURN
  END SUBROUTINE DEEMAT


  SUBROUTINE DERIVATIVES(points,i,CONSTANTS,SIXTEEND)
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
    REAL(iwp),INTENT(OUT)::SIXTEEND(:,:)
    REAL,PARAMETER::zero=0.0_iwp
    INTEGER::m,ndim,n
    xi=points(i,1)
    eta=points(i,2)
    ndim=UBOUND(CONSTANTS,1)
    IF(ndim==16)THEN
      DO n=1,ndim
        DO m=1,2
          IF(m==1) SIXTEEND(m,n)=CONSTANTS(2,n)+eta*CONSTANTS(4,n)+2*xi*CONSTANTS(5,n)+ &
            2*CONSTANTS(7,n)*xi*eta+CONSTANTS(8,n)*eta**2+2*CONSTANTS(9,n)*xi*eta**2+   &
            3*CONSTANTS(10,n)*xi**2+3*CONSTANTS(12,n)*xi**2*eta+CONSTANTS(13,n)*eta**3+ &
            3*CONSTANTS(14,n)*xi**2*eta**2+2*CONSTANTS(15,n)*eta**3*xi+                 &
            3*CONSTANTS(16,n)*xi**2*eta**3
          IF(m==2)SIXTEEND(m,n)=CONSTANTS(3,n)+xi*CONSTANTS(4,n)+2*eta*CONSTANTS(6,n)+  &
            CONSTANTS(7,n)*xi**2+2*CONSTANTS(8,n)*eta*xi+2*CONSTANTS(9,n)*xi**2*eta+    &
            3*CONSTANTS(11,n)*eta**2+CONSTANTS(12,n)*xi**3+3*CONSTANTS(13,n)*eta**2*xi+ &
            2*CONSTANTS(14,n)*xi**3*eta+3*CONSTANTS(15,n)*eta**2*xi**2+                 &
            3*CONSTANTS(16,n)*xi**3*eta**2
        END DO
      END DO
    END IF
    IF(ndim==9)THEN
      DO n=1,ndim
        DO m=1,2
          IF(m==1)SIXTEEND(m,n)=CONSTANTS(2,n)+eta*CONSTANTS(4,n)+2*xi*CONSTANTS(5,n)+ &
            2*CONSTANTS(7,n)*xi*eta+CONSTANTS(8,n)*eta**2+2*CONSTANTS(9,n)*xi*eta**2
          IF(m==2)SIXTEEND(m,n)=CONSTANTS(3,n)+xi*CONSTANTS(4,n)+2*eta*CONSTANTS(6,n)+ &
            CONSTANTS(7,n)*xi**2+2*CONSTANTS(8,n)*eta*xi+2*CONSTANTS(9,n)*xi**2*eta
        END DO
      END DO
    END IF
    RETURN
  END SUBROUTINE DERIVATIVES


  SUBROUTINE FORMNF(nf)
    !
    ! This subroutine forms the nf matrix.
    !
    IMPLICIT NONE
    INTEGER,INTENT(IN OUT)::nf(:,:)
    INTEGER::i,j,m
    m=0
    DO j=1,UBOUND(nf,2)
      DO i=1,UBOUND(nf,1)
        IF(nf(i,j)/=0)THEN
          m=m+1
          nf(i,j)=m
        END IF
      END DO
    END DO
    RETURN
  END SUBROUTINE FORMNF


  SUBROUTINE NUM_TO_G(num,nf,g)
    !
    ! This subroutine finds the g vector from num and nf.
    !
    IMPLICIT NONE
    INTEGER,INTENT(IN)::num(:),nf(:,:)
    INTEGER,INTENT(OUT)::g(:)
    INTEGER::i,k,nod,nodof
    nod=UBOUND(num,1)
    nodof=UBOUND(nf,1)
    DO i=1,nod
      k=i*nodof
      g(k-nodof+1:k)=nf(:,num(i))
    END DO
    RETURN
  END SUBROUTINE NUM_TO_G


  SUBROUTINE SAMPLE_GAUSS_IP(element,s,wt)
    !
    ! This subroutine returns the local coordinates and weighting coefficients
    ! of the integrating points.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(OUT)::s(:,:)
    REAL(iwp),INTENT(OUT),OPTIONAL::wt(:)
    CHARACTER(*),INTENT(IN)::element
    INTEGER::nip
    REAL(iwp)::root3,r15,w(3),v(9),b,c
    root3=1.0_iwp/SQRT(3.0_iwp)
    r15=0.2_iwp*SQRT(15.0_iwp)
    nip=UBOUND(s,1)
    w=(/5.0_iwp/9.0_iwp,8.0_iwp/9.0_iwp,5.0_iwp/9.0_iwp/)
    v=(/5.0_iwp/9.0_iwp*w,8.0_iwp/9.0_iwp*w,5.0_iwp/9.0_iwp*w/)
    SELECT CASE(element)
    CASE('line')
      SELECT CASE(nip) ! number of integrating points
      CASE(1)
        s(1,1)=0.0_iwp
        wt(1) =2.0_iwp
      CASE(2)
        s(1,1)=-0.577350269189626_iwp
        s(2,1)= 0.577350269189626_iwp
        wt(1) = 1.000000000000000_iwp
        wt(2) = 1.000000000000000_iwp
      CASE(3)
        s(1,1)=-0.774596669241484_iwp
        s(2,1)= 0.000000000000000_iwp
        s(3,1)= 0.774596669241484_iwp
        wt(1) = 0.555555555555556_iwp
        wt(2) = 0.888888888888889_iwp
        wt(3) = 0.555555555555556_iwp
      CASE(4)
        s(1,1)=-0.861136311594053_iwp
        s(2,1)=-0.339981043584856_iwp
        s(3,1)= 0.339981043584856_iwp
        s(4,1)= 0.861136311594053_iwp
        wt(1) = 0.347854845137454_iwp
        wt(2) = 0.652145154862546_iwp
        wt(3) = 0.652145154862546_iwp
        wt(4) = 0.347854845137454_iwp
      CASE(5)
        s(1,1)=-0.906179845938664_iwp
        s(2,1)=-0.538469310105683_iwp
        s(3,1)= 0.000000000000000_iwp
        s(4,1)= 0.538469310105683_iwp
        s(5,1)= 0.906179845938664_iwp
        wt(1) = 0.236926885056189_iwp
        wt(2) = 0.478628670499366_iwp
        wt(3) = 0.568888888888889_iwp
        wt(4) = 0.478628670499366_iwp
        wt(5) = 0.236926885056189_iwp
      CASE(6)
        s(1,1)=-0.932469514203152_iwp
        s(2,1)=-0.661209386466265_iwp
        s(3,1)=-0.238619186083197_iwp
        s(4,1)= 0.238619186083197_iwp
        s(5,1)= 0.661209386466265_iwp
        s(6,1)= 0.932469514203152_iwp
        wt(1) = 0.171324492379170_iwp
        wt(2) = 0.360761573048139_iwp
        wt(3) = 0.467913934572691_iwp
        wt(4) = 0.467913934572691_iwp
        wt(5) = 0.360761573048139_iwp
        wt(6) = 0.171324492379170_iwp
      CASE(7)
        s(1,1)=-0.9491079123427585245261897_iwp
        s(2,1)=-0.7415311855993944398638648_iwp
        s(3,1)=-0.4058451513773971669066064_iwp
        s(4,1)= 0.000000000000000_iwp
        s(5,1)= 0.4058451513773971669066064_iwp
        s(6,1)= 0.7415311855993944398638648_iwp
        s(7,1)= 0.9491079123427585245261897_iwp
        wt(1) = 0.1294849661688696932706114_iwp
        wt(2) = 0.2797053914892766679014678_iwp
        wt(3) = 0.3818300505051189449503698_iwp
        wt(4) = 0.4179591836734693877551020_iwp
        wt(5) = 0.3818300505051189449503698_iwp
        wt(6) = 0.2797053914892766679014678_iwp
        wt(7) = 0.1294849661688696932706114_iwp
      CASE(8)
        s(1,1)=-0.9602898564975362316835609_iwp
        s(2,1)=-0.7966664774136267395915539_iwp
        s(3,1)=-0.5255324099163289858177390_iwp
        s(4,1)=-0.1834346424956498049394761_iwp
        s(5,1)= 0.1834346424956498049394761_iwp
        s(6,1)= 0.5255324099163289858177390_iwp
        s(7,1)= 0.7966664774136267395915539_iwp
        s(8,1)= 0.9602898564975362316835609_iwp
        wt(1) = 0.1012285362903762591525314_iwp
        wt(2) = 0.2223810344533744705443560_iwp
        wt(3) = 0.3137066458778872873379622_iwp
        wt(4) = 0.3626837833783619829651504_iwp
        wt(5) = 0.3626837833783619829651504_iwp
        wt(6) = 0.3137066458778872873379622_iwp
        wt(7) = 0.2223810344533744705443560_iwp
        wt(8) = 0.1012285362903762591525314_iwp
      CASE(9)
        s(1,1)=-0.9681602395076260898355762_iwp
        s(2,1)=-0.8360311073266357942994298_iwp
        s(3,1)=-0.6133714327005903973087020_iwp
        s(4,1)=-0.3242534234038089290385380_iwp
        s(5,1)= 0.000000000000000_iwp
        s(6,1)= 0.3242534234038089290385380_iwp
        s(7,1)= 0.6133714327005903973087020_iwp
        s(8,1)= 0.8360311073266357942994298_iwp
        s(9,1)= 0.9681602395076260898355762_iwp
        wt(1) = 0.0812743883615744119718922_iwp
        wt(2) = 0.1806481606948574040584720_iwp
        wt(3) = 0.2606106964029354623187429_iwp
        wt(4) = 0.3123470770400028400686304_iwp
        wt(5) = 0.3302393550012597631645251_iwp
        wt(6) = 0.3123470770400028400686304_iwp
        wt(7) = 0.2606106964029354623187429_iwp
        wt(8) = 0.1806481606948574040584720_iwp
        wt(9) = 0.0812743883615744119718922_iwp
      CASE(10)
        s(1,1)=-0.9739065285171717200779640_iwp
        s(2,1)=-0.8650633666889845107320967_iwp
        s(3,1)=-0.6794095682990244062343274_iwp
        s(4,1)=-0.4333953941292471907992659_iwp
        s(5,1)=-0.1488743389816312108848260_iwp
        s(6,1)= 0.1488743389816312108848260_iwp
        s(7,1)= 0.4333953941292471907992659_iwp
        s(8,1)= 0.6794095682990244062343274_iwp
        s(9,1)= 0.8650633666889845107320967_iwp
        s(10,1)= 0.9739065285171717200779640_iwp
        wt(1) = 0.0666713443086881375935688_iwp
        wt(2) = 0.1494513491505805931457763_iwp
        wt(3) = 0.2190863625159820439955349_iwp
        wt(4) = 0.2692667193099963550912269_iwp
        wt(5) = 0.2955242247147528701738930_iwp
        wt(6) = 0.2955242247147528701738930_iwp
        wt(7) = 0.2692667193099963550912269_iwp
        wt(8) = 0.2190863625159820439955349_iwp
        wt(9) = 0.1494513491505805931457763_iwp
        wt(10) = 0.0666713443086881375935688_iwp
      CASE DEFAULT
        WRITE(*,*)"Wrong number of integrating points for a line"
      END SELECT
    CASE('triangle')
      SELECT CASE(nip)
      CASE(1)
        s(1,1)= 0.333333333333333_iwp
        s(1,2)= 0.333333333333333_iwp
        wt(1) = 0.500000000000000_iwp
      CASE(3)
        s(1,1)= 0.500000000000000_iwp
        s(1,2)= 0.500000000000000_iwp
        s(2,1)= 0.500000000000000_iwp
        s(2,2)= 0.000000000000000_iwp
        s(3,1)= 0.000000000000000_iwp
        s(3,2)= 0.500000000000000_iwp
        wt(1:3)=0.333333333333333_iwp
        wt=0.5_iwp*wt
      CASE(4)
        s(1,1)= 0.6_iwp
        s(1,2)= 0.2_iwp
        s(2,1)= 0.2_iwp
        s(2,2)= 0.6_iwp
        s(3,1)= 0.2_iwp
        s(3,2)= 0.2_iwp
        s(4,1)= 0.333333333333333_iwp
        s(4,2)= 0.333333333333333_iwp
        wt(1:3)= 0.520833333333333_iwp
        wt(4)=  -0.5625_iwp
        wt=0.5_iwp*wt
      CASE(6)
        s(1,1)= 0.816847572980459_iwp
        s(1,2)= 0.091576213509771_iwp
        s(2,1)= 0.091576213509771_iwp
        s(2,2)= 0.816847572980459_iwp
        s(3,1)= 0.091576213509771_iwp
        s(3,2)= 0.091576213509771_iwp
        s(4,1)= 0.108103018168070_iwp
        s(4,2)= 0.445948490915965_iwp
        s(5,1)= 0.445948490915965_iwp
        s(5,2)= 0.108103018168070_iwp
        s(6,1)= 0.445948490915965_iwp
        s(6,2)= 0.445948490915965_iwp
        wt(1:3)=0.109951743655322_iwp
        wt(4:6)=0.223381589678011_iwp
        wt=0.5_iwp*wt
      CASE(7)
        s(1,1)= 0.333333333333333_iwp
        s(1,2)= 0.333333333333333_iwp
        s(2,1)= 0.797426985353087_iwp
        s(2,2)= 0.101286507323456_iwp
        s(3,1)= 0.101286507323456_iwp
        s(3,2)= 0.797426985353087_iwp
        s(4,1)= 0.101286507323456_iwp
        s(4,2)= 0.101286507323456_iwp
        s(5,1)= 0.470142064105115_iwp
        s(5,2)= 0.059715871789770_iwp
        s(6,1)= 0.059715871789770_iwp
        s(6,2)= 0.470142064105115_iwp
        s(7,1)= 0.470142064105115_iwp
        s(7,2)= 0.470142064105115_iwp
        wt(1) = 0.225000000000000_iwp
        wt(2:4)=0.125939180544827_iwp
        wt(5:7)=0.132394152788506_iwp
        wt=0.5_iwp*wt
      CASE(12)
        s(1,1)= 0.873821971016996_iwp
        s(1,2)= 0.063089014491502_iwp
        s(2,1)= 0.063089014491502_iwp
        s(2,2)= 0.873821971016996_iwp
        s(3,1)= 0.063089014491502_iwp
        s(3,2)= 0.063089014491502_iwp
        s(4,1)= 0.501426509658179_iwp
        s(4,2)= 0.249286745170910_iwp
        s(5,1)= 0.249286745170910_iwp
        s(5,2)= 0.501426509658179_iwp
        s(6,1)= 0.249286745170910_iwp
        s(6,2)= 0.249286745170910_iwp
        s(7,1) =0.053145049844817_iwp
        s(7,2) =0.310352451033784_iwp
        s(8,1) =0.310352451033784_iwp
        s(8,2) =0.053145049844817_iwp
        s(9,1) =0.053145049844817_iwp
        s(9,2) =0.636502499121398_iwp
        s(10,1)=0.310352451033784_iwp
        s(10,2)=0.636502499121398_iwp
        s(11,1)=0.636502499121398_iwp
        s(11,2)=0.053145049844817_iwp
        s(12,1)=0.636502499121398_iwp
        s(12,2)=0.310352451033784_iwp
        wt(1:3)=0.050844906370207_iwp
        wt(4:6)=0.116786275726379_iwp
        wt(7:12)=0.082851075618374_iwp
        wt=0.5_iwp*wt
      CASE(16)
        s(1,1)=0.333333333333333_iwp
        s(1,2)=0.333333333333333_iwp
        s(2,1)=0.658861384496478_iwp
        s(2,2)=0.170569307751761_iwp
        s(3,1)=0.170569307751761_iwp
        s(3,2)=0.658861384496478_iwp
        s(4,1)=0.170569307751761_iwp
        s(4,2)=0.170569307751761_iwp
        s(5,1)=0.898905543365938_iwp
        s(5,2)=0.050547228317031_iwp
        s(6,1)=0.050547228317031_iwp
        s(6,2)=0.898905543365938_iwp
        s(7,1)=0.050547228317031_iwp
        s(7,2)=0.050547228317031_iwp
        s(8,1)=0.081414823414554_iwp
        s(8,2)=0.459292588292723_iwp
        s(9,1)=0.459292588292723_iwp
        s(9,2)=0.081414823414554_iwp
        s(10,1)=0.459292588292723_iwp
        s(10,2)=0.459292588292723_iwp
        s(11,1)=0.008394777409958_iwp
        s(11,2)=0.263112829634638_iwp
        s(12,1)=0.008394777409958_iwp
        s(12,2)=0.728492392955404_iwp
        s(13,1)=0.263112829634638_iwp
        s(13,2)=0.008394777409958_iwp
        s(14,1)=0.263112829634638_iwp
        s(14,2)=0.728492392955404_iwp
        s(15,1)=0.728492392955404_iwp
        s(15,2)=0.008394777409958_iwp
        s(16,1)=0.728492392955404_iwp
        s(16,2)=0.263112829634638_iwp
        wt(1)=0.144315607677787_iwp
        wt(2:4)=0.103217370534718_iwp
        wt(5:7)=0.032458497623198_iwp
        wt(8:10)=0.095091634267284_iwp
        wt(11:16)=0.027230314174435_iwp
        wt=0.5_iwp*wt
      CASE DEFAULT
        WRITE(*,*)"wrong number of integrating points for a triangle"
      END SELECT
    CASE('quadrilateral')
      SELECT CASE(nip)
      CASE(1)
        s(1,1)=0.0_iwp
        s(1,2)=0.0_iwp
        wt(1)=4.0_iwp
      CASE(4)
        s(1,1)=-root3
        s(1,2)= root3
        s(2,1)= root3
        s(2,2)= root3
        s(3,1)=-root3
        s(3,2)=-root3
        s(4,1)= root3
        s(4,2)=-root3
        wt=1.0_iwp
      CASE(9)
        s(1:7:3,1)=-r15
        s(2:8:3,1)=0.0_iwp
        s(3:9:3,1)=r15
        s(1:3,2)  =r15
        s(4:6,2)  =0.0_iwp
        s(7:9,2)  =-r15
        wt= v
      CASE(16)
        s(1:13:4,1)=-0.861136311594053_iwp
        s(2:14:4,1)=-0.339981043584856_iwp
        s(3:15:4,1)= 0.339981043584856_iwp
        s(4:16:4,1)= 0.861136311594053_iwp
        s(1:4,2)   = 0.861136311594053_iwp
        s(5:8,2)   = 0.339981043584856_iwp
        s(9:12,2)  =-0.339981043584856_iwp
        s(13:16,2) =-0.861136311594053_iwp
        wt(1)      = 0.121002993285602_iwp
        wt(4)      = wt(1)
        wt(13)     = wt(1)
        wt(16)     = wt(1)
        wt(2)      = 0.226851851851852_iwp
        wt(3)      = wt(2)
        wt(5)      = wt(2)
        wt(8)      = wt(2)
        wt(9)      = wt(2)
        wt(12)     = wt(2)
        wt(14)     = wt(2)
        wt(15)     = wt(2)
        wt(6)      = 0.425293303010694_iwp
        wt(7)      = wt(6)
        wt(10)     = wt(6)
        wt(11)     = wt(6)
      CASE(25)
        s(1:21:5,1)= 0.906179845938664_iwp
        s(2:22:5,1)= 0.538469310105683_iwp
        s(3:23:5,1)= 0.0_iwp
        s(4:24:5,1)=-0.538469310105683_iwp
        s(5:25:5,1)=-0.906179845938664_iwp
        s( 1: 5,2) = 0.906179845938664_iwp
        s( 6:10,2) = 0.538469310105683_iwp
        s(11:15,2) = 0.0_iwp
        s(16:20,2) =-0.538469310105683_iwp
        s(21:25,2) =-0.906179845938664_iwp
        wt(1) =0.056134348862429_iwp
        wt(2) =0.113400000000000_iwp
        wt(3) =0.134785072387521_iwp
        wt(4) =0.113400000000000_iwp
        wt(5) =0.056134348862429_iwp
        wt(6) =0.113400000000000_iwp
        wt(7) =0.229085404223991_iwp
        wt(8) =0.272286532550750_iwp
        wt(9) =0.229085404223991_iwp
        wt(10)=0.113400000000000_iwp
        wt(11)=0.134785072387521_iwp
        wt(12)=0.272286532550750_iwp
        wt(13)=0.323634567901235_iwp
        wt(14)=0.272286532550750_iwp
        wt(15)=0.134785072387521_iwp
        wt(16)=0.113400000000000_iwp
        wt(17)=0.229085404223991_iwp
        wt(18)=0.272286532550750_iwp
        wt(19)=0.229085404223991_iwp
        wt(20)=0.113400000000000_iwp
        wt(21)=0.056134348862429_iwp
        wt(22)=0.113400000000000_iwp
        wt(23)=0.134785072387521_iwp
        wt(24)=0.113400000000000_iwp
        wt(25)=0.056134348862429_iwp
      CASE DEFAULT
        WRITE(*,*)"wrong number of integrating points for a quadrilateral"
      END SELECT
    CASE('tetrahedron') ! for tetrahedra weights multiplied by 1/6
      SELECT CASE(nip)
      CASE(1)
        s(1,1)=0.25_iwp
        s(1,2)=0.25_iwp
        s(1,3)=0.25_iwp
        wt(1)=1.0_iwp/6.0_iwp
      CASE(4)
        s(1,1)=0.58541020_iwp
        s(1,2)=0.13819660_iwp
        s(1,3)=s(1,2)
        s(2,2)=s(1,1)
        s(2,3)=s(1,2)
        s(2,1)=s(1,2)
        s(3,3)=s(1,1)
        s(3,1)=s(1,2)
        s(3,2)=s(1,2)
        s(4,1)=s(1,2)
        s(4,2)=s(1,2)
        s(4,3)=s(1,2)
        wt(1:4)=0.25_iwp/6.0_iwp
      CASE(5)
        s(1,1)=0.25_iwp
        s(1,2)=0.25_iwp
        s(1,3)=0.25_iwp
        s(2,1)=0.5_iwp
        s(2,2)=1.0_iwp/6.0_iwp
        s(2,3)=s(2,2)
        s(3,2)=0.5_iwp
        s(3,3)=1.0_iwp/6.0_iwp
        s(3,1)=s(3,3)
        s(4,3)=0.5_iwp
        s(4,1)=1.0_iwp/6.0_iwp
        s(4,2)=s(4,1)
        s(5,1)=1.0_iwp/6.0_iwp
        s(5,2)=s(5,1)
        s(5,3)=s(5,1)
        wt(1)=-0.8_iwp
        wt(2)=9.0_iwp/20.0_iwp
        wt(3:5)=wt(2)
        wt=wt/6.0_iwp
      CASE DEFAULT
          WRITE(*,*)"wrong number of integrating points for a tetrahedron"
      END SELECT
    CASE('hexahedron')
      SELECT CASE(nip)
      CASE(1)
        s(1,1:3)=0.0_iwp
        wt(1)=8.0_iwp
      CASE(8)
        s(1,1)= root3
        s(1,2)= root3
        s(1,3)= root3
        s(2,1)= root3
        s(2,2)= root3
        s(2,3)=-root3
        s(3,1)= root3
        s(3,2)=-root3
        s(3,3)= root3
        s(4,1)= root3
        s(4,2)=-root3
        s(4,3)=-root3
        s(5,1)=-root3
        s(5,2)= root3
        s(5,3)= root3
        s(6,1)=-root3
        s(6,2)=-root3
        s(6,3)= root3
        s(7,1)=-root3
        s(7,2)= root3
        s(7,3)=-root3
        s(8,1)=-root3
        s(8,2)=-root3
        s(8,3)=-root3
        wt=1.0_iwp
      CASE(14)
        b=0.795822426_iwp
        c=0.758786911_iwp
        wt(1:6)=0.886426593_iwp
        wt(7:14)=0.335180055_iwp
        s(1,1)=-b
        s(2,1)=b
        s(3,2)=-b
        s(4,2)=b
        s(5,3)=-b
        s(6,3)=b
        s(7:,:)=c
        s(7,1)=-c
        s(7,2)=-c
        s(7,3)=-c
        s(8,2)=-c
        s(8,3)=-c
        s(9,1)=-c
        s(9,3)=-c
        s(10,3)=-c
        s(11,1)=-c
        s(11,2)=-c
        s(12,2)=-c
        s(13,1)=-c
      CASE(15)
        b=1.0_iwp
        c       =0.674199862_iwp
        wt(1)   =1.564444444_iwp
        wt(2:7) =0.355555556_iwp
        wt(8:15)=0.537777778_iwp
        s(2,1)=-b
        s(3,1)=b
        s(4,2)=-b
        s(5,2)=b
        s(6,3)=-b
        s(7,3)=b
        s(8:,:)=c
        s(8,1)=-c
        s(8,2)=-c
        s(8,3)=-c
        s(9,2)=-c
        s(9,3)=-c
        s(10,1)=-c
        s(10,3)=-c
        s(11,3)=-c
        s(12,1)=-c
        s(12,2)=-c
        s(13,2)=-c
        s(14,1)=-c
      CASE(27)
        wt=(/5.0_iwp/9.0_iwp*v,8.0_iwp/9.0_iwp*v,5.0_iwp/9.0_iwp*v/)
        s(1:7:3,1)=-r15
        s(2:8:3,1)=0.0_iwp
        s(3:9:3,1)=r15
        s(1:3,3)=r15
        s(4:6,3)=0.0_iwp
        s(7:9,3)=-r15
        s(1:9,2)=-r15
        s(10:16:3,1)=-r15
        s(11:17:3,1)=0.0_iwp
        s(12:18:3,1)=r15
        s(10:12,3)=r15
        s(13:15,3)=0.0_iwp
        s(16:18,3)=-r15
        s(10:18,2)=0.0_iwp
        s(19:25:3,1)=-r15
        s(20:26:3,1)=0.0_iwp
        s(21:27:3,1)=r15
        s(19:21,3)= r15
        s(22:24,3)=0.0_iwp
        s(25:27,3)=-r15
        s(19:27,2)= r15
      CASE DEFAULT
        WRITE(*,*)"wrong number of integrating points for a hexahedron"
      END SELECT
    CASE DEFAULT
      WRITE(*,*)"not a valid element type"
    END SELECT
    RETURN
  END SUBROUTINE SAMPLE_GAUSS_IP


  SUBROUTINE SHAPE_FUN(fun,points,i)
    !
    !   This subroutine computes the values of the shape functions.
    !   to local coordinates
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(in)::i
    REAL(iwp),INTENT(IN)::points(:,:)
    REAL(iwp),INTENT(OUT)::fun(:)
    REAL(iwp)::eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3
    REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9
    REAL(iwp)::zeta,xi0,eta0,zeta0
    INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
    REAL,PARAMETER::pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,pt75=0.75_iwp, &
      one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d8=8.0_iwp,d9=9.0_iwp,   &
      d16=16.0_iwp,d27=27.0_iwp,d32=32.0_iwp,d64=64.0_iwp,d128=128.0_iwp
    ndim=UBOUND(points,2)
    nod=UBOUND(fun,1)
    SELECT CASE(ndim)
    CASE(1) ! 1-Dimensional Case
      xi=points(i,1)
      SELECT CASE(nod)
      CASE(2)
        t1=-one-xi
        t2= one-xi
        fun(1)=t2/two
        fun(2)=-t1/two
      CASE(3)
        t1=-one-xi
        t2=-xi
        t3=one-xi
        fun(1)=t2*t3/two
        fun(2)=-t1*t3
        fun(3)=t1*t2/two
      CASE(4)
        t1=-one-xi
        t2=-one/d3-xi
        t3=one/d3-xi
        t4=one-xi
        fun(1)=t2*t3*t4*d9/d16
        fun(2)=-t1*t3*t4*d27/d16
        fun(3)=t1*t2*t4*d27/d16
        fun(4)=-t1*t2*t3*d9/d16
      CASE(5)
        t1=-one -xi
        t2=-pt5-xi
        t3=-xi
        t4=pt5-xi
        t5=one-xi
        fun(1)=t2*t3*t4*t5*two/d3
        fun(2)=-t1*t3*t4*t5*d8/d3
        fun(3)=t1*t2*t4*t5*d4
        fun(4)=-t1*t2*t3*t5*d8/d3
        fun(5)=t1*t2*t3*t4*two/d3
      CASE DEFAULT
        WRITE(*,*)"wrong number of nodes in shape_fun"
      END SELECT
    CASE(2) ! 2-Dimensional Case
      c1=points(i,1)
      c2=points(i,2)
      c3=one-c1-c2
      xi=points(i,1)
      eta=points(i,2)
      etam=pt25*(one-eta)
      etap=pt25*(one+eta)
      xim=pt25*(one-xi)
      xip=pt25*(one+xi)
      SELECT CASE(nod)
      CASE(3)
        fun = (/c1,c3,c2/)
      CASE(6)
        fun(1)=(two*c1-one)*c1
        fun(2)=d4*c3*c1
        fun(3)=(two*c3-one)*c3
        fun(4)=d4*c2*c3
        fun(5)=(two*c2-one)*c2
        fun(6)=d4*c1*c2
      CASE(10)
        fun(1)= ((d3*c1-one)*(d3*c1-two)*c1)/two
        fun(2)= -(d9*(d3*c1-one)*(c1+c2-one)*c1)/two
        fun(3)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c1)/two
        fun(4)=-((d3*c1+d3*c2-one)*(d3*c1+d3*c2-two)*(c1+c2-one))/two
        fun(5)=  (d9*(d3*c1+d3*c2-two)*(c1+c2-one)*c2)/two
        fun(6)= -(d9*(c1+c2-one)*(d3*c2-one)*c2)/two
        fun(7)= ((d3*c2-one)*(d3*c2-two)*c2)/two
        fun(8)=  (d9*(d3*c2-one)*c1*c2)/two
        fun(9)=  (d9*(d3*c1-one)*c1*c2)/two
        fun(10)=-d27*((c2-one)+c1)*c1*c2
      CASE(15)
        t1=c1-pt25
        t2=c1-pt5
        t3=c1-pt75
        t4=c2-pt25
        t5=c2-pt5
        t6=c2-pt75
        t7=c3-pt25
        t8=c3-pt5
        t9=c3-pt75
        fun(1)=d32/d3*c1*t1*t2*t3
        fun(2)=d128/d3*c3*c1*t1*t2
        fun(3)=d64*c3*c1*t1*t7
        fun(4)=d128/d3*c3*c1*t7*t8
        fun(5)=d32/d3*c3*t7*t8*t9
        fun(6)=d128/d3*c2*c3*t7*t8
        fun(7)=d64*c2*c3*t4*t7
        fun(8)=d128/d3*c2*c3*t4*t5
        fun(9)=d32/d3*c2*t4*t5*t6
        fun(10)=d128/d3*c1*c2*t4*t5
        fun(11)=d64*c1*c2*t1*t4
        fun(12)=d128/d3*c1*c2*t1*t2
        fun(13)=d128*c1*c2*t1*c3
        fun(15)=d128*c1*c2*c3*t4
        fun(14)=d128*c1*c2*c3*t7
      CASE(4)
        fun=(/d4*xim*etam,d4*xim*etap,d4*xip*etap,d4*xip*etam/)
      CASE(5)
        fun=(/d4*xim*etam-pt25*(one-xi**2)*(one-eta**2),             &
              d4*xim*etap-pt25*(one-xi**2)*(one-eta**2),             &
              d4*xip*etap-pt25*(one-xi**2)*(one-eta**2),             &
              d4*xip*etam-pt25*(one-xi**2)*(one-eta**2),             &
              (one-xi**2)*(one-eta**2)/)
      CASE(8)
        fun=(/d4*etam*xim*(-xi-eta-one),d32*etam*xim*etap,                   &
              d4*etap*xim*(-xi+eta-one),d32*xim*xip*etap,                    &
              d4*etap*xip*(xi+eta-one), d32*etap*xip*etam,                   &
              d4*xip*etam*(xi-eta-one), d32*xim*xip*etam/)
      CASE(9)
        etam=eta-one
        etap=eta+one
        xim=xi-one
        xip=xi+one
        fun=(/pt25*xi*xim*eta*etam,-pt5*xi*xim*etap*etam,                    &
              pt25*xi*xim*eta*etap,-pt5*xip*xim*eta*etap,                    &
              pt25*xi*xip*eta*etap,-pt5*xi*xip*etap*etam,                    &
              pt25*xi*xip*eta*etam,-pt5*xip*xim*eta*etam,                    &
              xip*xim*etap*etam/)
      CASE DEFAULT
        WRITE(*,*)"wrong number of nodes in shape_fun"
      END SELECT

    CASE(3) ! 3-Dimension Case
      xi=points(i,1)
      eta=points(i,2)
      zeta=points(i,3)
      etam=one-eta
      xim=one-xi
      zetam=one-zeta
      etap=eta+one
      xip=xi+one
      zetap=zeta+one
      SELECT CASE(nod)
      CASE(4)
        fun(1)=xi
        fun(2)=eta
        fun(3)=zeta
        fun(4)=one-fun(1)-fun(2)-fun(3)
      CASE(8)
        fun=(/pt125*xim*etam*zetam,pt125*xim*etam*zetap,                     &
              pt125*xip*etam*zetap,pt125*xip*etam*zetam,                     &
              pt125*xim*etap*zetam,pt125*xim*etap*zetap,                     &
              pt125*xip*etap*zetap,pt125*xip*etap*zetam/)
      CASE(14) ! type 6 element
        fun(1) = (xi*eta+xi*zeta+two*xi+eta*zeta+two*eta+two*zeta+two)*      &
          (xi-one)*(eta-one)*(zeta-one)/d8
        fun(2) =-(xi*eta-xi*zeta+two*xi-eta*zeta+two*eta-two*zeta+two)*      &
          (xi-one)*(eta-one)*(zeta+one)/d8
        fun(3) =-(xi*eta-xi*zeta+two*xi+eta*zeta-two*eta+two*zeta-two)*      &
          (xi+one)*(eta-one)*(zeta+one)/d8
        fun(4) = (xi*eta+xi*zeta+two*xi-eta*zeta-two*eta-two*zeta-two)*      &
          (xi+one)*(eta-one)*(zeta-one)/d8
        fun(5) =-(xi+one)*(xi-one)*(eta-one)*(zeta+one)*(zeta-one)/two
        fun(6) =-(xi-one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
        fun(7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta+one)/two
        fun(8) = (xi+one)*(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
        fun(9) =-(xi+one)*(xi-one)*(eta+one)*(eta-one)*(zeta-one)/two
        fun(10)= (xi*eta-xi*zeta-two*xi+eta*zeta+two*eta-two*zeta-two)*      &
          (xi-one)*(eta+one)*(zeta-one)/d8
        fun(11)=-(xi*eta+xi*zeta-two*xi-eta*zeta+two*eta+two*zeta-two)*      &
          (xi-one)*(eta+one)*(zeta+one)/d8
        fun(12)=-(xi*eta+xi*zeta-two*xi+eta*zeta-two*eta-two*zeta+two)*      &
          (xi+one)*(eta+one)*(zeta+one)/d8
        fun(13)= (xi*eta-xi*zeta-two*xi-eta*zeta-two*eta+two*zeta+two)*      &
          (xi+one)*(eta+one)*(zeta-one)/d8
        fun(14)= (xi+one)*(xi-one)*(eta+one)*(zeta+one)*(zeta-one)/two
      CASE(20)
        xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
        etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
        zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
        DO l=1,20
          xi0=xi*xii(l)
          eta0=eta*etai(l)
          zeta0=zeta*zetai(l)
          IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
            fun(l)=pt25*(one-xi*xi)*(one+eta0)*(one+zeta0)
          ELSE IF(l>=9.AND.l<=12)THEN
            fun(l)=pt25*(one+xi0)*(one-eta*eta)*(one+zeta0)
          ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18)THEN
            fun(l)=pt25*(one+xi0)*(one+eta0)*(one-zeta*zeta)
          ELSE
            fun(l)=pt125*(one+xi0)*(one+eta0)*(one+zeta0)*(xi0+eta0+zeta0-2)
          END IF
        END DO
      CASE DEFAULT
        WRITE(*,*)"wrong number of nodes in shape_fun"
      END SELECT
    CASE DEFAULT
      WRITE(*,*)"wrong number of dimensions in shape_fun"
    END SELECT
    RETURN
  END SUBROUTINE SHAPE_FUN


  SUBROUTINE SHAPE_DER(der,points,i)
    !
    !   This subroutine produces derivatives of shape functions withe respect
    !   to local coordinates.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::i
    REAL(iwp),INTENT(IN)::points(:,:)
    REAL(iwp),INTENT(OUT)::der(:,:)
    REAL(iwp)::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3
    REAL(iwp)::t1,t2,t3,t4,t5,t6,t7,t8,t9,x2p1,x2m1,e2p1,e2m1,zetam,zetap
    REAL,PARAMETER::zero=0.0_iwp,pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,     &
        pt75=0.75_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp, &
        d6=6.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d10=10.0_iwp,d11=11.0_iwp,             &
        d12=12.0_iwp,d16=16.0_iwp,d18=18.0_iwp,d27=27.0_iwp,d32=32.0_iwp,       &
        d36=36.0_iwp,d54=54.0_iwp,d64=64.0_iwp,d128=128.0_iwp
    INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
    ndim=UBOUND(der,1)
    nod= UBOUND(der,2)
    SELECT CASE(ndim)
    CASE(1)   ! one dimensional elements
    xi=points(i,1)
    SELECT CASE(nod)
        CASE(2)
            der(1,1)=-pt5
            der(1,2)= pt5
        CASE(3)
            t1=-one-xi
            t2=-xi
            t3=one-xi
            der(1,1)=-(t3+t2)/two
            der(1,2)=(t3+t1)
            der(1,3)=-(t2+t1)/two
        CASE(4)
            t1=-one-xi
            t2=-one/d3-xi
            t3=one/d3-xi
            t4=one-xi
            der(1,1)=-(t3*t4+t2*t4+t2*t3)*d9/d16
            der(1,2)=(t3*t4+t1*t4+t1*t3)*d27/d16
            der(1,3)=-(t2*t4+t1*t4+t1*t2)*d27/d16
            der(1,4)=(t2*t3+t1*t3+t1*t2)*d9/d16
        CASE(5)
            t1=-one-xi
            t2=-pt5-xi
            t3=-xi
            t4=pt5-xi
            t5=one-xi
            der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*two/d3
            der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*d8/d3
            der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*d4
            der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*d8/d3
            der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*two/d3
        CASE DEFAULT
            WRITE(*,*)"wrong number of nodes in shape_der"
        END SELECT
    CASE(2)      ! two dimensional elements
    xi=points(i,1)
    eta=points(i,2)
    c1=xi
    c2=eta
    c3=one-c1-c2
    etam=pt25*(one-eta)
    etap=pt25*(one+eta)
    xim= pt25*(one-xi)
    xip= pt25*(one+xi)
    x2p1=two*xi+one
    x2m1=two*xi-one
    e2p1=two*eta+one
    e2m1=two*eta-one
    SELECT CASE(nod)
        CASE(3)
            der(1,1)=one
            der(1,3)=zero
            der(1,2)=-one
            der(2,1)=zero
            der(2,3)=one
            der(2,2)=-one
        CASE(6)
            der(1,1)=d4*c1-one
            der(1,6)=d4*c2
            der(1,5)=zero
            der(1,4)=-d4*c2
            der(1,3)=-(d4*c3-one)
            der(1,2)=d4*(c3-c1)
            der(2,1)=zero
            der(2,6)=d4*c1
            der(2,5)=d4*c2-one
            der(2,4)=d4*(c3-c2)
            der(2,3)=-(d4*c3-one)
            der(2,2)=-d4*c1
        CASE(10)
            der(1,1)=(d27*c1**2-d18*c1+two)/two
            der(1,9)=(d9*(d6*c1-one)*c2)/two
            der(1,8)=(d9*(d3*c2-one)*c2)/two
            der(1,7)=zero
            der(1,6)=-(d9*(d3*c2-one)*c2)/two
            der(1,5)= (d9*(d6*c1+d6*c2-d5)*c2)/two
            der(1,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
            der(1,3)= (d9*(d9*c1**2+d12*c1*c2-d10*c1+d3*c2**2-d5*c2+two))/two
            der(1,2)=-(d9*(d9*c1**2+d6*c1*c2-d8*c1-c2+one))/two
            der(1,10)=-d27*(((c2-one)+c1)+c1)*c2
            der(2,1)=zero
            der(2,9)= (d9*(d3*c1-one)*c1)/two
            der(2,8)= (d9*(d6*c2-one)*c1)/two
            der(2,7)=(d27*c2**2-d18*c2+two)/two
            der(2,6)=-(d9*((c1+c2-one)*(d6*c2-one)+(d3*c2-one)*c2))/two
            der(2,5)= (d9*(d3*c1**2+d12*c1*c2-d5*c1+d9*c2**2-d10*c2+two))/two
            der(2,4)=-(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
            der(2,3)= (d9*(d6*c1+d6*c2-d5)*c1)/two
            der(2,2)=-(d9*(d3*c1-one)*c1)/two
            der(2,10)=-d27*(((c2-one)+c1)+c2)*c1
        CASE(15)
            t1=c1-pt25
            t2=c1-pt5
            t3=c1-pt75
            t4=c2-pt25
            t5=c2-pt5
            t6=c2-pt75
            t7=c3-pt25
            t8=c3-pt5
            t9=c3-pt75
            der(1,1)=d32/d3*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
            der(1,12)=d128/d3*c2*(t2*(t1+c1)+c1*t1)
            der(1,11)=d64*c2*t4*(t1+c1)
            der(1,10)=d128/d3*c2*t4*t5
            der(1,9)=zero
            der(1,8)=-d128/d3*c2*t4*t5
            der(1,7)=-d64*c2*t4*(t7+c3)
            der(1,6)=-d128/d3*c2*(t8*(t7+c3)+c3*t7)
            der(1,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
            der(1,4)=d128/d3*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
            der(1,3)=d64*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
            der(1,2)=d128/d3*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
            der(1,13)=d128*c2*(c3*(t1+c1)-c1*t1)
            der(1,15)=d128*c2*t4*(c3-c1)
            der(1,14)=d128*c2*(c3*t7-c1*(t7+c3))
            der(2,1)=zero
            der(2,12)=d128/d3*c1*t1*t2
            der(2,11)=d64*c1*t1*(t4+c2)
            der(2,10)=d128/d3*c1*(t5*(t4+c2)+c2*t4)
            der(2,9)=d32/d3*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
            der(2,8)=d128/d3*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
            der(2,7)=d64*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
            der(2,6)=d128/d3*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
            der(2,5)=-d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
            der(2,4)=-d128/d3*c1*(t8*(t7+c3)+c3*t7)
            der(2,3)=-d64*c1*t1*(t7+c3)
            der(2,2)=-d128/d3*c1*t1*t2
            der(2,13)=d128*c1*t1*(c3-c2)
            der(2,15)=d128*c1*(c3*(t4+c2)-c2*t4)
            der(2,14)=d128*c1*(c3*t7-c2*(c3+t7))
        CASE (4)
            der(1,1)=-etam
            der(1,2)=-etap
            der(1,3)=etap
            der(1,4)=etam
            der(2,1)=-xim
            der(2,2)=xim
            der(2,3)=xip
            der(2,4)=-xip
        CASE(5)
            der(1,1)=-etam+pt5*xi*(one-eta**2)
            der(1,2)=-etap+pt5*xi*(one-eta**2)
            der(1,3)=etap+pt5*xi*(one-eta**2)
            der(1,4)=etam+pt5*xi*(one-eta**2)
            der(1,5)=-two*xi*(one-eta**2)
            der(2,1)=-xim+pt5*eta*(one-xi**2)
            der(2,2)=xim+pt5*eta*(one-xi**2)
            der(2,3)=xip+pt5*eta*(one-xi**2)
            der(2,4)=-xip+pt5*eta*(one-xi**2)
            der(2,5)=-two*eta*(one-xi**2)
        CASE(8)
            der(1,1)=etam*(two*xi+eta)
            der(1,2)=-d8*etam*etap
            der(1,3)=etap*(two*xi-eta)
            der(1,4)=-d4*etap*xi
            der(1,5)=etap*(two*xi+eta)
            der(1,6)=d8*etap*etam
            der(1,7)=etam*(two*xi-eta)
            der(1,8)=-d4*etam*xi
            der(2,1)=xim*(xi+two*eta)
            der(2,2)=-d4*xim*eta
            der(2,3)=xim*(two*eta-xi)
            der(2,4)=d8*xim*xip
            der(2,5)=xip*(xi+two*eta)
            der(2,6)=-d4*xip*eta
            der(2,7)=xip*(two*eta-xi)
            der(2,8)=-d8*xim*xip
        CASE(9)
            etam=eta-one
            etap=eta+one
            xim=xi-one
            xip=xi+one
            der(1,1)=pt25*x2m1*eta*etam
            der(1,2)=-pt5*x2m1*etap*etam
            der(1,3)=pt25*x2m1*eta*etap
            der(1,4)=-xi*eta*etap
            der(1,5)=pt25*x2p1*eta*etap
            der(1,6)=-pt5*x2p1*etap*etam
            der(1,7)=pt25*x2p1*eta*etam
            der(1,8)=-xi*eta*etam
            der(1,9)=two*xi*etap*etam
            der(2,1)=pt25*xi*xim*e2m1
            der(2,2)=-xi*xim*eta
            der(2,3)=pt25*xi*xim*e2p1
            der(2,4)=-pt5*xip*xim*e2p1
            der(2,5)=pt25*xi*xip*e2p1
            der(2,6)=-xi*xip*eta
            der(2,7)=pt25*xi*xip*e2m1
            der(2,8)=-pt5*xip*xim*e2m1
            der(2,9)=two*xip*xim*eta
        CASE DEFAULT
            WRITE(*,*)"wrong number of nodes in shape_der"
        END SELECT
    CASE(3)  ! d3 dimensional elements
    xi=points(i,1)
    eta=points(i,2)
    zeta=points(i,3)
    etam=one-eta
    xim=one-xi
    zetam=one-zeta
    etap=eta+one
    xip=xi+one
    zetap=zeta+one
    SELECT CASE(nod)
        CASE(4)
            der(1:3,1:4)=zero
            der(1,1)=one
            der(2,2)=one
            der(3,3)=one
            der(1,4)=-one
            der(2,4)=-one
            der(3,4)=-one
        CASE(8)
            der(1,1)=-pt125*etam*zetam
            der(1,2)=-pt125*etam*zetap
            der(1,3)= pt125*etam*zetap
            der(1,4)= pt125*etam*zetam
            der(1,5)=-pt125*etap*zetam
            der(1,6)=-pt125*etap*zetap
            der(1,7)= pt125*etap*zetap
            der(1,8)= pt125*etap*zetam
            der(2,1)=-pt125*xim*zetam
            der(2,2)=-pt125*xim*zetap
            der(2,3)=-pt125*xip*zetap
            der(2,4)=-pt125*xip*zetam
            der(2,5)= pt125*xim*zetam
            der(2,6)= pt125*xim*zetap
            der(2,7)= pt125*xip*zetap
            der(2,8)= pt125*xip*zetam
            der(3,1)=-pt125*xim*etam
            der(3,2)= pt125*xim*etam
            der(3,3)= pt125*xip*etam
            der(3,4)=-pt125*xip*etam
            der(3,5)=-pt125*xim*etap
            der(3,6)= pt125*xim*etap
            der(3,7)= pt125*xip*etap
            der(3,8)=-pt125*xip*etap
        CASE(14) ! type 6 element
            der(1,1)= (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*         &
            (eta-one)*(zeta-one)/d8
            der(1,2)=-(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*         &
            (eta-one)*(zeta+one)/d8
            der(1,3)=-(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*         &
            (eta-one)*(zeta+one)/d8
            der(1,4)= (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*         &
            (eta-one)*(zeta-one)/d8
            der(1,5)= -(eta-one)*(zeta+one)*(zeta-one)*xi
            der(1,6)=-(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
            der(1,7)=  (eta+one)*(eta-one)*(zeta+one)*xi
            der(1,8)= (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
            der(1,9)= -(eta+one)*(eta-one)*(zeta-one)*xi
            der(1,10)= (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*        &
            (eta+one)*(zeta-one)/d8
            der(1,11)=-(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*        &
            (eta+one)*(zeta+one)/d8
            der(1,12)=-(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*        &
            (eta+one)*(zeta+one)/d8
            der(1,13)= (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*        &
            (eta+one)*(zeta-one)/d8
            der(1,14)=  (eta+one)*(zeta+one)*(zeta-one)*xi
            der(2,1)= (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*         &
            (xi-one)*(zeta-one)/d8
            der(2,2)=-(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*         &
            (xi-one)*(zeta+one)/d8
            der(2,3)=-(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*         &
            (xi+one)*(zeta+one)/d8
            der(2,4)= (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*         &
            (xi+one)*(zeta-one)/d8
            der(2,5)=-(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
            der(2,6)= -(xi-one)*(zeta+one)*(zeta-one)*eta
            der(2,7)=  (xi+one)*(xi-one)*(zeta+one)*eta
            der(2,8)=  (xi+one)*(zeta+one)*(zeta-one)*eta
            der(2,9)= -(xi+one)*(xi-one)*(zeta-one)*eta
            der(2,10)= (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*        &
            (xi-one)*(zeta-one)/d8
            der(2,11)=-(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*        &
            (xi-one)*(zeta+one)/d8
            der(2,12)=-(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*        &
            (xi+one)*(zeta+one)/d8
            der(2,13)= (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*        &
            (xi+one)*(zeta-one)/d8
            der(2,14)= (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
            der(3,1)= (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*         &
            (xi-one)*(eta-one)/d8
            der(3,2)=-(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*         &
            (xi-one)*(eta-one)/d8
            der(3,3)=-(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*         &
            (xi+one)*(eta-one)/d8
            der(3,4)= (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*         &
            (xi+one)*(eta-one)/d8
            der(3,5)= -(xi+one)*(xi-one)*(eta-one)*zeta
            der(3,6)= -(xi-one)*(eta+one)*(eta-one)*zeta
            der(3,7)= (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
            der(3,8)=  (xi+one)*(eta+one)*(eta-one)*zeta
            der(3,9)=-(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
            der(3,10)= (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*        &
            (xi-one)*(eta+one)/d8
            der(3,11)=-(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*        &
            (xi-one)*(eta+one)/d8
            der(3,12)=-(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*        &
            (xi+one)*(eta+one)/d8
            der(3,13)= (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*        &
            (xi+one)*(eta+one)/d8
            der(3,14)=  (xi+one)*(xi-one)*(eta+one)*zeta
        CASE(20)
            xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
            etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
            zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
            DO l=1,20
            xi0=xi*xii(l)
            eta0=eta*etai(l)
            zeta0=zeta*zetai(l)
            IF(l==4.OR.l==8.OR.l==16.OR.l==20)THEN
                der(1,l)=-pt5*xi*(one+eta0)*(one+zeta0)
                der(2,l)=pt25*etai(l)*(one-xi*xi)*(one+zeta0)
                der(3,l)=pt25*zetai(l)*(one-xi*xi)*(one+eta0)
            ELSE IF(l>=9.AND.l<=12)THEN
                der(1,l)=pt25*xii(l)*(one-eta*eta)*(one+zeta0)
                der(2,l)=-pt5*eta*(one+xi0)*(one+zeta0)
                der(3,l)=pt25*zetai(l)*(one+xi0)*(one-eta*eta)
            ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
                der(1,l)=pt25*xii(l)*(one+eta0)*(one-zeta*zeta)
                der(2,l)=pt25*etai(l)*(one+xi0)*(one-zeta*zeta)
                der(3,l)=-pt5*zeta*(one+xi0)*(one+eta0)
            ELSE
                der(1,l)=pt125*xii(l)*(one+eta0)*(one+zeta0)*                   &
                (two*xi0+eta0+zeta0-one)
                der(2,l)=pt125*etai(l)*(one+xi0)*(one+zeta0)*                   &
                (xi0+two*eta0+zeta0-one)
                der(3,l)=pt125*zetai(l)*(one+xi0)*(one+eta0)*                   &
                (xi0+eta0+two*zeta0-one)
            END IF
            END DO
        CASE DEFAULT
            WRITE(*,*)"wrong number of nodes in shape_der"
        END SELECT
    CASE DEFAULT
        WRITE(*,*)"wrong number of dimensions in shape_der"
    END SELECT
    RETURN
  END SUBROUTINE SHAPE_DER


  SUBROUTINE INVAR(stress,sigm,dsbar,theta)
    !
    ! This subroutine forms the stress INVARiants in 2- or 3-d. iii
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::stress(:)
    REAL(iwp),INTENT(OUT),OPTIONAL::sigm,dsbar,theta
    REAL(iwp)::sx,sy,sz,txy,dx,dy,dz,xj3,sine,s1,s2,s3,s4,s5,s6,ds1,ds2,ds3, &
    d2,d3,sq3,zero=0.0_iwp,small=1.e-10_iwp,one=1.0_iwp,two=2.0_iwp,       &
    three=3.0_iwp,six=6.0_iwp,thpt5=13.5_iwp
    INTEGER::nst
    nst=UBOUND(stress,1)
    SELECT CASE(nst)
    CASE(4)
      sx=stress(1)
      sy=stress(2)
      txy=stress(3)
      sz=stress(4)
      sigm=(sx+sy+sz)/three
      dsbar=SQRT((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+six*txy**2)/SQRT(two)
      IF(dsbar<small)THEN
        theta=zero
      ELSE
        dx=(two*sx-sy-sz)/three
        dy=(two*sy-sz-sx)/three
        dz=(two*sz-sx-sy)/three
        xj3=dx*dy*dz-dz*txy**2
        sine=-thpt5*xj3/dsbar**3
        IF(sine>=one)sine=one
        IF(sine<-one)sine=-one
        theta=ASIN(sine)/three
      END IF
    CASE DEFAULT
      WRITE(*,*)"wrong size for nst in invar"
    END SELECT

    RETURN
  END SUBROUTINE invar 


  SUBROUTINE INVERT(matrix)
    !
    ! This subroutine inverts a small square matrix onto itself.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN OUT)::matrix(:,:)
    REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
    INTEGER::ndim,i,k
    ndim=UBOUND(matrix,1)
    IF(ndim==2)THEN
        det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
        j11=matrix(1,1)
        matrix(1,1)=matrix(2,2)
        matrix(2,2)=j11
        matrix(1,2)=-matrix(1,2)
        matrix(2,1)=-matrix(2,1)
        matrix=matrix/det
    ELSE IF(ndim==3)THEN
        det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
        det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
        det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
        j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
        j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
        j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
        j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
        j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
        j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
        j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
        j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
        j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
        matrix(1,1)=j11
        matrix(1,2)=j12
        matrix(1,3)=j13
        matrix(2,1)=j21
        matrix(2,2)=j22
        matrix(2,3)=j23
        matrix(3,1)=j31
        matrix(3,2)=j32
        matrix(3,3)=j33
        matrix=matrix/det
    ELSE
        DO k=1,ndim
            con=matrix(k,k)
            matrix(k,k)=1.0_iwp
            matrix(k,:)=matrix(k,:)/con
            DO i=1,ndim
                IF(i/=k)THEN
                    con=matrix(i,k)
                    matrix(i,k)=0.0_iwp
                    matrix(i,:)=matrix(i,:)-matrix(k,:)*con
                END IF
            END DO
        END DO
    END IF
    RETURN
  END SUBROUTINE INVERT


  
!                                                                              !
!**************************** MATERIAL FUNCTIONS ******************************!
!                                                                              !

  SUBROUTINE PLASTICMOD(k,Maxplastic,SModulus,h)
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::k,Maxplastic,SModulus
    REAL(iwp),INTENT(OUT)::h

    IF(k>=Maxplastic) h=0.0
    IF(k<Maxplastic)  h=SModulus
    RETURN
  END SUBROUTINE PLASTICMOD


  SUBROUTINE VMFLOW(stress,dsbar,vmfl)
    !
    ! This subroutine forms the von-Mises flow vector.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::stress(:),dsbar
    REAL(iwp),INTENT(OUT)::vmfl(:)
    REAL(iwp)::sigm,onept5=1.5_iwp,two=2.0_iwp,d3=3.0_iwp
    sigm=(stress(1)+stress(2)+stress(4))/d3
    vmfl(1)=stress(1)-sigm
    vmfl(2)=stress(2)-sigm
    vmfl(3)=stress(3)*two
    vmfl(4)=stress(4)-sigm
    vmfl=vmfl*onept5/dsbar
    RETURN
  END SUBROUTINE VMFLOW

  
!                                                                              !
!******************************* MPM FUNCTIONS ********************************!
!                                                                              !

  SUBROUTINE LOCATE_MP_LOCAL(mpoints,coord,gm_coord)
    !
    ! local coordinate of the material points
    !
    ! PARAMETER
    ! ---------
    !   mpoints: real(iwp)
    !     Material point local coordinate. Output value.
    !   coord : real(iwp) (node per element, dimension)
    !     Reference element global coords. 1st column x-ord, 2nd column y-ord 
    !   gm_coord : real(iwp) (number of mps, dimension)
    !     Material point global coordinate.
    !
    IMPLICIT NONE
    INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
    REAL(iwp), INTENT(IN)  :: coord(:,:), gm_coord(:,:)
    REAL(iwp), INTENT(OUT) :: mpoints(:,:)
    REAL(iwp), ALLOCATABLE:: fun(:), der(:,:),jac(:,:)
    REAL(iwp), ALLOCATABLE:: guess_local(:,:), guess_global(:,:)
    INTEGER   :: count, count_limit, nod, ndim, nmps
    REAL(iwp) :: tolerance1, tolerance2, beta
    REAL(iwp) :: zero=0.0_iwp, pt9=0.9_iwp
    LOGICAL:: find_local
        
    nod  = UBOUND(coord,1) 
    ndim = UBOUND(coord,2)
    nmps = UBOUND(gm_coord,2)

    ALLOCATE(fun(nod),der(ndim,nod),jac(ndim,ndim),                             &
             guess_local(1,ndim),guess_global(1,ndim))

    find_local  = .FALSE.
    guess_local = zero
    CALL SHAPE_FUN(fun, guess_local, 1)
    guess_global(1,:)=MATMUL(fun,coord)

    count=0
    count_limit = 100
    DO WHILE (count < count_limit)
      count = count+1
      beta  = pt9
      ! Evaluate Error of Global Coordinates
      guess_global(1,:)=gm_coord(:,1)-guess_global(1,:)

      ! Adjust Guess_local with the gradient of the error
      CALL SHAPE_DER(der,guess_local,1)
      jac = MATMUL(der,coord)
      CALL INVERT(jac)
      guess_local(1,:)=guess_local(1,:)+beta*MATMUL(TRANSPOSE(jac),guess_global(1,:))
      
      ! Re-interpolate Guess_Global with the new Guess_local Value
      CALL shape_fun(fun,guess_local,1)
      guess_global(1,:)=MATMUL(fun,coord)

      ! Check Convergence
      tolerance1=ABS(MAXVAL(guess_global(1,:)-gm_coord(:,1)))
      tolerance2=ABS(MINVAL(guess_global(1,:)-gm_coord(:,1)))
      IF(tolerance1<=1.0D-10.and.tolerance2<=1.0D-10) EXIT
    END DO
    mpoints(1,:)=guess_local(1,:)

    DEALLOCATE(fun,guess_local,guess_global,der,jac)
    RETURN
  END SUBROUTINE LOCATE_MP_LOCAL


  SUBROUTINE COUNT_MP_IN_ELEMENT(c_ele,k_ele,etype,nels,nmps,a_ele)
    !
    ! count the number of material points inside an element.
    !
    ! PARAMETER
    ! ---------
    !   c_ele: List of number of MPs inside a given element
    !   k_ele: accumulation of c_ele
    !   etype: List of element indicators. 1: ele has no MP, 2: ele has MP 
    !   a_ele: List of 
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::a_ele(:)
    INTEGER,INTENT(IN)::nels, nmps
    INTEGER,INTENT(OUT)::c_ele(:), k_ele(0:), etype(:)
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
  END SUBROUTINE COUNT_MP_IN_ELEMENT
  

!                                                                              !
!******************************* GIMP FUNCTIONS *******************************!
!                                                                              !

  SUBROUTINE GIMP_GET_SUPPORT_ELEMENTS(member_elements,g_coord,gm_coord,       &
    c_ele,nf,g_num,nels,lp,tol,neighbour)
    !
    ! This subroutine save all the elements inside 
    ! the support domain of each material point
    !
    IMPLICIT NONE
    INTEGER, PARAMETER    :: iwp=SELECTED_REAL_KIND(15)
    INTEGER, INTENT(IN)   :: nels
    INTEGER, INTENT(IN)   :: g_num(:,:), nf(:,:), c_ele(:)
    REAL(iwp), INTENT(IN) :: g_coord(:,:), gm_coord(:,:), lp(:,:)
    INTEGER, INTENT(OUT)  :: member_elements(:,:) ! shape(n_max_element, nmps)
    ! Optional Arguments
    REAL(iwp), INTENT(IN), OPTIONAL :: tol
    INTEGER, INTENT(IN), OPTIONAL :: neighbour ! Number of external elements 
                                               ! affected by characteristic 
                                               ! functions
    ! Local Parameters
    INTEGER   :: i,k,nmps,count
    REAL(iwp) :: lpx, lpy ! particle size in both directions 
                          ! (Bardenhagen et al., 2004) page 485
    ! Default Arguments
    INTEGER :: d_neighbour=2
    REAL(iwp) :: d_tol = 1.0e-5
    IF(present(neighbour)) d_neighbour = neighbour
    IF(present(tol)) d_tol = tol

    nmps=UBOUND(gm_coord,2)
    member_elements=0 

    SELECT CASE(d_neighbour) 
    CASE(1) ! No neighbour
      ! Commented out until a locating algorithm without relying on 
      ! cell geometry is available. One options would be using LOCATE_MP_LOCAL 
      ! but this means high computational cost when doing the search (AS) 
      !
      ! IF (present(cellsize)) THEN
      !   count=1
      !   DO k=1,nmps
      !     colx=gm_coord(1,k)/cellsize+1.0
      !     rowy=ABS(gm_coord(2,k))/cellsize+2.0
      !     ielloc=(rowy-1.0)*nx1+colx
      !     member_elements(k,count)=ielloc
      !   END DO  
      ! END IF

    CASE(2) ! Immediate Cell Neighbour
      DO k=1,nmps
        lpx=lp(1,k)
        lpy=lp(2,k)
        count=1
        Ele:DO i=1,nels
          ! Check if the element is active (4 nodes of the element free at 
          ! least in one direction or with a material point inside)
          ACTIVE_ELEMENT: IF(                                                  &
             ((nf(1,g_num(1,i))>0.or.nf(2,g_num(1,i))>0).and.                  &
              (nf(1,g_num(2,i))>0.or.nf(2,g_num(2,i))>0).and.                  &
              (nf(1,g_num(3,i))>0.or.nf(2,g_num(3,i))>0).and.                  &
              (nf(1,g_num(4,i))>0.or.nf(2,g_num(4,i))>0))                      &
             .or. c_ele(i)>0 ) THEN

            ! Check whether the MP is inside the x support domain of one element
            IF((gm_coord(1,k)<MAXVAL(g_coord(1,g_num(:,i)))+lpx+d_tol).and.      &
               (gm_coord(1,k)>MINVAL(g_coord(1,g_num(:,i)))-lpx-d_tol)) THEN 
              ! Now check in y-dir support domain
              IF((gm_coord(2,k)<MAXVAL(g_coord(2,g_num(:,i)))+lpy+d_tol).and.    &
                 (gm_coord(2,k)>MINVAL(g_coord(2,g_num(:,i)))-lpy-d_tol))THEN
                ! Element in Support Domain! put its index on member_elements 
                member_elements(k,count)=i
                count=count+1
              END IF
            END IF
          END IF ACTIVE_ELEMENT
          IF (count > 4) EXIT ! Early exit
        END DO Ele
      END DO

    CASE DEFAULT
      WRITE(*,*)'wrong number of neigbouring elements for support elements'
    END SELECT
    RETURN
  END SUBROUTINE GIMP_GET_SUPPORT_ELEMENTS


  SUBROUTINE GIMP_GET_PARTICLE_G(s,g_vector,nf,nodes_v,nodes_n)
    !
    ! Get the steering vector of the nodes within a material point
    ! support domain
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::s, nf(:,:)
    INTEGER,INTENT(IN)::nodes_v(:,:), nodes_n
    INTEGER,INTENT(OUT)::g_vector(:)
    INTEGER::i,n
    g_vector=0
    n=1
    DO i=1, nodes_n
        g_vector(n:n+1)=nf(:,nodes_v(i,s))
        n=n+2
    END DO
    RETURN 
  END SUBROUTINE GIMP_GET_PARTICLE_G

  
  SUBROUTINE GIMP_GET_SHAPE_FUNCTION_AND_DERIVATIVES(func,deriv,s,iel,lm_coord, &
    lp,cellsize,g_coord,gm_coord,support_nodes,a_ele,gimptol)
    !     
    ! Subroutine to compute implicit GIMP shape functions and its derivatives
    ! for a given particle `s` w.r.t. a given element `iel`.
    !
    ! iGIMP_funder3 uses a constant gradient at the center of the nodes 
    ! and it reduces as the support domain ends
    !
    ! TODO: Shape Function Derivatives is still not functional
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::s, iel
    INTEGER,INTENT(IN)::support_nodes(:,:), a_ele(:)
    REAL(iwp),INTENT(IN)::cellsize(:,:),lp(:,:)
    REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),lm_coord(:,:)
    REAL(iwp),INTENT(IN),OPTIONAL::gimptol
    REAL(iwp),INTENT(OUT)::func(:),deriv(:,:)
    ! local variable
    REAL(iwp)::node_coordx,node_coordy,mp_coordx,mp_coordy
    REAL(iwp)::Sr,Sz,dSr,dSz,factx,facty
    REAL(iwp)::offset
    REAL(iwp)::lpx,lpy
    INTEGER::i,j,k,m,n,nod,nodes
    ! Constants
    REAL(iwp)::one=1.0_iwp,two=2.0_iwp
    REAL(iwp)::tol=1.0e-5
    if(present(gimptol)) tol=gimptol
    func=0.0; deriv=0.0 ! Initiate global shape function and its derivates
    i=1; k=1; m=1

    ! Get Material Point Coordinate
    mp_coordx=gm_coord(1,s) !--x coordinate of current material point 's'
    mp_coordy=gm_coord(2,s) !--y coordinate of material point 's'
    
    ! Determine Particle Size
    nodes=UBOUND(deriv,2)
    lpx=lp(1,s)*(two/cellsize(1,s))
    lpy=lp(2,s)*(two/cellsize(2,s))

    Fun_Deriv:DO n=1,nodes
        Sr=0.0; Sz=0.0      ! Initiate local shape function
        dSr=0.0; dSz=0.0    ! Initiate local shape function derivatives
        nod=support_nodes(n,1) ! Get nodes within the MP support domain

        ! Obtain nodal coordinates
        node_coordx=g_coord(1,nod)
        node_coordy=g_coord(2,nod)

        !*** LOCAL COORDINATES ***!
        ! Evaluate the local coordinate of the material point inside 
        ! and outside the current evaluated affected element in x-dir
        IF(a_ele(s)==iel)THEN 
          ! Material point is inside the element
          ! set factx from local mp coodinate relative to the element
          factx = lm_coord(s,1) 
        ELSE     
          ! Material point is outside the element
          IF(mp_coordx<node_coordx+tol)THEN          
            ! point outside the element in the x direction (left side)
            IF(n==1.or.n==2) offset=node_coordx-mp_coordx
            IF(n==3.or.n==4) offset=(node_coordx-cellsize(1,s))-mp_coordx
            ! local coordinate point at the left side of the element
            ! calculated as the ratio of the distance between `the center of the
            ! element to point` and `the center of the element to the element boundary`
            factx = ((cellsize(1,s)/two+offset)/(cellsize(1,s)/two))*(-one)   
          ELSE IF(mp_coordx>node_coordx-tol)THEN     
            !-- point outside the element in the x direction (right side)
            IF(n==1.or.n==2) offset=mp_coordx-(node_coordx+cellsize(1,s))
            IF(n==3.or.n==4) offset=mp_coordx-(node_coordx)
            ! local coordinate point at the right side of the element
            factx = ((cellsize(1,s)/two+offset)/(cellsize(1,s)/two))       
          ELSE
            WRITE(*, *) "cannot locate particle. Abort program"
            PAUSE
            CALL EXIT(1)
          END IF
        END IF

        !-- Evaluate local coordinate again in y-dir
        IF(a_ele(s)==iel)THEN
            facty=lm_coord(s,2)
        ELSE 
            IF(mp_coordy>=node_coordy)THEN        
              ! point outside the element in the y direction (above the element)
              IF(n==1.or.n==4)offset=((node_coordy+cellsize(2,s))-mp_coordy)*(-one)
              IF(n==2.or.n==3)offset=(node_coordy-mp_coordy)*(-one)
              facty=((cellsize(2,s)/two+offset)/(cellsize(2,s)/two))
            ELSE IF(mp_coordy<=node_coordy)THEN     
              ! point outside the element in the y direction (below the element)
              IF(n==2.or.n==3)offset=(mp_coordy-(node_coordy-cellsize(2,s)))*(-one)
              IF(n==1.or.n==4)offset=(mp_coordy-node_coordy)*(-one)
              facty=((cellsize(2,s)/two+offset)/(cellsize(2,s)/two))*(-one)
            END IF  
        END IF

        
        !*** Shape functions ***!
        IF(n==1)THEN        !Comparing with node 1
          IF(factx<=-0.50_iwp)THEN
              Sr=0.5_iwp*(factx+lpx)*(1.0_iwp-0.5_iwp*(factx+lpx))+0.75_iwp
          ELSE IF(factx>=-0.50_iwp.and.factx<=0.5_iwp)THEN    
              Sr=0.5_iwp*(factx+lpx)-0.25_iwp*(factx+lpx)**2-0.5_iwp*(factx-lpx)+0.25_iwp*(factx-lpx)**2
          ELSE
              Sr=0.25_iwp-0.5_iwp*(factx-lpx)+0.25_iwp*(factx-lpx)**2
          END IF   

          IF(facty<=-0.50_iwp)THEN
              Sz=0.5_iwp*(facty+lpy)*(1.0_iwp-0.5_iwp*(facty+lpy))+0.75_iwp
          ELSE IF(facty>=-0.50_iwp.and.facty<=0.5_iwp)THEN    
              Sz=0.5_iwp*(facty+lpy)-0.25_iwp*(facty+lpy)**2-0.5_iwp*(facty-lpy)+0.25_iwp*(facty-lpy)**2
          ELSE
              Sz=0.25_iwp-0.5_iwp*(facty-lpy)+0.25_iwp*(facty-lpy)**2
          END IF 

        ELSE IF(n==2)THEN   !Comparing with node 2
          IF(factx<=-0.50_iwp)THEN
              Sr=0.5_iwp*(factx+lpx)*(1.0_iwp-0.5_iwp*(factx+lpx))+0.75_iwp
          ELSE IF(factx>=-0.50_iwp.and.factx<=0.5_iwp)THEN    
              Sr=0.5_iwp*(factx+lpx)-0.25_iwp*(factx+lpx)**2-0.5_iwp*(factx-lpx)+0.25_iwp*(factx-lpx)**2
          ELSE
              Sr=0.25_iwp-0.5_iwp*(factx-lpx)+0.25_iwp*(factx-lpx)**2
          END IF 

          IF(facty<=-0.50_iwp)THEN
              Sz=0.5_iwp*(facty+lpy)*(1.0_iwp+0.5_iwp*(facty+lpy))+0.25_iwp
          ELSE IF(facty>=-0.50_iwp.and.facty<=0.5_iwp)THEN    
              Sz=0.5_iwp*(facty+lpy)+0.25_iwp*(facty+lpy)**2-0.5_iwp*(facty-lpy)-0.25_iwp*(facty-lpy)**2
          ELSE
              Sz=0.75_iwp-0.5_iwp*(facty-lpy)-0.25_iwp*(facty-lpy)**2
          END IF       

        ELSE IF(n==3)THEN   !Comparing with node 3
          IF(factx<=-0.50_iwp)THEN
              Sr=0.5_iwp*(factx+lpx)*(1.0_iwp+0.5_iwp*(factx+lpx))+0.25_iwp
          ELSE IF(factx>=-0.50_iwp.and.factx<=0.5_iwp)THEN    
              Sr=0.5_iwp*(factx+lpx)+0.25_iwp*(factx+lpx)**2-0.5_iwp*(factx-lpx)-0.25_iwp*(factx-lpx)**2
          ELSE
              Sr=0.75_iwp-0.5_iwp*(factx-lpx)-0.25_iwp*(factx-lpx)**2
          END IF 

          IF(facty<=-0.50_iwp)THEN
              Sz=0.5_iwp*(facty+lpy)*(1.0_iwp+0.5_iwp*(facty+lpy))+0.25_iwp
          ELSE IF(facty>=-0.50_iwp.and.facty<=0.5_iwp)THEN    
              Sz=0.5_iwp*(facty+lpy)+0.25_iwp*(facty+lpy)**2-0.5_iwp*(facty-lpy)-0.25_iwp*(facty-lpy)**2
          ELSE
              Sz=0.75_iwp-0.5_iwp*(facty-lpy)-0.25_iwp*(facty-lpy)**2
          END IF 

        ELSE IF(n==4)THEN   !Comparing with node 4
          IF(factx<=-0.50_iwp)THEN
              Sr=0.5_iwp*(factx+lpx)*(1.0_iwp+0.5_iwp*(factx+lpx))+0.25_iwp
          ELSE IF(factx>=-0.50_iwp.and.factx<=0.5_iwp)THEN    
              Sr=0.5_iwp*(factx+lpx)+0.25_iwp*(factx+lpx)**2-0.5_iwp*(factx-lpx)-0.25_iwp*(factx-lpx)**2
          ELSE
              Sr=0.75_iwp-0.5_iwp*(factx-lpx)-0.25_iwp*(factx-lpx)**2
          END IF 

          IF(facty<=-0.50_iwp)THEN
              Sz=0.5_iwp*(facty+lpy)*(1.0_iwp-0.5_iwp*(facty+lpy))+0.75_iwp
          ELSE IF(facty>=-0.50_iwp.and.facty<=0.5_iwp)THEN    
              Sz=0.5_iwp*(facty+lpy)-0.25_iwp*(facty+lpy)**2-0.5_iwp*(facty-lpy)+0.25_iwp*(facty-lpy)**2
          ELSE
              Sz=0.25_iwp-0.5_iwp*(facty-lpy)+0.25_iwp*(facty-lpy)**2
          END IF 
        END IF  

        !*** Shape Function Derivatives ***!
        ! ...
        
        ! Assign Calculated Shape Function and its Derivatives to Array
        func(m)=Sr*Sz
        ! IF(abs(func(m))>2)THEN
        !     Sr=Sr
        !     PAUSE
        ! END IF    
        ! IF((Sr<=0.0.or.Sz<=0.0))THEN
        !     deriv(1,m)=0.0
        !     deriv(2,m)=0.0
        ! ELSE  
        !     deriv(1,m)=dSr*Sz
        !     deriv(2,m)=dSz*Sr
        ! END IF

        m=m+1

    END DO Fun_Deriv

    RETURN
  END SUBROUTINE GIMP_GET_SHAPE_FUNCTION_AND_DERIVATIVES


!                                                                              !
!******************************* CMPM FUNCTIONS *******************************!
!                                                                              !
  
  SUBROUTINE CONST(EQUATIONS,bound,m)

    ! this soubrutine is to compute the constants of a 9x9 or 16x16
    ! matrix, considering also the position of the element

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(OUT)::EQUATIONS(:,:)
    INTEGER,INTENT(IN)::bound,m
    INTEGER::x,y,fil,col

    equations=0.0

    IF(bound==5)THEN
      x=-3;y=3;fil=1;col=1
      DO fil=1,16
        DO col=1,16
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
          IF(col==10)EQUATIONS(fil,col)=x**3
          IF(col==11)EQUATIONS(fil,col)=y**3
          IF(col==12)EQUATIONS(fil,col)=x**3*y
          IF(col==13)EQUATIONS(fil,col)=x*y**3
          IF(col==14)EQUATIONS(fil,col)=x**3*y**2
          IF(col==15)EQUATIONS(fil,col)=x**2*y**3
          IF(col==16)EQUATIONS(fil,col)=x**3*y**3
        END DO
        x=x+2
        IF(x>3)THEN
          x=-3
          y=y-2
        END IF
      END DO
    END IF

    IF((bound==1.or.bound==4.or.bound==2).and.m==1)THEN
      x=-1;y=1;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>3)THEN
          x=-1
          y=y-2
        END IF
      END DO
    END IF

    IF(bound==4.and.m==2)THEN
      x=-1;y=3;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>3)THEN
          x=-1
          y=y-2
        END IF
      END DO
    END IF

    IF(bound==2.and.m==2)THEN
      x=-3;y=1;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>1)THEN
          x=-3
          y=y-2
        END IF
      END DO
    END IF

    IF((bound==3.or.bound==6).and.m==1)THEN
      x=-3;y=1;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>1)THEN
          x=-3
          y=y-2
        END IF
      END DO
    END IF

    IF(bound==6.and.m==2)THEN
      x=-3;y=3;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>1)THEN
          x=-3
          y=y-2
        END IF
      END DO
    END IF

    IF((bound==7.or.bound==8).and.m==1)THEN
      x=-1;y=3;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>3)THEN
          x=-1
          y=y-2
        END IF
      END DO
    END IF

    IF(bound==8.and.m==2)THEN
      x=-3;y=3;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>1)THEN
          x=-3
          y=y-2
        END IF
      END DO
    END IF

    IF(bound==9)THEN
      x=-3;y=3;fil=1;col=1
      DO fil=1,9
        DO col=1,9
          IF(col==1)EQUATIONS(fil,col)=1.0
          IF(col==2)EQUATIONS(fil,col)=x
          IF(col==3)EQUATIONS(fil,col)=y
          IF(col==4)EQUATIONS(fil,col)=y*x
          IF(col==5)EQUATIONS(fil,col)=x**2
          IF(col==6)EQUATIONS(fil,col)=y**2
          IF(col==7)EQUATIONS(fil,col)=y*x**2
          IF(col==8)EQUATIONS(fil,col)=x*y**2
          IF(col==9)EQUATIONS(fil,col)=x**2*y**2
        END DO
        x=x+2
        IF(x>1)THEN
          x=-3
          y=y-2
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE CONST


!                                                                              !
!*************************** SPARSE MATRIX ROUTINES ***************************!
!                                                                              !

  SUBROUTINE FKDIAG(kdiag,g)
    !
    ! This subroutine computes the skyline profile.
    !
    IMPLICIT NONE
    INTEGER,INTENT(IN)::g(:)
    INTEGER,INTENT(IN OUT)::kdiag(:)
    INTEGER::idof,i,iwp1,j,im,k
    idof=SIZE(g)
    DO i=1,idof
      iwp1=1
      IF(g(i)/=0)THEN
        DO j=1,idof
          IF(g(j)/=0)THEN
            im=g(i)-g(j)+1
            IF(im>iwp1)iwp1=im
          END IF
        END DO
        k=g(i)
        IF(iwp1>kdiag(k))kdiag(k)=iwp1
      END IF
    END DO
    RETURN
  END SUBROUTINE FKDIAG


  SUBROUTINE FSPARV(kv,km,g,kdiag)
    !
    ! This subroutine assles element matrices into a symmetric skyline
    ! global matrix.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::g(:),kdiag(:)
    REAL(iwp),INTENT(IN)::km(:,:)
    REAL(iwp),INTENT(IN OUT)::kv(:)
    INTEGER::i,idof,k,j,iw,ival
    idof=UBOUND(g,1)
    DO i=1,idof
      k=g(i)
      IF(k/=0)THEN
        DO j=1,idof
          IF(g(j)/=0)THEN
            iw=k-g(j)
            IF(iw>=0)THEN
              ival=kdiag(k)-iw
              kv(ival)=kv(ival)+km(i,j)
            END IF
          END IF
        END DO
      END IF
    END DO
    RETURN
  END SUBROUTINE FSPARV


  SUBROUTINE LINMUL_SKY(kv,disps,loads,kdiag)
    !
    ! This subroutine forms the product of symmetric matrix stored as
    ! a skyline and a vector.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::kv(:),disps(0:)
    REAL(iwp),INTENT(OUT)::loads(0:)
    INTEGER,INTENT(IN)::kdiag(:)
    INTEGER::n,i,j,low,lup,k
    REAL(iwp)::x,zero=0.0_iwp
    n=UBOUND(disps,1)
    DO i=1,n
      x=zero
      lup=kdiag(i)
      IF(i==1) low=lup
      IF(i/=1) low=kdiag(i-1)+1
      DO j=low,lup
        x=x+kv(j)*disps(i+j-lup)
      END DO
      loads(i)=x
      IF(i==1) CYCLE ! Go to next iteration
      lup=lup-1
      DO j=low,lup
        k=i+j-lup-1
        loads(k)=loads(k)+kv(j)*disps(i)
      END DO
    END DO
    RETURN
  END SUBROUTINE LINMUL_SKY


  SUBROUTINE SPABAC(kv,loads,kdiag)
    !
    ! This subroutine performs Cholesky forward and back-substitution
    ! on a symmetric skyline global matrix.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::kv(:)
    REAL(iwp),INTENT(IN OUT)::loads(0:)
    INTEGER,INTENT(IN)::kdiag(:)
    INTEGER::n,i,ki,l,m,j,it,k
    REAL(iwp)::x
    n=UBOUND(kdiag,1)
    loads(1)=loads(1)/kv(1)
    DO i=2,n
      ki=kdiag(i)-i
      l=kdiag(i-1)-ki+1
      x=loads(i)
      IF(l/=i)THEN
        m=i-1
        DO j=l,m
          x=x-kv(ki+j)*loads(j)
        END DO
      END IF
      loads(i)=x/kv(ki+i)
    END DO
    DO it=2,n
      i=n+2-it
      ki=kdiag(i)-i
      x=loads(i)/kv(ki+i)
      loads(i)=x
      l=kdiag(i-1)-ki+1
      IF(l/=i)THEN
        m=i-1
        DO k=l,m
        loads(k)=loads(k)-x*kv(ki+k)
        END DO
      END IF
    END DO
    loads(1)=loads(1)/kv(1)
    RETURN
  END SUBROUTINE SPABAC


  SUBROUTINE SPARIN(kv,kdiag)
    !
    ! This subroutine performs Cholesky factorisation on a symmetric
    ! skyline global matrix.
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN OUT)::kv(:)
    INTEGER,INTENT(IN)::kdiag(:)
    INTEGER::n,i,ki,l,kj,j,ll,m,k
    REAL(iwp)::x
    n=UBOUND(kdiag,1)
    kv(1)=SQRT(kv(1))
    DO i=2,n
      ki=kdiag(i)-i
      l=kdiag(i-1)-ki+1
        DO j=l,i
          x=kv(ki+j)
          kj=kdiag(j)-j
          IF(j/=1)THEN
            ll=kdiag(j-1)-kj+1
            ll=max(l,ll)
            IF(ll/=j)THEN
              m=j-1
              DO k=ll,m
                x=x-kv(ki+k)*kv(kj+k)
              END DO
            END IF
          END IF
          kv(ki+j)=x/kv(kj+j)
        END DO
      kv(ki+i)=SQRT(x)
    END DO

    RETURN
  END SUBROUTINE SPARIN


!                                                                              !
!****************************** UNSORTED & DRAFTS *****************************!
!                                                                              !


!                                                                              !
!******************************************************************************!
END MODULE FUNCTIONS