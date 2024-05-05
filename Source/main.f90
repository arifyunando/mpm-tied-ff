module main

     SAVE

     CONTAINS

     SUBROUTINE num_to_g(num,nf,g)
     !elmat
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
     END SUBROUTINE num_to_g

SUBROUTINE elmat(area,rho,emm,lumpval)
!
! This subroutine forms the "analytical" lumped mass matrix for
! quadrilateral 4- or 8-node plane strain elements.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::area,rho,lumpval
 REAL(iwp),INTENT(OUT)::emm(:,:)
 REAL(iwp)::zero=0.0_iwp,pt2=0.2_iwp,pt25=0.25_iwp
 INTEGER::i,ndof
 ndof=UBOUND(emm,1)
 emm=zero
 SELECT CASE(ndof)
 CASE(8)
   DO i=1,8
     emm(i,i)=pt25*area*rho
   END DO
 CASE(16)
   DO i=1,16
     emm(i,i)=lumpval*area*rho
   END DO
   DO i=1,13,4
     emm(i,i)=(pt25-lumpval)*area*rho
   END DO
   DO i=2,14,4
     emm(i,i)=(pt25-lumpval)*area*rho
   END DO
 CASE DEFAULT
   WRITE(*,*)"Wrong number of nodes for rectangular element"
 END SELECT
RETURN
END SUBROUTINE elmat

SUBROUTINE elmat2(emm)
!
! This subroutine forms the "analytical" lumped mass matrix for
! quadrilateral 4- or 8-node plane strain elements.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(OUT)::emm(:,:)
 REAL(iwp)::zero=0.0_iwp,pt2=0.2_iwp,pt25=0.25_iwp
 INTEGER::i,ndof
 ndof=UBOUND(emm,1)
 emm=zero
 SELECT CASE(ndof)
 CASE(8)
   DO i=1,8
     emm(i,i)=pt25
   END DO
 CASE(16)
   DO i=1,16
     emm(i,i)=pt2
   END DO
   DO i=1,13,4
     emm(i,i)=pt25*emm(3,3)
   END DO
   DO i=2,14,4
     emm(i,i)=pt25*emm(3,3)
   END DO
 CASE DEFAULT
   WRITE(*,*)"Wrong number of nodes for rectangular element"
 END SELECT
RETURN
END SUBROUTINE elmat2

SUBROUTINE elmatvis(rho,emm)
!
! This subroutine forms the "analytical" lumped mass matrix for
! quadrilateral 4- or 8-node plane strain elements.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::rho
 REAL(iwp),INTENT(OUT)::emm(:,:)
 REAL(iwp)::zero=0.0_iwp,pt2=0.2_iwp,pt25=0.25_iwp
 INTEGER::i,ndof
 ndof=UBOUND(emm,1)
 emm=zero
 SELECT CASE(ndof)
 CASE(8)
   DO i=1,8
     emm(i,i)=pt25*rho
   END DO
 CASE(16)
   DO i=1,16
     emm(i,i)=pt2*rho
   END DO
   DO i=1,13,4
     emm(i,i)=pt25*emm(3,3)
   END DO
   DO i=2,14,4
     emm(i,i)=pt25*emm(3,3)
   END DO
 CASE DEFAULT
   WRITE(*,*)"Wrong number of nodes for rectangular element"
 END SELECT
RETURN
END SUBROUTINE elmatvis

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



     SUBROUTINE fkdiag(kdiag,g)
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
     END SUBROUTINE fkdiag

     SUBROUTINE linmul_sky(kv,disps,loads,kdiag)
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
      IF(i==1)low=lup
      IF(i/=1)low=kdiag(i-1)+1
      DO j=low,lup
        x=x+kv(j)*disps(i+j-lup)
      END DO
      loads(i)=x
      IF(i==1)CYCLE
      lup=lup-1
      DO j=low,lup
        k=i+j-lup-1
        loads(k)=loads(k)+kv(j)*disps(i)
      END DO
    END DO
   RETURN
   END SUBROUTINE linmul_sky

     SUBROUTINE formnf(nf)
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
     END SUBROUTINE formnf

          SUBROUTINE formnf_tied(nf,nex,ney)
     !
     ! This subroutine forms the nf matrix.
     !
      IMPLICIT NONE
      INTEGER,INTENT(IN OUT)::nf(:,:),nex,ney
      INTEGER::i,j,m,nodx,nody,k,n,col_ini,col_end,nod
      m=0
      k=0
      j=0
      col_ini=0
      col_end=nex*2+1
      nodx=nex*2+1
      nody=ney*2+1

      DO i=1,nodx
          k=k+1
          IF(k==1)THEN
              col_ini=col_ini+1
              nod=(col_ini-1)*nody+1
              k=k+1
          ELSE
              col_end=col_end-1
              nod=col_end*nody+1
              k=0
          END IF


      DO j=1,nody
              DO n=1,UBOUND(nf,1)
              IF(nf(n,nod)/=0)THEN
                m=m+1
                nf(n,nod)=m
              END IF
              END DO
              nod=nod+1
      END DO
      END DO

      !DO j=1,UBOUND(nf,2)
      !  DO i=1,UBOUND(nf,1)
      !    IF(nf(i,j)/=0)THEN
      !      m=m+1
      !      nf(i,j)=m
      !    END IF
      !  END DO
      !END DO

     RETURN
     END SUBROUTINE formnf_tied

     SUBROUTINE mesh(g_coord,g_num,argv,nlen,ips)
     !
     ! This subroutine produces a PostScript output file "*.msh" displaying
     ! the undeformed finite element mesh.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::g_coord(:,:)
      INTEGER,INTENT(IN)::g_num(:,:),ips,nlen
      CHARACTER(*),INTENT(IN)::argv
      REAL(iwp)::xmin,xmax,ymin,ymax,width,height,scale=72,sxy,xo,yo,x,y,      &
        pt5=0.5_iwp,opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,ept5=8.5_iwp,         &
        d11=11.0_iwp
      INTEGER::i,ii,j,jj,nn,nod,nel
      OPEN(ips,FILE=argv(1:nlen)//'.msh')
     !
     !                       compute size of mesh
     !
      nn=UBOUND(g_coord,2)
      xmin=g_coord(1,1)
      xmax=g_coord(1,1)
      ymin=g_coord(2,1)
      ymax=g_coord(2,1)
      DO i=2,nn
        IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
        IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
        IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
        IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
      END DO
      width =xmax-xmin
      height=ymax-ymin
     !
     !                       allow 1.5" margin minimum on each side of figure
     !
      IF(height.GE.d11/ept5*width)THEN
     !
     !                       height governs the scale
     !
        sxy=scale*d8/height
        xo=scale*pt5*(ept5-d8*width/height)
        yo=scale*opt5
      ELSE
     !
     !                       width governs the scale
     !
        sxy=scale*fpt5/width
        xo=scale*opt5
        yo=scale*pt5*(d11-fpt5*height/width)
      END IF
     !
     !                       start PostScript output
     !
      WRITE(ips,'(a)')'%!PS-Adobe-1.0'
      WRITE(ips,'(a)')'%%DocumentFonts: none'
      WRITE(ips,'(a)')'%%Pages: 1'
      WRITE(ips,'(a)')'%%EndComments'
      WRITE(ips,'(a)')'/m {moveto} def'
      WRITE(ips,'(a)')'/l {lineto} def'
      WRITE(ips,'(a)')'/s {stroke} def'
      WRITE(ips,'(a)')'/c {closepath} def'
      WRITE(ips,'(a)')'%%EndProlog'
      WRITE(ips,'(a)')'%%Page: 0 1'
      WRITE(ips,'(a)')'gsave'
      WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
      WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
     !
     !                       draw the mesh
     !
      nod=UBOUND(g_num,1)
      nel=UBOUND(g_num,2)
      IF(nod==5)nod=4
      IF(nod==9)nod=8
      IF(nod==10)nod=9
      IF(nod==15)nod=12
      DO i=1,nel
        ii=g_num(1,i)
        IF(ii==0)CYCLE
        x=sxy*(g_coord(1,ii)-xmin)
        y=sxy*(g_coord(2,ii)-ymin)
        WRITE(ips,'(2f9.2,a)')x,y,' m'
        DO j=2,nod
          jj=g_num(j,i)
          x=sxy*(g_coord(1,jj)-xmin)
          y=sxy*(g_coord(2,jj)-ymin)
          WRITE(ips,'(2f9.2,a)') x, y,' l'
        END DO
        WRITE(ips,'(a)')'c s'
      END DO
     !
     !                       close output file
     !
      WRITE(ips,'(a)')'grestore'
      WRITE(ips,'(a)')'showpage'
      CLOSE(ips)
     !
     RETURN
     END SUBROUTINE mesh


     SUBROUTINE sample2(element,s,wt)
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

      nip=UBOUND(s,1)

      SELECT CASE(element)

      CASE('quadrilateral')
        SELECT CASE(nip)
        CASE(1)
          s(1,1)=0.0_iwp
          s(1,2)=0.0_iwp
          wt(1)=4.0_iwp
         CASE(4)
           s(1,1)=-0.50_iwp
           s(1,2)= 0.50_iwp
           s(2,1)= 0.50_iwp
           s(2,2)= 0.50_iwp
           s(3,1)=-0.50_iwp
           s(3,2)=-0.50_iwp
           s(4,1)= 0.50_iwp
           s(4,2)=-0.50_iwp
           wt=1.0_iwp
!$$$$$$          CASE(4)
!$$$$$$            s(1,1)=-1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(1,2)= 1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(2,1)= 1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(2,2)= 1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(3,1)=-1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(3,2)=-1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(4,1)= 1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            s(4,2)=-1.0_iwp/sqrt(3.0_iwp)
!$$$$$$            wt=1.0_iwp
        CASE(9)
          s(1:7:3,1)=-sqrt(5.0_iwp)/3.0_iwp!0.666666666_iwp
          s(2:8:3,1)=0.0_iwp
          s(3:9:3,1)=sqrt(5.0_iwp)/3.0_iwp
          s(1:3,2)  =sqrt(5.0_iwp)/3.0_iwp
          s(4:6,2)  =0.0_iwp
          s(7:9,2)  =-sqrt(5.0_iwp)/3.0_iwp
          wt= 0.4444444444_iwp
        CASE(16)
          s(1:13:4,1)=-0.75_iwp
          s(2:14:4,1)=-0.25_iwp
          s(3:15:4,1)= 0.25_iwp
          s(4:16:4,1)= 0.75_iwp
          s(1:4,2)   = 0.75_iwp
          s(5:8,2)   = 0.25_iwp
          s(9:12,2)  =-0.25_iwp
          s(13:16,2) =-0.75_iwp
          wt = 0.250_iwp

        END SELECT
      END SELECT
     RETURN
     END SUBROUTINE sample2

     SUBROUTINE shape_fun_1D(fun,funf,points,i)
     !
     !   This subroutine computes the values of the shape functions.
     !   to local coordinates
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(in)::i
      REAL(iwp),INTENT(IN)::points(:,:)
      REAL(iwp),INTENT(OUT)::fun(:,:),funf(:,:)
      REAL(iwp)::xi
      INTEGER::ndim,nod
      REAL,PARAMETER::one=1.0_iwp,two=2.0_iwp,zero=0.0_iwp
      ndim=1
      nod=UBOUND(fun,1)
      fun=zero
      funf=zero
      SELECT CASE(ndim)
      CASE(1) ! one dimensional case
        xi=points(i,1)
        SELECT CASE(nod)
        CASE(2)

          fun(1,1)=0.5
          fun(2,1)=0.5
        CASE(3)
          fun(1,1)=(xi**2-xi)/two
          fun(2,1)=one-xi**2
          fun(3,1)=(xi**2+xi)/two
          funf(1,1)=(one-xi)/two
          funf(2,1)=(one+xi)/two
        CASE(4)


        CASE(5)

        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      CASE(2) ! two dimensional case

        SELECT CASE(nod)


        CASE(9)

        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      CASE(3) ! d3 dimensional case

        SELECT CASE(nod)


        CASE(20)


        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      CASE DEFAULT
        WRITE(*,*)"wrong number of dimensions in shape_fun"
      END SELECT
     RETURN
     END SUBROUTINE shape_fun_1D

      SUBROUTINE shape_der_1D(der,points,i)
     !
     !   This subroutine computes the values of the shape functions.
     !   to local coordinates
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(in)::i
      REAL(iwp),INTENT(IN)::points(:,:)
      REAL(iwp),INTENT(OUT)::der(:,:)
      REAL(iwp)::xi
      INTEGER::ndim,nod
      REAL,PARAMETER::one=1.0_iwp,two=2.0_iwp,zero=0.0_iwp
      ndim=1
      nod=UBOUND(der,2)
      der=zero
      SELECT CASE(ndim)
      CASE(1) ! one dimensional case
        xi=points(i,1)
        SELECT CASE(nod)
        CASE(2)
          der(1,1)=0.5
          der(2,1)=0.5
        CASE(3)
          der(1,1)=(xi)-one/two
          der(1,2)=-two*xi
          der(1,3)=(xi)+one/two
        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      
      CASE DEFAULT
        WRITE(*,*)"wrong number of dimensions in shape_fun"
      END SELECT
     RETURN
      END SUBROUTINE shape_der_1D
      
      SUBROUTINE shape_der_1D_left(der,points,i)
     !
     !   This subroutine computes the values of the shape functions.
     !   to local coordinates
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(in)::i
      REAL(iwp),INTENT(IN)::points(:,:)
      REAL(iwp),INTENT(OUT)::der(:,:)
      REAL(iwp)::xi
      INTEGER::ndim,nod
      REAL,PARAMETER::one=1.0_iwp,two=2.0_iwp,zero=0.0_iwp
      ndim=1
      nod=UBOUND(der,2)
      der=zero
      SELECT CASE(ndim)
      CASE(1) ! one dimensional case
        xi=points(i,1)
        SELECT CASE(nod)
        CASE(2)
          der(1,1)=0.5
          der(2,1)=0.5
        CASE(3)
          der(1,1)=(xi)+one/two
          der(1,2)=-two*xi
          der(1,3)=(xi)-one/two
        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      
      CASE DEFAULT
        WRITE(*,*)"wrong number of dimensions in shape_fun"
      END SELECT
     RETURN
      END SUBROUTINE shape_der_1D_left
      
    SUBROUTINE shape_der_int(der,points,i)
     !
     !   This subroutine computes the values of the shape functions.
     !   to local coordinates
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(in)::i
      REAL(iwp),INTENT(IN)::points(:,:)
      REAL(iwp),INTENT(OUT)::der(:,:)
      REAL(iwp)::xi
      INTEGER::ndim,nod
      REAL,PARAMETER::one=1.0_iwp,two=2.0_iwp,zero=0.0_iwp
      ndim=1
      nod=UBOUND(der,2)
      der=zero
      SELECT CASE(ndim)
      CASE(1) ! one dimensional case
        xi=points(i,1)
        SELECT CASE(nod)
        CASE(2)

        CASE(3)
          der(1,1)=(xi)-one/two
          der(1,2)=-two*xi
          der(1,3)=(xi)+one/two
    

        CASE(5)

        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      CASE(2) ! two dimensional case

        SELECT CASE(nod)


        CASE(9)

        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      CASE(3) ! d3 dimensional case

        SELECT CASE(nod)


        CASE(20)


        CASE DEFAULT
          WRITE(*,*)"wrong number of nodes in shape_fun"
        END SELECT
      CASE DEFAULT
        WRITE(*,*)"wrong number of dimensions in shape_fun"
      END SELECT
     RETURN
      END SUBROUTINE shape_der_int

       SUBROUTINE formdee(Young,poisson,dee)
!
! This subroutine forms the stiffness matrix of a
! beam element (bending only).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::Young,poisson
 REAL(iwp),INTENT(OUT)::dee(:,:)
 REAL(iwp)::two=2.0_iwp,lamda,mu,one=1.0_iwp,zero=0.0_iwp

 dee=zero
 lamda=Young*poisson/((one+poisson)*(one-two*poisson))
 mu=Young/(two*(one+poisson))

 dee(1,1)=lamda+two*mu
 dee(2,1)=lamda
 dee(1,2)=dee(2,1)
 dee(2,2)=dee(1,1)
 dee(3,3)=mu

RETURN
       END SUBROUTINE formdee

       SUBROUTINE formbee(bee,size,nod,deriv)
!
! This subroutine forms the stiffness matrix of a
! beam element (bending only).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::size,deriv(:,:)
 REAL(iwp),INTENT(OUT)::bee(:,:)
 INTEGER,INTENT(IN)::nod
 REAL(iwp)::two=2.0_iwp,one=1.0_iwp,zero=0.0_iwp

 bee=zero

 IF(nod==2)THEN
 bee(3,1)=-one
 bee(2,2)=-one
 bee(3,3)=one
 bee(2,4)=one
 bee=bee/size
 ELSE IF(nod==3)THEN
    bee(3,1)=deriv(1,1)
    bee(2,2)=deriv(1,1)
    bee(3,3)=deriv(1,2)
    bee(2,4)=deriv(1,2)
    bee(3,5)=deriv(1,3)
    bee(2,6)=deriv(1,3)
 END IF

RETURN
       END SUBROUTINE formbee
       
              SUBROUTINE formbee_int(bee,size,nod,deriv)
!
! This subroutine forms the stiffness matrix of a
! beam element (bending only).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::size,deriv(:,:)
 REAL(iwp),INTENT(OUT)::bee(:,:)
 INTEGER,INTENT(IN)::nod
 REAL(iwp)::two=2.0_iwp,one=1.0_iwp,zero=0.0_iwp

 bee=zero

 IF(nod==2)THEN

 ELSE IF(nod==3)THEN
    bee(1,1)=deriv(1,1)
    bee(3,1)=deriv(1,1)
    bee(2,2)=deriv(1,1)
    bee(3,3)=deriv(1,2)
    bee(1,3)=deriv(1,2)
    bee(2,4)=deriv(1,2)
    bee(3,5)=deriv(1,3)
    bee(1,5)=deriv(1,3)
    bee(2,6)=deriv(1,3)
 END IF

RETURN
    END SUBROUTINE formbee_int

     SUBROUTINE sample3(element,s,wt)
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
      REAL::three=3.0_iwp,five=5.0_iwp

      nip=UBOUND(s,1)

      SELECT CASE(element)

      CASE('quadrilateral')
        SELECT CASE(nip)
        CASE(1)
          s(1,1)=1.0_iwp
          s(1,2)=1.0_iwp
          wt(1)=4.0_iwp
         CASE(3)
          s(1,1)=sqrt(three/five)
          s(1,2)=0.0_iwp
          s(2,1)=0.0_iwp
          s(2,2)=0.0_iwp
          s(3,1)=-sqrt(three/five)
          s(3,2)=0.0_iwp
          wt(1) = 0.555555555555556_iwp
          wt(2) = 0.888888888888889_iwp
          wt(3) = 0.555555555555556_iwp

        CASE(4)
          s(1,1)=-0.00_iwp
          s(1,2)= 1.00_iwp
          s(2,1)= 1.00_iwp
          s(2,2)= -0.00_iwp
          s(3,1)=-0.415_iwp
          s(3,2)=-0.415_iwp
          s(4,1)= 1.00_iwp
          s(4,2)=1.00_iwp
          wt=1.0_iwp
        CASE(9)
          nip=6
          s(1:7:3,1)=-0.666666666_iwp
          s(2:8:3,1)=0.0_iwp
          s(3:9:3,1)=0.666666666_iwp
          s(1:3,2)  =0.666666666_iwp
          s(4:6,2)  =0.0_iwp
          s(7:9,2)  =-0.666666666_iwp
          wt= 0.4444444444_iwp
        CASE(16)
          nip=10
          s(1:13:4,1)=-0.75_iwp
          s(2:14:4,1)=-0.25_iwp
          s(3:15:4,1)= 0.25_iwp
          s(4:16:4,1)= 0.75_iwp
          s(1:4,2)   = 0.75_iwp
          s(5:8,2)   = 0.25_iwp
          s(9:12,2)  =-0.25_iwp
          s(13:16,2) =-0.75_iwp
          wt = 0.250_iwp

        END SELECT
      END SELECT
     RETURN
     END SUBROUTINE sample3


     SUBROUTINE deemat(dee,e,v)
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
     END SUBROUTINE deemat

     SUBROUTINE shape_fun(fun,points,i)
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
      CASE(1) ! one dimensional case
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
      CASE(2) ! two dimensional case
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
      CASE(3) ! d3 dimensional case
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
     END SUBROUTINE shape_fun

     SUBROUTINE bee8(bee,coord,xi,eta,det)
     !
     ! Analytical version of the bee matrix for an 8-node quad
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::coord(:,:),xi,eta
      REAL(iwp),INTENT(OUT)::bee(:,:),det
      REAL::x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8,xi2,xi3,et2,et3,   &
        xi2et2,xi1et1,xi2et1,xi1et2,xy12,xy13,xy14,xy15,xy16,xy17,xy18,xy23,   &
        xy24,xy25,xy26,xy27,xy28,xy34,xy35,xy36,xy37,xy38,xy45,xy46,xy47,xy48, &
        xy56,xy57,xy58,xy67,xy68,xy78,xcu,xsq,xon,ecu,esq,eon,x2e2,x1e1,x2e1,  &
        x1e2,cons,bot,top,dn1dx,dn2dx,dn3dx,dn4dx,dn5dx,dn6dx,dn7dx,dn8dx,     &
        dn1dy,dn2dy,dn3dy,dn4dy,dn5dy,dn6dy,dn7dy,dn8dy,zero=0.0_iwp,          &
        one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,d6=6.0_iwp,   &
        d8=8.0_iwp,pt125=0.125_iwp
      x1=coord(1,1)
      x2=coord(2,1)
      x3=coord(3,1)
      x4=coord(4,1)
      x5=coord(5,1)
      x6=coord(6,1)
      x7=coord(7,1)
      x8=coord(8,1)
      y1=coord(1,2)
      y2=coord(2,2)
      y3=coord(3,2)
      y4=coord(4,2)
      y5=coord(5,2)
      y6=coord(6,2)
      y7=coord(7,2)
      y8=coord(8,2)
      xi2=xi*xi
      xi3=xi2*xi
      et2=eta*eta
      et3=et2*eta
      xi2et2=xi2*et2
      xi1et1=xi*eta
      xi2et1=xi2*eta
      xi1et2=xi*et2
      xy12=x1*y2-y1*x2
      xy13=x1*y3-y1*x3
      xy14=x1*y4-y1*x4
      xy15=x1*y5-y1*x5
      xy16=x1*y6-y1*x6
      xy17=x1*y7-y1*x7
      xy18=x1*y8-y1*x8
      xy23=x2*y3-y2*x3
      xy24=x2*y4-y2*x4
      xy25=x2*y5-y2*x5
      xy26=x2*y6-y2*x6
      xy27=x2*y7-y2*x7
      xy28=x2*y8-y2*x8
      xy34=x3*y4-y3*x4
      xy35=x3*y5-y3*x5
      xy36=x3*y6-y3*x6
      xy37=x3*y7-y3*x7
      xy38=x3*y8-y3*x8
      xy45=x4*y5-y4*x5
      xy46=x4*y6-y4*x6
      xy47=x4*y7-y4*x7
      xy48=x4*y8-y4*x8
      xy56=x5*y6-y5*x6
      xy57=x5*y7-y5*x7
      xy58=x5*y8-y5*x8
      xy67=x6*y7-y6*x7
      xy68=x6*y8-y6*x8
      xy78=x7*y8-y7*x8
      xcu=-d8*xy48+d4*(-xy14+xy38+xy47+xy58)+two*(xy13+xy15-xy37-xy57)
      xsq=two*(-xy13+xy14-xy17+xy18+xy24-xy28-xy34+xy35-xy38-xy45+xy46+xy47-   &
        xy57+xy58+xy68-xy78)+one*(-xy12+xy16-xy23-xy25+xy27-xy36-xy56-xy67)
      xon=d8*xy48+two*(xy14-xy18+xy34-xy38-xy45-xy47-xy58-xy78)+xy12-xy16+xy23-&
        xy25+xy27+xy36-xy56-xy67
      ecu=-d8*xy26+d4*(xy16+xy25+xy36+xy27)+two*(-xy15-xy17-xy35-xy37)
      esq=two*(-xy12+xy13-xy16+xy17-xy23+xy24+xy25-xy27-xy28-xy35+xy36+xy46-   &
        xy56+xy57-xy67+xy68)-xy14+xy18-xy34+xy38-xy45-xy47-xy58-xy78
      eon=d8*xy26+two*( xy12-xy16-xy23-xy25-xy27-xy36-xy56+xy67)+xy14-xy18-    &
        xy34+xy38-xy45+xy47-xy58+xy78
      x2e2=d6*(xy24-xy28+xy46+xy68)+d3*(-xy12+xy13-xy14+xy16-xy17+xy18-xy23-   &
        xy25+xy27-xy34+xy35-xy36+xy38-xy45-xy47-xy56+xy57-xy58-xy67-xy78)
      x1e1=d8*(-xy24-xy28+xy46-xy68)+d6*(-xy12+xy18+xy23+xy34-xy45-xy56+xy67+  &
        xy78)+two*(xy14-xy16+xy25+xy27-xy36+xy38-xy47+xy58)
      x2e1=d8*(xy24+xy28+xy46-xy68)+d5*( xy17-xy18-xy34+xy35-xy45+xy78)+d4*    &
        (xy12-xy16-xy23-xy25-xy27-xy36-xy56+xy67)+d3*(-xy14+xy15+xy37-xy38-    &
        xy47+xy58)
      x1e2=d8*(-xy24+xy28+xy46+xy68)+d5*(xy12-xy13+xy23+xy57-xy56-xy67)+d4*    &
        (xy14-xy18+xy34-xy38-xy45-xy47-xy58-xy78)+d3*(-xy15+xy16+xy25-xy27-    &
        xy36+xy37)
      cons= two*(-xy24+xy28-xy46-xy68)
     !
      bot=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      bot=bot+xi2et2*x2e2+xi1et1*x1e1
      bot=bot+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      det=pt125*bot
     !
      xcu=two*y3-d4*y4+two*y5
      xsq=-y2-two*y3+two*y4+y6-two*y7+two*y8
      xon=y2+two*y4 -y6-two*y8
      ecu=-two*y5+d4*y6-two*y7
      esq=-two*y2+two*y3-y4-two*y6+two*y7+y8
      eon=two*y2+y4-two*y6-y8
      x2e2=-d3*y2+d3*y3-d3*y4+d3*y6-d3*y7+d3*y8
      x1e1=-d6*y2+two*y4-two*y6+d6*y8
      x2e1=d4*y2-d3*y4+d3*y5-d4*y6+d5*y7-d5*y8
      x1e2=d5*y2-d5*y3+d4*y4-d3*y5+d3*y6-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn1dx=top/bot
     !
      xcu=zero
      xsq=y1-y3+two*y4-y5+y7-two*y8
      xon=-y1+y3-y5+y7
      ecu=d4*y5-d8*y6+d4*y7
      esq=two*y1-two*y3+two*y4+two*y5-two*y7-two*y8
      eon=-two*y1-two*y3-two*y5+d8*y6-two*y7
      x2e2= d3*y1-d3*y3+d6*y4-d3*y5+d3*y7-d6*y8
      x1e1= d6*y1+d6*y3-d8*y4+two*y5+two*y7-d8*y8
      x2e1=-d4*y1-d4*y3+d8*y4-d4*y5-d4*y7+d8*y8
      x1e2=-d5*y1+d5*y3-d8*y4+d3*y5-d3*y7+d8*y8
      cons=-two*y4+two*y8
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn2dx=top/bot
     !
      xcu=-two*y1-two*y7+d4*y8
      xsq=two*y1+y2-two*y4+two*y5-y6-two*y8
      xon=-y2+two*y4+y6-two*y8
      ecu=-two*y5+d4*y6-two*y7
      esq=-two*y1+two*y2-y4-two*y5+two*y6+y8
      eon=two*y2-y4-two*y6+y8
      x2e2=-d3*y1+d3*y2-d3*y4+d3*y5-d3*y6+d3*y8
      x1e1=-d6*y2+d6*y4-two*y6+two*y8
      x2e1=+d4*y2-d5*y4+d5*y5-d4*y6+d3*y7-d3*y8
      x1e2=d5*y1-d5*y2+d4*y4-d3*y6+d3*y7-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn3dx=top/bot
     !
      xcu=d4*y1+d4*y7-d8*y8
      xsq=-two*y1-two*y2+two*y3-two*y5+two*y6+two*y7
      xon=-two*y1-two*y3-two*y5-two*y7+d8*y8
      ecu=zero
      esq=y1-two*y2+y3-y5+two*y6-y7
      eon=-y1+y3-y5+y7
      x2e2=d3*y1-d6*y2+d3*y3-d3*y5+d6*y6-d3*y7
      x1e1=-two*y1+d8*y2-d6*y3-d6*y5+d8*y6-two*y7
      x2e1=d3*y1-d8*y2+d5*y3-d5*y5+d8*y6-d3*y7
      x1e2=-d4*y1+d8*y2-d4*y3-d4*y5+d8*y6-d4*y7
      cons=two*y2-two*y6
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn4dx=top/bot
     !
      xcu=-two*y1-two*y7+d4*y8
      xsq=y2-two*y3+two*y4-y6-two*y7+two*y8
      xon=y2 +two*y4-y6-two*y8
      ecu=two*y1-d4*y2+two*y3
      esq=-two*y2+two*y3+y4-two*y6+two*y7-y8
      eon=two*y2+y4-two*y6-y8
      x2e2=d3*y2-d3*y3+d3*y4-d3*y6+d3*y7-d3*y8
      x1e1=-two*y2+d6*y4-d6*y6+two*y8
      x2e1=-d3*y1+d4*y2-d5*y3+d5*y4-d4*y6+d3*y8
      x1e2=d3*y1-d3*y2+d4*y4-d5*y6+d5*y7-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn5dx=top/bot
     !
      xcu=zero
      xsq=-y1+y3-two*y4+y5-y7+two*y8
      xon=y1-y3+y5-y7
      ecu=-d4*y1+d8*y2-d4*y3
      esq=two*y1-two*y3-two*y4+two*y5-two*y7+two*y8
      eon=two*y1-d8*y2+two*y3+two*y5+two*y7
      x2e2=-d3*y1+d3*y3-d6*y4+d3*y5-d3*y7+d6*y8
      x1e1=two*y1+two*y3-d8*y4+d6*y5+d6*y7-d8*y8
      x2e1=d4*y1+d4*y3-d8*y4+d4*y5+d4*y7-d8*y8
      x1e2=-d3*y1+d3*y3-d8*y4+d5*y5-d5*y7+d8*y8
      cons=two*y4-two*y8
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn6dx=top/bot
     !
      xcu=two*y3-d4*y4+two*y5
      xsq=two*y1-y2-two*y4+two*y5+y6-two*y8
      xon=-y2+two*y4+y6-two*y8
      ecu=two*y1-d4*y2+two*y3
      esq=-two*y1+two*y2+y4-two*y5+two*y6-y8
      eon=two*y2-y4-two*y6+y8
      x2e2=d3*y1-d3*y2+d3*y4-d3*y5+d3*y6-d3*y8
      x1e1=-two*y2+two*y4-d6*y6+d6*y8
      x2e1=-d5*y1+d4*y2-d3*y3+d3*y4-d4*y6+d5*y8
      x1e2=d3*y2-d3*y3+d4*y4-d5*y5+d5*y6-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn7dx=top/bot
     !
      xcu=-d4*y3+d8*y4-d4*y5
      xsq=-two*y1+two*y2+two*y3-two*y5-two*y6+two*y7
      xon=two*y1+two*y3-d8*y4+two*y5+two*y7
      ecu=zero
      esq=-y1+two*y2-y3+y5-two*y6+y7
      eon=y1-y3+y5-y7
      x2e2=-d3*y1+d6*y2-d3*y3+d3*y5-d6*y6+d3*y7
      x1e1=-d6*y1+d8*y2-two*y3-two*y5+d8*y6-d6*y7
      x2e1=d5*y1-d8*y2+d3*y3-d3*y5+d8*y6-d5*y7
      x1e2=d4*y1-d8*y2+d4*y3+d4*y5-d8*y6+d4*y7
      cons=-two*y2+two*y6
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn8dx=top/bot
     !
      y1=-x1
      y2=-x2
      y3=-x3
      y4=-x4
      y5=-x5
      y6=-x6
      y7=-x7
      y8=-x8
     !
      xcu=two*y3-d4*y4+two*y5
      xsq=-y2-two*y3+two*y4+y6-two*y7+two*y8
      xon=y2+two*y4-y6-two*y8
      ecu=-two*y5+d4*y6-two*y7
      esq=-two*y2+two*y3-y4-two*y6+two*y7+y8
      eon=two*y2+y4-two*y6-y8
      x2e2=-d3*y2+d3*y3-d3*y4 +d3*y6-d3*y7+d3*y8
      x1e1=-d6*y2+two*y4-two*y6+d6*y8
      x2e1=d4*y2-d3*y4+d3*y5-d4*y6+d5*y7-d5*y8
      x1e2=d5*y2-d5*y3+d4*y4-d3*y5+d3*y6-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn1dy=top/bot
     !
      xcu=zero
      xsq=y1-y3+two*y4-y5+y7-two*y8
      xon=-y1+y3-y5+y7
      ecu=d4*y5-d8*y6+d4*y7
      esq=two*y1-two*y3+two*y4+two*y5-two*y7-two*y8
      eon=-two*y1-two*y3-two*y5+d8*y6-two*y7
      x2e2=d3*y1-d3*y3+d6*y4-d3*y5+d3*y7-d6*y8
      x1e1=d6*y1+d6*y3-d8*y4+two*y5+two*y7-d8*y8
      x2e1=-d4*y1-d4*y3+d8*y4-d4*y5-d4*y7+d8*y8
      x1e2=-d5*y1+d5*y3-d8*y4+d3*y5-d3*y7+d8*y8
      cons=-two*y4+two*y8
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn2dy=top/bot
     !
      xcu=-two*y1-two*y7+d4*y8
      xsq=two*y1+y2-two*y4+two*y5-y6-two*y8
      xon=-y2+two*y4+y6-two*y8
      ecu=-two*y5+d4*y6-two*y7
      esq=-two*y1+two*y2-y4-two*y5+two*y6+y8
      eon=two*y2-y4-two*y6+y8
      x2e2=-d3*y1+d3*y2-d3*y4+d3*y5-d3*y6+d3*y8
      x1e1=-d6*y2+d6*y4-two*y6+two*y8
      x2e1=d4*y2-d5*y4+d5*y5-d4*y6+d3*y7-d3*y8
      x1e2=d5*y1-d5*y2+d4*y4-d3*y6+d3*y7-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn3dy=top/bot
     !
      xcu=d4*y1+d4*y7-d8*y8
      xsq=-two*y1-two*y2+two*y3-two*y5+two*y6+two*y7
      xon=-two*y1-two*y3-two*y5-two*y7+d8*y8
      ecu=zero
      esq=y1-two*y2+y3-y5+two*y6-y7
      eon=-y1+y3-y5+y7
      x2e2=d3*y1-d6*y2+d3*y3-d3*y5+d6*y6-d3*y7
      x1e1=-two*y1+d8*y2-d6*y3-d6*y5+d8*y6-two*y7
      x2e1=d3*y1-d8*y2+d5*y3-d5*y5+d8*y6-d3*y7
      x1e2=-d4*y1+d8*y2-d4*y3-d4*y5+d8*y6-d4*y7
      cons=two*y2-two*y6
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn4dy=top/bot
     !
      xcu=-two*y1-two*y7+d4*y8
      xsq=y2-two*y3+two*y4-y6-two*y7+two*y8
      xon=y2+two*y4-y6-two*y8
      ecu=two*y1-d4*y2+two*y3
      esq=-two*y2+two*y3+y4-two*y6+two*y7-y8
      eon=two*y2+y4-two*y6-y8
      x2e2=d3*y2-d3*y3+d3*y4-d3*y6+d3*y7-d3*y8
      x1e1=-two*y2+d6*y4-d6*y6+two*y8
      x2e1=-d3*y1+d4*y2-d5*y3+d5*y4-d4*y6 +d3*y8
      x1e2=d3*y1-d3*y2+d4*y4-d5*y6+d5*y7-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn5dy=top/bot
     !
      xcu=zero
      xsq=-y1+y3-two*y4+y5-y7+two*y8
      xon=y1-y3+y5-y7
      ecu=-d4*y1+d8*y2-d4*y3
      esq=two*y1-two*y3-two*y4+two*y5-two*y7+two*y8
      eon=two*y1-d8*y2+two*y3+two*y5+two*y7
      x2e2=-d3*y1+d3*y3-d6*y4+d3*y5-d3*y7+d6*y8
      x1e1=two*y1+two*y3-d8*y4+d6*y5+d6*y7-d8*y8
      x2e1=d4*y1+d4*y3-d8*y4+d4*y5+d4*y7-d8*y8
      x1e2=-d3*y1+d3*y3-d8*y4+d5*y5-d5*y7+d8*y8
      cons=two*y4-two*y8
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn6dy=top/bot
     !
      xcu=two*y3-d4*y4+two*y5
      xsq=two*y1-y2-two*y4+two*y5+y6-two*y8
      xon=-y2+two*y4+y6-two*y8
      ecu=two*y1-d4*y2+two*y3
      esq=-two*y1+two*y2+y4-two*y5+two*y6-y8
      eon=two*y2-y4-two*y6+y8
      x2e2=d3*y1-d3*y2+d3*y4-d3*y5+d3*y6-d3*y8
      x1e1=-two*y2+two*y4-d6*y6+d6*y8
      x2e1=-d5*y1+d4*y2-d3*y3+d3*y4-d4*y6+d5*y8
      x1e2=d3*y2-d3*y3+d4*y4-d5*y5+d5*y6-d4*y8
      cons=zero
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn7dy=top/bot
     !
      xcu=-d4*y3+d8*y4-d4*y5
      xsq=-two*y1+two*y2+two*y3-two*y5-two*y6+two*y7
      xon=two*y1+two*y3-d8*y4+two*y5+two*y7
      ecu=zero
      esq=-y1+two*y2-y3+y5-two*y6+y7
      eon=y1-y3+y5-y7
      x2e2=-d3*y1+d6*y2-d3*y3+d3*y5-d6*y6+d3*y7
      x1e1=-d6*y1+d8*y2-two*y3-two*y5+d8*y6-d6*y7
      x2e1=d5*y1-d8*y2+d3*y3-d3*y5+d8*y6-d5*y7
      x1e2=d4*y1-d8*y2+d4*y3+d4*y5-d8*y6+d4*y7
      cons=-two*y2+two*y6
     !
      top=xi3*xcu+xi2*xsq+xi*xon+et3*ecu+et2*esq+eta*eon
      top=top+xi2et2*x2e2+xi1et1*x1e1
      top=top+xi2et1*x2e1+xi1et2*x1e2+cons
     !
      dn8dy=top/bot
     !
      bee=zero
      bee(1,1)=dn1dx
      bee(1,3)=dn2dx
      bee(1,5)=dn3dx
      bee(1,7)=dn4dx
      bee(1,9)=dn5dx
      bee(1,11)=dn6dx
      bee(1,13)=dn7dx
      bee(1,15)=dn8dx
      bee(2,2)=dn1dy
      bee(2,4)=dn2dy
      bee(2,6)=dn3dy
      bee(2,8)=dn4dy
      bee(2,10)=dn5dy
      bee(2,12)=dn6dy
      bee(2,14)=dn7dy
      bee(2,16)=dn8dy
      bee(3,1)=dn1dy
      bee(3,3)=dn2dy
      bee(3,5)=dn3dy
      bee(3,7)=dn4dy
      bee(3,9)=dn5dy
      bee(3,11)=dn6dy
      bee(3,13)=dn7dy
      bee(3,15)=dn8dy
      bee(3,2)=dn1dx
      bee(3,4)=dn2dx
      bee(3,6)=dn3dx
      bee(3,8)=dn4dx
      bee(3,10)=dn5dx
      bee(3,12)=dn6dx
      bee(3,14)=dn7dx
      bee(3,16)=dn8dx
      RETURN
     END SUBROUTINE bee8

     SUBROUTINE fsparv(kv,km,g,kdiag)
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
     END SUBROUTINE fsparv

   SUBROUTINE fsparv_CSR(kv,km,g,ia,ja,nnz,ielacum,nnzaum)
     !
     ! This subroutine assles element matrices into a symmetric skyline
     ! global matrix.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::g(:),ielacum,nnzaum
      INTEGER,INTENT(INOUT)::nnz,ia(:),ja(:)
      REAL(iwp),INTENT(IN)::km(:,:)
      REAL(iwp),INTENT(IN OUT)::kv(:)
      INTEGER::i,idof,k,j,iw,ival,size,aval,m,diag,res,nzero
      INTEGER::gg(36)
      idof=UBOUND(g,1)
      size=UBOUND(km,1)
      res=0
      nzero=0
      aval=1
      diag=0
      DO i=1,20
          DO j=1+diag,idof
              IF(g(j)/=0.and.g(i)/=0)THEN
              IF(km(j,i).ne.abs(0.0_iwp))THEN
                nnz=nnz+1
                kv(nnz)=km(j,i)
                ja(nnz)=g(j)
                ia(nnz)=g(i)
              END IF
         END IF
          END DO
          diag=diag+1
          IF(diag>=20)diag=20
      END DO

      !iw=1
      !diag=0
      !DO i=1,20
      !   k=g(i)
      !   IF(k>0)THEN
      !    DO j=1+diag,idof
      !        m=g(j)
      !        IF(g(j)/=0)THEN
      !        IF(km(i,j).ne.abs(0.0_iwp))THEN
      !          nnz=nnz+1
      !          kv(aval+ielacum*nnzaum)=km(i,j)
      !          ja(aval+ielacum*nnzaum)=g(j)
      !          IF(iw==i)THEN
      !                !res=res+1
      !                ia(i-res+ielacum*nnzaum)=nnz
      !                iw=iw+1
      !          END IF
      !          aval=aval+1
      !        END IF
      !        END IF
      !    END DO
      !   END IF
      !    diag=diag+1
      !    IF(k==0)res=res+1
      !    IF(k==0)iw=iw+1
      !END DO
      !ia(i-res+ielacum*nnzaum)=nnz+1

     RETURN
     END SUBROUTINE fsparv_CSR

     SUBROUTINE sparin(kv,kdiag)
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
     END SUBROUTINE sparin

     SUBROUTINE spabac(kv,loads,kdiag)
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
     END SUBROUTINE spabac

     SUBROUTINE checon(DuB1,DuB2,UB1,UB2,converged,tol)
     !
     ! This subroutine sets converged to .FALSE. if relative change in loads
     ! and oldlds is greater than tol and updates oldlds.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::DuB1(0:),DuB2(0:),UB1(0:),UB2(0:)
      REAL(iwp),INTENT(IN)::tol
      !REAL(iwp),INTENT(IN OUT)::oldlds(0:)
      INTEGER::i
      REAL(iwp)::VAL,upvalB1,upvalB2,downval,energy_acum,energy_step,work
      LOGICAL,INTENT(OUT)::converged
      CONVERGED=.TRUE.
        upvalB1=0.0;upvalB2=0.0;downval=0.0;energy_acum=0.0;energy_step=0.0
      DO i=1,SIZE(DuB1)-1
          upvalB1=upvalB1+DuB1(i)**2
          !downval=downval+oldlds(i)**2
          !energy_acum=energy_acum+ABS(forces_acum(i)*oldlds(i))
          !energy_step=energy_step+ABS(forces(i)*loads(i))
      END DO

        DO i=1,SIZE(DuB2)-1
          upvalB2=upvalB2+DuB2(i)**2
          !downval=downval+oldlds(i)**2
          !energy_acum=energy_acum+ABS(forces_acum(i)*oldlds(i))
          !energy_step=energy_step+ABS(forces(i)*loads(i))
      END DO

      upvalB1=SQRT(upvalB1)
      upvalB2=SQRT(upvalB2)
      CONVERGED=(upvalB1+upvalB2<=tol)

      !work=energy_step/energy_acum
      !upval=SQRT(upval)
      !downval=SQRT(downval)
      !IF(upval/downval<=tol.and.work<=tol)THEN
      !CONVERGED=.true.
      !END IF
      !VAL=upval/downval
      !CONVERGED=(MAXVAL(ABS(loads-oldlds))/MAXVAL(ABS(loads))<=tol)
      !VAL=MAXVAL(ABS(loads-oldlds))/MAXVAL(ABS(loads))

      !oldlds=loads
     RETURN
     END SUBROUTINE checon

    SUBROUTINE checon_1(DuC,DuP,converged,tol,neq)
     !
     ! This subroutine sets converged to .FALSE. if relative change in loads
     ! and oldlds is greater than tol and updates oldlds.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::DuC(:),DuP(:)
      REAL(iwp),INTENT(IN)::tol
      INTEGER,INTENT(IN)::neq
      !REAL(iwp),INTENT(IN OUT)::oldlds(0:)
      INTEGER::i
      REAL(iwp)::VAL,upvalB1,upvalB2,downval,energy_acum,energy_step,work
      LOGICAL,INTENT(OUT)::converged
      converged=.FALSE.
        upvalB1=0.0;upvalB2=0.0;downval=0.0;energy_acum=0.0;energy_step=0.0
        !upvalB1=MAXVAL(DuC)
        !upvalB2=MAXVAL(DuP)
      DO i=1,neq
          upvalB1=upvalB1+DuC(i)**2
      END DO

    DO i=1,neq
          upvalB2=upvalB2+DuP(i)**2
      END DO

      upvalB1=SQRT(upvalB1)
      upvalB2=SQRT(upvalB2)

      VAL=(ABS(upvalB2))/(ABS(upvalB1))
      CONVERGED=(VAL<=tol)


     RETURN
     END SUBROUTINE checon_1

     SUBROUTINE checon2(gravlo,ddylds,bdylds,tol,converged)
     !
     ! This subroutine sets converged to .FALSE. if relative change in loads
     ! and oldlds is greater than tol and updates oldlds.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::tol,gravlo(0:),ddylds(0:),bdylds(0:)
      REAL(iwp)::FGEN,FEXT,zero
      LOGICAL,INTENT(OUT)::converged
      INTEGER:: i,neq
      CONVERGED=.false.
      zero=0.0_iwp
      neq=UBOUND(gravlo,1)
      FEXT=zero
      FGEN=zero
      DO i=1,neq
        FEXT=FEXT+gravlo(i)**2
        FGEN=FGEN+(gravlo(i)-(ddylds(i)+bdylds(i)))**2
      END DO
        FEXT=sqrt(FEXT)
        FGEN=sqrt(FGEN)
        IF((FEXT-FGEN)/FEXT<=tol)THEN
          CONVERGED=.true.
        END IF
      PRINT*,(FEXT-FGEN)/FEXT
     RETURN
     END SUBROUTINE checon2

     SUBROUTINE invar(stress,sigm,dsbar,theta)
     !
     ! This subroutine forms the stress invariants in 2- or 3-d. iii
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

     SUBROUTINE invar_MC(stress,sigm,dsbar,theta)
     !
     ! This subroutine forms the stress invariants in 2- or 3-d. iii
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
        dx=(two*sx-sy-sz)/three
        dy=(two*sy-sz-sx)/three
        dz=(two*sz-sx-sy)/three
        xj3=dx*dy*dz-dz*txy**2
        dsbar=SQRT(0.5_iwp*(dx**2+dy**2+dz**2)+txy**2)
        sine=-(xj3/dsbar**3)*three*SQRT(three)/two
        theta=ASIN(sine)/three
        !IF(dsbar<small)THEN
        !  theta=zero
        !ELSE
        !  dx=(two*sx-sy-sz)/three
        !  dy=(two*sy-sz-sx)/three
        !  dz=(two*sz-sx-sy)/three
        !  xj3=dx*dy*dz-dz*txy**2
        !  sine=-thpt5*xj3/dsbar**3
        !  IF(sine>=one)sine=one
        !  IF(sine<-one)sine=-one
        !  theta=ASIN(sine)/three
        !END IF
      CASE DEFAULT
        WRITE(*,*)"wrong size for nst in invar"
      END SELECT

     RETURN
     END SUBROUTINE invar_MC

     SUBROUTINE invar2(stress,dsbar,theta)
     !
     ! This subroutine forms the stress invariants in 2- or 3-d. iii
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::stress(:)
      REAL(iwp),INTENT(OUT),OPTIONAL::dsbar,theta
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
        dsbar=SQRT((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+six*txy**2)/SQRT(three)
        dsbar=dsbar*SQRT(three/two)
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
     END SUBROUTINE invar2

     SUBROUTINE mocouf(phi,c,sigm,dsbar,theta,f)
     !
     ! This subroutine calculates the value of the yield function
     ! for a Mohr-Coulomb material (phi in degrees).
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::phi,c,sigm,dsbar,theta
      REAL(iwp),INTENT(OUT)::f
      REAL(iwp)::phir,snph,csph,csth,snth,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,   &
        d180=180.0_iwp
      phir=phi*d4*ATAN(one)/d180
      snph=SIN(phir)
      csph=COS(phir)
      csth=COS(theta)
      snth=SIN(theta)
      f=snph*sigm+dsbar*(csth/SQRT(d3)-snth*snph/d3)-c*csph
     RETURN
     END SUBROUTINE mocouf

     SUBROUTINE mocouf_Sloan(phi,c,sigm,dsbar,theta,f,a_slo,psi,Tload,KT)
     !
     ! This subroutine calculates the value of the yield function
     ! for a Mohr-Coulomb material (phi in degrees).
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::phi,c,sigm,dsbar,theta,a_slo,psi,Tload
      REAL(iwp),INTENT(OUT)::f,KT
      REAL(iwp)::phir,snph,csph,csth,snth,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,   &
        d180=180.0_iwp,fprev,LOADT,LOADTRad,KTDeg,pi,p3,i3=3.0_iwp,            &
		FACTOR,NEWKT,ayield,Bval,Aval,theta_sign,cspsi,snpsi,zero=0.0_iwp,     &
        psirad,afact,cotpsi


    p3=one/3.0_iwp
    pi=ACOS(-one)
    LOADTRad=Tload*pi/d180
    NEWKT=zero
    KT=zero

	phir=phi*d4*ATAN(one)/d180 !-Friction angle in radians
	csth=COS(theta) !Theta is already in radians
	snth=SIN(theta)
    psirad=psi*pi/d180 !psi in radians
	snpsi=SIN(psirad)
	cspsi=COS(psirad)
	fprev=psirad*sigm+dsbar*(csth/SQRT(d3)-snth*snpsi/d3)-c*cspsi
	KTDeg=theta*d180/pi
    cotpsi=one/TAN(psirad)
    afact=a_slo*c*cotpsi

    IF(ABS(theta)<=LOADTRad)THEN
       KT=(csth/SQRT(d3)-snpsi*snth/d3)
    ELSE
      IF(theta>=0.0)theta_sign=one
      IF(theta<0.0)theta_sign=-one
      Aval=p3*COS(LOADTRad)*(i3+TAN(LOADTRad)*TAN(i3*LOADTRad)+    &
           one/SQRT(i3)*theta_sign*(TAN(i3*LOADTRad)-3*TAN(LOADTRad))*snpsi)
      Bval=(one/(i3*COS(i3*LOADTRad)))*(theta_sign*SIN(LOADTRad)+ &
            (one/SQRT(i3))*snpsi*COS(LOADTRad))
      KT=(Aval-Bval*SIN(i3*theta))*(1.0/SQRT(3.0))

    END IF
    FACTOR=(dsbar**2)*(KT**2)+(afact**2)*(snpsi**2)
    f=sigm*snpsi+SQRT(FACTOR)-c*cspsi
   RETURN
   END SUBROUTINE mocouf_Sloan

     SUBROUTINE mocouq(psi,dsbar,theta,dq1,dq2,dq3)
     !
     ! This subroutine forms the derivatives of a Mohr-Coulomb potential
     ! function with respect to the three invariants (psi in degrees).
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::psi,dsbar,theta
      REAL(iwp),INTENT(OUT)::dq1,dq2,dq3
      REAL(iwp)::psir,snth,snps,sq3,c1,csth,cs3th,tn3th,tnth,zero=0.0_iwp,     &
        pt49=0.49_iwp,pt5=0.5_iwp,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,           &
        d180=180.0_iwp
      psir=psi*d4*ATAN(one)/d180
      snth=SIN(theta)
      snps=SIN(psir)
      sq3=SQRT(d3)
      dq1=snps
      if(ABS(snth).GT.pt49)THEN
        c1=one
        IF(snth.LT.zero)c1=-one
        dq2=(sq3*pt5-c1*snps*pt5/sq3)*sq3*pt5/dsbar
        dq3=zero
      ELSE
        csth=COS(theta)
        cs3th=COS(d3*theta)
        tn3th=TAN(d3*theta)
        tnth=snth/csth
        dq2=sq3*csth/dsbar*((one+tnth*tn3th)+snps*(tn3th-tnth)/sq3)*pt5
        dq3=pt5*d3*(sq3*snth+snps*csth)/(cs3th*dsbar*dsbar)
      END IF
     RETURN
     END SUBROUTINE mocouq

     SUBROUTINE formm(stress,m1,m2,m3)
     !
     ! This subroutine forms the derivatives of the invariants with respect to
     ! stress in 2- or 3-d.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::stress(:)
      REAL(iwp),INTENT(OUT)::m1(:,:),m2(:,:),m3(:,:)
      REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
        two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
      INTEGER::nst,i,j
      nst=UBOUND(stress,1)
      SELECT CASE(nst)
      CASE(4)
        sx=stress(1)
        sy=stress(2)
        txy=stress(3)
        sz=stress(4)
        dx=(two*sx-sy-sz)/three
        dy=(two*sy-sz-sx)/three
        dz=(two*sz-sx-sy)/three
        sigm=(sx+sy+sz)/three
        m1=zero
        m2=zero
        m3=zero
        m1(1,1:2)=one
        m1(2,1:2)=one
        m1(4,1:2)=one
        m1(1,4)=one
        m1(4,4)=one
        m1(2,4)=one
        IF(sigm==0.0)THEN
          m1=0.0
        ELSE
        m1=m1/nine/sigm
        END IF
        m2(1,1)=two/three
        m2(2,2)=two/three
        m2(4,4)= two/three
        m2(2,4)=-one/three
        m2(4,2)=-one/three
        m2(1,2)=-one/three
        m2(2,1)=-one/three
        m2(1,4)=-one/three
        m2(4,1)=-one/three
        m2(3,3)=two
        m3(3,3)=-dz
        m3(1:2,3)=txy/three
        m3(3,1:2)=txy/three
        m3(3,4)=-two*txy/three
        m3(4,3)=-two*txy/three
        m3(1,1)=dx/three
        m3(2,4)=dx/three
        m3(4,2)=dx/three
        m3(2,2)=dy/three
        m3(1,4)=dy/three
        m3(4,1)=dy/three
        m3(4,4)=dz/three
        m3(1,2)=dz/three
        m3(2,1)=dz/three
      CASE(6)
        sx=stress(1)
        sy=stress(2)
        sz=stress(3)
        txy=stress(4)
        tyz=stress(5)
        tzx=stress(6)
        sigm=(sx+sy+sz)/three
        dx=sx-sigm
        dy=sy-sigm
        dz=sz-sigm
        m1=zero
        m2=zero
        m1(1:3,1:3)=one/(three*sigm)
        DO i=1,3
          m2(i,i)=two
          m2(i+3,i+3)=six
        END DO
        m2(1,2)=-one
        m2(1,3)=-one
        m2(2,3)=-one
        m3(1,1)=dx
        m3(1,2)=dz
        m3(1,3)=dy
        m3(1,4)=txy
        m3(1,5)=-two*tyz
        m3(1,6)=tzx
        m3(2,2)=dy
        m3(2,3)=dx
        m3(2,4)=txy
        m3(2,5)=tyz
        m3(2,6)=-two*tzx
        m3(3,3)=dz
        m3(3,4)=-two*txy
        m3(3,5)=tyz
        m3(3,6)=tzx
        m3(4,4)=-three*dz
        m3(4,5)=three*tzx
        m3(4,6)=three*tyz
        m3(5,5)=-three*dx
        m3(5,6)=three*txy
        m3(6,6)=-three*dy
        DO i=1,6
          DO j=i+1,6
            m1(j,i)=m1(i,j)
            m2(j,i)=m2(i,j)
            m3(j,i)=m3(i,j)
          END DO
        END DO
        m1=m1/three
        m2=m2/three
        m3=m3/three
      CASE DEFAULT
        WRITE(*,*)"nst size not recognised in formm"
      END SELECT
     RETURN
     END SUBROUTINE formm

     SUBROUTINE formm_Sloan(stress,dsbar,m1,m2,m3)
     !
     ! This subroutine forms the derivatives of the invariants with respect to
     ! stress in 2- or 3-d.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::stress(:),dsbar
      REAL(iwp),INTENT(OUT)::m1(:),m2(:),m3(:)
      REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
        two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
      INTEGER::nst,i,j
      nst=UBOUND(stress,1)
      SELECT CASE(nst)
      CASE(4)
        sx=stress(1)
        sy=stress(2)
        txy=stress(3)
        sz=stress(4)
        dx=(two*sx-sy-sz)/three
        dy=(two*sy-sz-sx)/three
        dz=(two*sz-sx-sy)/three
        sigm=(sx+sy+sz)/three
        m1=zero
        m2=zero
        m3=zero
        m1(1)=one/three
        m1(2)=one/three
        m1(3)=zero
        m1(4)=one/three

        m2(1)=dx
        m2(2)=dy
        m2(3)=two*txy
        m2(4)=dz
       m2=m2/(two*dsbar)

        m3(1)=dz*dy-txy**2+dsbar**2.0_iwp/three
        m3(2)=dz*dx-txy**2+dsbar**2.0_iwp/three
        m3(3)=two*(txy**2.0_iwp-dz*txy)
        m3(4)=dx*dy-txy**2+dsbar**2.0_iwp/three




      CASE DEFAULT
        WRITE(*,*)"nst size not recognised in formm"
      END SELECT
     RETURN
     END SUBROUTINE formm_Sloan

     SUBROUTINE dismsh(loads,nf,ratmax,g_coord,g_num,argv,nlen,ips)
     !
     ! This subroutine produces a PostScript output file "*.dis" displaying
     ! the deformed finite element mesh.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::g_coord(:,:),loads(0:),ratmax
      INTEGER,INTENT(IN)::g_num(:,:),ips,nf(:,:),nlen
      CHARACTER(*),INTENT(IN)::argv
      REAL(iwp)::width,height,scale=72,sxy,xo,yo,x,y,dismag,vmax
      REAL(iwp)::xmin,xmax,ymin,ymax,dmax,zero=0.0_iwp,pt5=0.5_iwp,            &
        opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,ept5=8.5_iwp,d11=11.0_iwp
      INTEGER::i,ii,j,jj,nn,nel,nod
      OPEN(ips,FILE=argv(1:nlen)//'.dis')
     !
      nn=UBOUND(nf,2)
      xmin=g_coord(1,1)
      xmax=g_coord(1,1)
      ymin=g_coord(2,1)
      ymax=g_coord(2,1)
      DO i=2,nn
        IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
        IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
        IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
        IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
      END DO
      width=xmax-xmin
      height=ymax-ymin
      dmax=ratmax*width
      IF(height>width)dmax=ratmax*height
     !
      vmax=zero
      DO i=1,nn
        DO j=1,2
          IF(ABS(loads(nf(j,i)))>vmax)vmax=ABS(loads(nf(j,i)))
        END DO
      END DO
      dismag=dmax/vmax
     !
      xmin=g_coord(1,1)
      xmax=g_coord(1,1)
      ymin=g_coord(2,1)
      ymax=g_coord(2,1)
      DO i=1,nn
        IF(g_coord(1,i)+dismag*loads(nf(1,i))<xmin)                            &
          xmin=g_coord(1,i)+dismag*loads(nf(1,i))
        IF(g_coord(1,i)+dismag*loads(nf(1,i))>xmax)                            &
          xmax=g_coord(1,i)+dismag*loads(nf(1,i))
        IF(g_coord(2,i)+dismag*loads(nf(2,i))<ymin)                            &
          ymin=g_coord(2,i)+dismag*loads(nf(2,i))
        IF(g_coord(2,i)+dismag*loads(nf(2,i))>ymax)                            &
          ymax=g_coord(2,i)+dismag*loads(nf(2,i))
     !
        IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
        IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
        IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
        IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
      END DO
     !
      width =xmax-xmin
      height=ymax-ymin
     !
     !                       allow 1.5" margin minimum on each side of figure
     !
     !                       portrait mode
     !
      IF(height.GE.d11/ept5*width)THEN
     !
     !                       height governs the scale
     !
        sxy=scale*d8/height
        xo=scale*pt5*(ept5-d8*width/height)
        yo=scale*opt5
      ELSE
     !
     !                       width governs the scale
     !
        sxy=scale*fpt5/width
        xo=scale*opt5
        yo=scale*pt5*(d11-fpt5*height/width)
      END IF
     !
      nel=UBOUND(g_num,2)
      nod=UBOUND(g_num,1)
     !
     !                       start PostScript output
     !
      WRITE(ips,'(a)')'%!PS-Adobe-1.0'
      WRITE(ips,'(a)')'%%DocumentFonts: none'
      WRITE(ips,'(a)')'%%Pages: 1'
      WRITE(ips,'(a)')'%%EndComments'
      WRITE(ips,'(a)')'/m {moveto} def'
      WRITE(ips,'(a)')'/l {lineto} def'
      WRITE(ips,'(a)')'/s {stroke} def'
      WRITE(ips,'(a)')'/c {closepath} def'
      WRITE(ips,'(a)')'%%EndProlog'
      WRITE(ips,'(a)')'%%Page: 0 1'
      WRITE(ips,'(a)')'gsave'
     !
     !                       draw the deformed mesh
     !
      WRITE(ips,'(2f9.2,a)') xo, yo, ' translate'
      WRITE(ips,'(f9.2,a)') 0.5, ' setlinewidth'
      IF(nod==5)nod=4
      IF(nod==9)nod=8
      IF(nod==10)nod=9
      IF(nod==15)nod=12
      DO i=1,nel
        ii=g_num(1,i)
        IF(ii==0)CYCLE
        x=sxy*(g_coord(1,ii)+dismag*loads(nf(1,ii))-xmin)
        y=sxy*(g_coord(2,ii)+dismag*loads(nf(2,ii))-ymin)
        WRITE(ips,'(2f9.2,a)') x, y,' m'
        DO j=2,nod
          jj=g_num(j,i)
          x=sxy*(g_coord(1,jj)+dismag*loads(nf(1,jj))-xmin)
          y=sxy*(g_coord(2,jj)+dismag*loads(nf(2,jj))-ymin)
          WRITE(ips,'(2f9.2,a)') x, y,' l'
        END DO
        WRITE(ips,'(a)')'c s'
      END DO
     !
      WRITE(ips,'(a)')'grestore'
      WRITE(ips,'(a)')'showpage'
      CLOSE(ips)
     !
     RETURN
     END SUBROUTINE dismsh



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
          
SUBROUTINE ecmat_mf(ecm,fun,ndof,nodof,side)
!
! This subroutine forms the element consistent mass matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::fun(:)
 REAL(iwp),INTENT(OUT)::ecm(:,:)
 INTEGER,INTENT(IN)::nodof,ndof,side
 INTEGER::nod,i,j
 REAL,ALLOCATABLE::nt(:,:),tn(:,:)
 REAL::zero=0.0_iwp
 ecm=zero
 nod=ndof/nodof
 ALLOCATE(nt(ndof,nodof),tn(nodof,6))
 nt=zero
 tn=zero
 DO i=1,nod
   DO j=1,nodof
     nt((i-1)*nodof+j,j)=fun(i)
     !tn(j,(i-1)*nodof+j)=fun(i)
   END DO
 END DO

 IF(side==1)THEN
    tn(:,1:2)= TRANSPOSE(nt(13:14,:))
    tn(:,3:4)= TRANSPOSE(nt(11:12,:))
    tn(:,5:6)= TRANSPOSE(nt(9:10,:))
 ELSE
     tn=TRANSPOSE(nt(1:6,:))
 END IF
 
 
 ecm=MATMUL(nt,tn)
RETURN
END SUBROUTINE ecmat_mf
          
          
 SUBROUTINE ecmat_vis(ecm,fun,ndof,nodof,cc)
!
! This subroutine forms the element consistent mass matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::fun(:),cc(:,:)
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
 tn=MATMUL(cc,tn)
 ecm=MATMUL(nt,tn)
RETURN
END SUBROUTINE ecmat_vis

     SUBROUTINE vecmsh(loads,nf,ratmax,cutoff,g_coord,g_num,argv,nlen,ips)
     !
     ! This subroutine produces a PostScript output file "*.vec" displaying
     ! the nodal displacement vectors.
     !
     IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::g_coord(:,:),loads(0:),ratmax,cutoff
      INTEGER,INTENT(IN)::g_num(:,:),ips,nf(:,:),nlen
      REAL(iwp)::width,height,scale=72,sxy,xo,yo,x1,y1,x2,y2,dismag,           &
        zero=0.0_iwp,pt5=0.5_iwp,opt5=1.5_iwp,fpt5=5.5_iwp,d8=8.0_iwp,         &
        ept5=8.5_iwp,d11=11.0_iwp,xmin,xmax,ymin,ymax,dmax,vlen,vmax
      INTEGER::i,j,k,l,nn,nels,nod,ns,i1,i2,j1,j2
      INTEGER,ALLOCATABLE::corner(:,:)
      CHARACTER(*),INTENT(IN)::argv
      LOGICAL::draw
     !                       formats
      OPEN(ips,FILE=argv(1:nlen)//'.vec')
     !                       open output file and compute scale factors
      nn=UBOUND(nf,2)
     !
      xmin=g_coord(1,1)
      xmax=g_coord(1,1)
      ymin=g_coord(2,1)
      ymax=g_coord(2,1)
      DO i=2,nn
        IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
        IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
        IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
        IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
      END DO
      width =xmax-xmin
      height=ymax-ymin
      dmax=ratmax*width
      IF(height>width)dmax=ratmax*height
     !
      vmax=zero
      DO i=1,nn
        DO j=1,2
          IF(ABS(loads(nf(j,i)))>vmax)vmax=ABS(loads(nf(j,i)))
        END DO
      END DO
      dismag=dmax/vmax
     !
      xmin=g_coord(1,1)
      xmax=g_coord(1,1)
      ymin=g_coord(2,1)
      ymax=g_coord(2,1)
     !
      DO i=1,nn
        IF(g_coord(1,i)+dismag*loads(nf(1,i)) < xmin)                          &
          xmin=g_coord(1,i)+dismag*loads(nf(1,i))
        IF(g_coord(1,i)+dismag*loads(nf(1,i)) > xmax)                          &
          xmax=g_coord(1,i)+dismag*loads(nf(1,i))
        IF(g_coord(2,i)+dismag*loads(nf(2,i)) < ymin)                          &
          ymin=g_coord(2,i)+dismag*loads(nf(2,i))
        IF(g_coord(2,i)+dismag*loads(nf(2,i)) > ymax)                          &
          ymax=g_coord(2,i)+dismag*loads(nf(2,i))
     !
        IF(g_coord(1,i)<xmin)xmin=g_coord(1,i)
        IF(g_coord(1,i)>xmax)xmax=g_coord(1,i)
        IF(g_coord(2,i)<ymin)ymin=g_coord(2,i)
        IF(g_coord(2,i)>ymax)ymax=g_coord(2,i)
      END DO
     !
      width=xmax-xmin
      height=ymax-ymin
     !
     !                       allow 1.5" margin minimum on each side of figure
     !
     !                       Portrait mode
     !
      IF(height.GE.d11/ept5*width)THEN
     !
     !                       height governs the scale
     !
        sxy=scale*d8/height
        xo=scale*pt5*(ept5-d8*width/height)
        yo=scale*opt5
      ELSE
     !
     !                       width governs the scale
     !
        sxy=scale*fpt5/width
        xo=scale*opt5
        yo=scale*pt5*(d11-fpt5*height/width)
      END IF
     !
      WRITE(ips,'(a)')'%!PS-Adobe-1.0'
      WRITE(ips,'(a)')'%%DocumentFonts: none'
      WRITE(ips,'(a)')'%%Pages: 1'
      WRITE(ips,'(a)')'%%EndComments'
      WRITE(ips,'(a)')'/m {moveto} def'
      WRITE(ips,'(a)')'/l {lineto} def'
      WRITE(ips,'(a)')'/s {stroke} def'
      WRITE(ips,'(a)')'/c {closepath} def'
      WRITE(ips,'(a)')'/edef {exch def} bind def'
      WRITE(ips,'(a)')                                                         &
        '/arrow {/@to_y edef /@to_x edef /@from_y edef /@from_x edef'
      WRITE(ips,'(a)')'/@dx @to_x @from_x sub def /@dy @to_y @from_y sub def'
      WRITE(ips,'(a)')'/@length @dx @dx mul @dy @dy mul add sqrt def'
      WRITE(ips,'(a)')'/@angle @dy @dx atan def'
      WRITE(ips,'(a)')'gsave @from_x @from_y translate @angle rotate'
      WRITE(ips,'(a)')                                                         &
        '0 0 moveto @length 0 lineto currentpoint stroke newpath moveto'
      WRITE(ips,'(a)')'-4 -2 rlineto @length 0 moveto'
      WRITE(ips,'(a)')'-4  2 rlineto stroke grestore'
      WRITE(ips,'(a)')'} def'
      WRITE(ips,'(a)')'/*sf {'
      WRITE(ips,'(a)')'exch findfont exch'
      WRITE(ips,'(a)')                                                         &
        'dup type /arraytype eq {makefont}{scalefont} ifelse setfont'
      WRITE(ips,'(a)')'} bind def'
      WRITE(ips,'(a)')'/languagelevel where'
      WRITE(ips,'(a)')'{pop languagelevel} {1} ifelse'
      WRITE(ips,'(a)')'2 lt { % ifelse'
      WRITE(ips,'(a)')'/sf /*sf load def'
      WRITE(ips,'(a)')'} { % else'
      WRITE(ips,'(a)')'/sf /selectfont load def'
      WRITE(ips,'(a)')'} ifelse'
      WRITE(ips,'(a)')'%%EndProlog'
      WRITE(ips,'(a)')'%%Page: 0 1'
      WRITE(ips,'(a)')'gsave'
     !
      WRITE(ips,'(2f9.2,a)')xo,yo, ' translate'
      WRITE(ips,'(f9.2,a)')0.5,' setlinewidth'
     !
     !                       draw the displacement vectors
     !
      vmax=zero
      DO i=1,nn
        vlen=loads(nf(1,i))**2+loads(nf(2,i))**2
        IF(vlen>vmax)vmax=vlen
      END DO
      vmax=SQRT(vmax)*cutoff
      DO i=1,nn
        vlen=SQRT(loads(nf(1,i))**2+loads(nf(2,i))**2)
        x1=sxy*(g_coord(1,i)-xmin)
        y1=sxy*(g_coord(2,i)-ymin)
        x2=sxy*(g_coord(1,i)+dismag*loads(nf(1,i))-xmin)
        y2=sxy*(g_coord(2,i)+dismag*loads(nf(2,i))-ymin)
        IF(vlen>vmax)THEN
          WRITE(ips,'(2f9.2,a,2f9.2,a)') x1, y1,' ', x2, y2, ' arrow'
          WRITE(ips,'(a)') 's'
        END IF
      END DO
     !
     !                       draw the mesh border
     !
      nels=UBOUND(g_num,2)
      nod=UBOUND(g_num,1)
      IF(nod==3.OR.nod==6.OR.nod==10.OR.nod==15)ns=3
      IF(nod==4.OR.nod==5.OR.nod==8.OR.nod==9)ns=4
      ALLOCATE(corner(ns,2))
      IF(nod== 3)corner=RESHAPE((/1,2,3,2,3,1/),(/3,2/))
      IF(nod== 6)corner=RESHAPE((/1,3,5,3,5,1/),(/3,2/))
      IF(nod==10)corner=RESHAPE((/1,4,7,4,7,1/),(/3,2/))
      IF(nod==15)corner=RESHAPE((/1,5,9,5,9,1/),(/3,2/))
      IF(nod== 4)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
      IF(nod== 5)corner=RESHAPE((/1,2,3,4,2,3,4,1/),(/4,2/))
      IF(nod== 8)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
      IF(nod== 9)corner=RESHAPE((/1,3,5,7,3,5,7,1/),(/4,2/))
      DO i=1,nels
        DO j=1,ns
          draw=.TRUE.
          i1=g_num(corner(j,1),i)
          i2=g_num(corner(j,2),i)
          DO k=1,nels
            DO l=1,ns
              j1=g_num(corner(l,1),k)
              j2=g_num(corner(l,2),k)
              IF((i1==j2).AND.(i2==j1))THEN
                draw=.FALSE.
                EXIT
              END IF
            END DO
            IF(.NOT.draw)EXIT
          END DO
          IF(draw)THEN
            x1=sxy*(g_coord(1,i1)-xmin)
            y1=sxy*(g_coord(2,i1)-ymin)
            WRITE(ips,'(2f9.2,a)')x1, y1,' m'
            x1=sxy*(g_coord(1,i2)-xmin)
            y1=sxy*(g_coord(2,i2)-ymin)
            WRITE(ips,'(2f9.2,a)')x1, y1,' l'
            WRITE(ips,'(a)')' s'
          END IF
        END DO
      END DO
     !                       close output file?
      WRITE(ips,'(a)')'grestore'
      WRITE(ips,'(a)')'showpage'
      CLOSE(ips)
     RETURN
     END SUBROUTINE vecmsh

     SUBROUTINE shape_der(der,points,i)
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
      REAL,PARAMETER::zero=0.0_iwp,pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,  &
        pt75=0.75_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,&
        d6=6.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d10=10.0_iwp,d11=11.0_iwp,            &
        d12=12.0_iwp,d16=16.0_iwp,d18=18.0_iwp,d27=27.0_iwp,d32=32.0_iwp,      &
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
          der(1,1)= (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*          &
            (eta-one)*(zeta-one)/d8
          der(1,2)=-(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*          &
            (eta-one)*(zeta+one)/d8
          der(1,3)=-(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*          &
            (eta-one)*(zeta+one)/d8
          der(1,4)= (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*          &
            (eta-one)*(zeta-one)/d8
          der(1,5)= -(eta-one)*(zeta+one)*(zeta-one)*xi
          der(1,6)=-(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
          der(1,7)=  (eta+one)*(eta-one)*(zeta+one)*xi
          der(1,8)= (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
          der(1,9)= -(eta+one)*(eta-one)*(zeta-one)*xi
          der(1,10)= (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*         &
            (eta+one)*(zeta-one)/d8
          der(1,11)=-(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*         &
            (eta+one)*(zeta+one)/d8
          der(1,12)=-(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*         &
            (eta+one)*(zeta+one)/d8
          der(1,13)= (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*         &
            (eta+one)*(zeta-one)/d8
          der(1,14)=  (eta+one)*(zeta+one)*(zeta-one)*xi
          der(2,1)= (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*          &
            (xi-one)*(zeta-one)/d8
          der(2,2)=-(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*          &
            (xi-one)*(zeta+one)/d8
          der(2,3)=-(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*          &
            (xi+one)*(zeta+one)/d8
          der(2,4)= (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*          &
            (xi+one)*(zeta-one)/d8
          der(2,5)=-(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
          der(2,6)= -(xi-one)*(zeta+one)*(zeta-one)*eta
          der(2,7)=  (xi+one)*(xi-one)*(zeta+one)*eta
          der(2,8)=  (xi+one)*(zeta+one)*(zeta-one)*eta
          der(2,9)= -(xi+one)*(xi-one)*(zeta-one)*eta
          der(2,10)= (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*         &
            (xi-one)*(zeta-one)/d8
          der(2,11)=-(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*         &
            (xi-one)*(zeta+one)/d8
          der(2,12)=-(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*         &
            (xi+one)*(zeta+one)/d8
          der(2,13)= (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*         &
            (xi+one)*(zeta-one)/d8
          der(2,14)= (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
          der(3,1)= (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*          &
            (xi-one)*(eta-one)/d8
          der(3,2)=-(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*          &
            (xi-one)*(eta-one)/d8
          der(3,3)=-(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*          &
            (xi+one)*(eta-one)/d8
          der(3,4)= (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*          &
            (xi+one)*(eta-one)/d8
          der(3,5)= -(xi+one)*(xi-one)*(eta-one)*zeta
          der(3,6)= -(xi-one)*(eta+one)*(eta-one)*zeta
          der(3,7)= (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
          der(3,8)=  (xi+one)*(eta+one)*(eta-one)*zeta
          der(3,9)=-(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
          der(3,10)= (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*         &
            (xi-one)*(eta+one)/d8
          der(3,11)=-(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*         &
            (xi-one)*(eta+one)/d8
          der(3,12)=-(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*         &
            (xi+one)*(eta+one)/d8
          der(3,13)= (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*         &
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
              der(1,l)=pt125*xii(l)*(one+eta0)*(one+zeta0)*                    &
                (two*xi0+eta0+zeta0-one)
              der(2,l)=pt125*etai(l)*(one+xi0)*(one+zeta0)*                    &
                (xi0+two*eta0+zeta0-one)
              der(3,l)=pt125*zetai(l)*(one+xi0)*(one+eta0)*                    &
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
     END SUBROUTINE shape_der

     SUBROUTINE shape_der2(der,points,i)
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
      REAL,PARAMETER::zero=0.0_iwp,pt125=0.125_iwp,pt25=0.25_iwp,pt5=0.5_iwp,  &
        pt75=0.75_iwp,one=1.0_iwp,two=2.0_iwp,d3=3.0_iwp,d4=4.0_iwp,d5=5.0_iwp,&
        d6=6.0_iwp,d8=8.0_iwp,d9=9.0_iwp,d10=10.0_iwp,d11=11.0_iwp,            &
        d12=12.0_iwp,d16=16.0_iwp,d18=18.0_iwp,d27=27.0_iwp,d32=32.0_iwp,      &
        d36=36.0_iwp,d54=54.0_iwp,d64=64.0_iwp,d128=128.0_iwp
      INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
      ndim=UBOUND(der,1)
      nod= UBOUND(der,2)
      SELECT CASE(ndim)
      CASE(1)   ! one dimensional elements

      CASE(2)      ! two dimensional elements
        xi=points(i,1)
        eta=points(i,2)
        c1=xi
        c2=eta
        c3=one-c1-c2
        etam=(one-eta)
        etap=(one+eta)
        xim= (one-xi)
        xip= (one+xi)
        x2p1=two*xi+one
        x2m1=two*xi-one
        e2p1=two*eta+one
        e2m1=two*eta-one
        SELECT CASE(nod)

        CASE (4)

          der(1,1)=-0.25_iwp
          der(1,2)=-0.25_iwp
          der(1,3)=0.25_iwp
          der(1,4)=0.25_iwp
          der(2,1)=-0.25_iwp
          der(2,2)=0.25_iwp
          der(2,3)=0.25_iwp
          der(2,4)=-0.25_iwp

        END SELECT
      END SELECT

     RETURN
     END SUBROUTINE shape_der2

     FUNCTION determinant(jac)RESULT(det)
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
      SELECT CASE(it)
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
     END FUNCTION determinant

     SUBROUTINE invert(matrix)
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
     END SUBROUTINE invert

     SUBROUTINE beemat(bee,deriv)
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
      SELECT CASE (ih)
      CASE(3,4)
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
      CASE(6)
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
        WRITE(*,*)'wrong dimension for nst in bee matrix'
      END SELECT
     RETURN
     END SUBROUTINE beemat

     SUBROUTINE vonmisf(c,dsbar,f)
     !
     ! This subroutine calculates the value of the yield function
     ! for a Mohr-Coulomb material (phi in degrees).
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::c,dsbar
      REAL(iwp),INTENT(OUT)::f
      REAL(iwp)::d3=3.0_iwp,three,sx,sy,txy,sz,six=6.0_iwp,zero=0.0_iwp

        !sx=stress(1)
        !sy=stress(2)
        !txy=stress(3)
        !sz=stress(4)
        three=sqrt(3.0_iwp)
      f=dsbar-(three*c)
     RETURN
     END SUBROUTINE vonmisf


    SUBROUTINE vonmisq(stress,dq1,dq2,dq3)
     !
     ! This subroutine forms the derivatives of a Mohr-Coulomb potential
     ! function with respect to the three riants (psi in degrees).
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::stress(:)
      REAL(iwp),INTENT(OUT)::dq1,dq2,dq3
      REAL(iwp)::three,sx,sy,txy,sz,six=6.0_iwp,zero=0.0_iwp,two,t

        sx=stress(1)
        sy=stress(2)
        txy=stress(3)
        sz=stress(4)
        three=sqrt(3.0_iwp)
        two=sqrt(2.0_iwp)
        t=(1.0_iwp/three)*(((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+six*txy**2)**0.5)

      dq1=zero
      dq2=(three/two)*(1.0_iwp/t)
      dq3=zero

     RETURN
     END SUBROUTINE vonmisq


      SUBROUTINE formmvm(stress,m1,m2,m3)
     !
     ! This subroutine forms the derivatives of the invariants with respect to
     ! stress in 2- or 3-d.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::stress(:)
      REAL(iwp),INTENT(OUT)::m1(:,:),m2(:,:),m3(:,:)
      REAL(iwp)::sx,sy,txy,tyz,tzx,sz,dx,dy,dz,sigm,zero=0.0_iwp,one=1.0_iwp,  &
        two=2.0_iwp,three=3.0_iwp,six=6.0_iwp,nine=9.0_iwp
      INTEGER::nst,i,j
      nst=UBOUND(stress,1)
      SELECT CASE(nst)
      CASE(4)
        sx=stress(1)
        sy=stress(2)
        txy=stress(3)
        sz=stress(4)
        dx=(two*sx-sy-sz)/three
        dy=(two*sy-sz-sx)/three
        dz=(two*sz-sx-sy)/three
        sigm=(sx+sy+sz)/three
        m1=zero
        m2=zero
        m3=zero
        m2(1,1)=two/three
        m2(2,2)=two/three
        m2(4,4)= two/three
        m2(2,4)=-one/three
        m2(4,2)=-one/three
        m2(1,2)=-one/three
        m2(2,1)=-one/three
        m2(1,4)=-one/three
        m2(4,1)=-one/three
        m2(3,3)=two
      CASE DEFAULT
        WRITE(*,*)"nst size not recognised in formm"
      END SELECT
     RETURN
      END SUBROUTINE formmvm

SUBROUTINE beam_mm_grav(mm,dens,size,nod,fun)
!
! This subroutine forms the consistent mass matrix of a beam element.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::nod
 REAL(iwp),INTENT(IN)::dens,size,fun(:,:)
 REAL(iwp),INTENT(OUT)::mm(:,:)
 REAL(iwp),ALLOCATABLE::nt(:,:),tn(:,:),mm_L(:,:)
 REAL(iwp)::two=2.0_iwp,zero=0.0_iwp,pt5=0.5_iwp,pt3=1.0_iwp/3_iwp,one=1.0_iwp
 INTEGER::i,m


 mm=zero
 IF(nod==2)THEN
 mm(1,1)=pt5
 mm(2,2)=pt5
 mm(3,3)=pt5
 mm(4,4)=pt5
 mm=mm*dens*size !the thickness and width of the element are considered as one
 ELSE
  ALLOCATE(nt(nod*2,2),tn(2,nod*2),mm_L(nod*2,nod*2))
  nt=zero;tn=zero;mm_L=zero

  DO i=1,nod
      DO m=1,2
         nt((i-1)*2+m,m)=fun(i,1)
      END DO
  END DO
  tn=TRANSPOSE(nt)
  mm=MATMUL(nt,tn)

   DO i=1,nod*2
    DO m=1,nod*2
      mm_L(i,i)=mm_L(i,i)+mm(i,m)
    END DO
   END DO
   mm=mm_L
  mm=mm*dens !the thickness and width of the element are considered as one
 END IF

RETURN
END SUBROUTINE beam_mm_grav

SUBROUTINE beam_mm_2(mm,dens,size,nod,fun)
!
! This subroutine forms the consistent mass matrix of a beam element.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::nod
 REAL(iwp),INTENT(IN)::dens,size,fun(:,:)
 REAL(iwp),INTENT(OUT)::mm(:,:)
 REAL(iwp),ALLOCATABLE::nt(:,:),tn(:,:),mm_L(:,:)
 REAL(iwp)::two=2.0_iwp,zero=0.0_iwp,pt5=0.5_iwp,pt3=1.0_iwp/3_iwp,one=1.0_iwp
 INTEGER::i,m


 mm=zero
 IF(nod==2)THEN
 mm(1,1)=pt5
 mm(2,2)=pt5
 mm(3,3)=pt5
 mm(4,4)=pt5
 mm=mm*dens*size !the thickness and width of the element are considered as one
 ELSE
  ALLOCATE(nt(nod*2,2),tn(2,nod*2),mm_L(nod*2,nod*2))
  nt=zero;tn=zero;mm_L=zero

  DO i=1,nod
      DO m=1,2
         nt((i-1)*2+m,m)=fun(i,1)
      END DO
  END DO
  tn=TRANSPOSE(nt)
  mm=MATMUL(nt,tn)*dens

   DO i=1,nod*2
    DO m=1,nod*2
      mm_L(i,i)=mm_L(i,i)+mm(i,m)
    END DO
   END DO
   mm=mm_L
   CONTINUE
  !mm=mm*dens !the thickness and width of the element are considered as one
 END IF

RETURN
END SUBROUTINE beam_mm_2

SUBROUTINE beam_mm_mf(mm,dens,size,nod,fun_1D,fun)
!
! This subroutine forms the consistent mass matrix of a beam element.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::nod
 REAL(iwp),INTENT(IN)::dens,size,fun_1D(:,:),fun(:)
 REAL(iwp),INTENT(OUT)::mm(:,:)
 REAL(iwp),ALLOCATABLE::nt(:,:),tn(:,:),mm_L(:,:)
 REAL(iwp)::two=2.0_iwp,zero=0.0_iwp,pt5=0.5_iwp,pt3=1.0_iwp/3_iwp,one=1.0_iwp
 INTEGER::i,m,nod_2D


 mm=zero
 IF(nod==2)THEN
 mm(1,1)=pt5
 mm(2,2)=pt5
 mm(3,3)=pt5
 mm(4,4)=pt5
 mm=mm*dens*size !the thickness and width of the element are considered as one
 ELSE
  nod_2D=UBOUND(fun,1)*2
  ALLOCATE(nt(nod*2,2),tn(2,nod_2D),mm_L(nod*2,nod*2))
  nt=zero;tn=zero;mm_L=zero

  DO i=1,nod
      DO m=1,2
         nt((i-1)*2+m,m)=fun_1D(i,1)
      END DO
  END DO
  
    DO i=1,UBOUND(fun,1)
      DO m=1,2
         tn(m,(i-1)*2+m)=fun(i)
      END DO
  END DO
  
  !tn=TRANSPOSE(nt)
  mm=MATMUL(TRANSPOSE(tn),TRANSPOSE(nt))*dens

 !  DO i=1,nod*2
 !   DO m=1,nod*2
 !     mm_L(i,i)=mm_L(i,i)+mm(i,m)
 !   END DO
 !  END DO
 !  mm=mm_L
 ! mm=mm*dens !the thickness and width of the element are considered as one
 END IF

RETURN
END SUBROUTINE beam_mm_mf

      SUBROUTINE beam_km_2(km,Young,poisson,size)
!
! This subroutine forms the stiffness matrix of a
! beam element (bending only).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::Young,poisson,size
 REAL(iwp),INTENT(OUT)::km(:,:)
 REAL(iwp)::two=2.0_iwp,lamda,mu,one=1.0_iwp,zero=0.0_iwp

 km=zero
 lamda=Young*poisson/((one+poisson)*(one-two*poisson))
 mu=Young/(two*(one+poisson))

 km(1,1)=mu
 km(3,1)=-mu
 km(2,2)=lamda+two*mu
 km(4,2)=-(lamda+two*mu)

 km(1,3)=-mu
 km(3,3)=mu
 km(2,4)=-(lamda+two*mu)
 km(4,4)=(lamda+two*mu)

 km=km/size

RETURN
END SUBROUTINE beam_km_2


     SUBROUTINE sample(element,s,wt)
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
   SELECT CASE(nip)
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
 CASE('tetrahedron')
!                       for tetrahedra weights multiplied by 1/6
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
END SUBROUTINE sample

SUBROUTINE vmdpl(dee,stress,pl)
!
! This subroutine forms the plastic stress/strain matrix
! for a von-Mises material.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:),dee(:,:)
 REAL(iwp),INTENT(OUT)::pl(:,:)
 REAL(iwp),ALLOCATABLE::dfds(:),ddfds(:)
 REAL(iwp)::t2,t6,t10,t14,t16,t17,t18,t19,t21,t22,t23,t25,t26,t30,sq2,sx, &
   sy,sz,txy,tyz,tzx,one=1.0_iwp,two=2.0_iwp,d4=4.0_iwp,d6=6.0_iwp
 REAL(iwp)::denom
 INTEGER::i,j,ih
 ih=SIZE(stress)
 ALLOCATE(dfds(ih),ddfds(ih))
 sq2=SQRT(two)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3)
   sz=stress(4)
   t2=sx**2
   t6=sy**2
   t10=sz**2
   t14=txy**2
   t17=SQRT(two*t2-two*sx*sy+two*t6-two*sy*sz+two*t10-two*sz*sx+d6*t14)
   t19=one/sq2/t17
   t21=two*sy
   t22=two*sz
   t26=two*sx
   dfds(1)=t19*(d4*sx-t21-t22)/two
   dfds(2)=t19*(-t26+d4*sy-t22)/two
   dfds(3)=t19*d6*txy
   dfds(4)=t19*(-t21+d4*sz-t26)/two
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)
   txy=stress(4)
   tyz=stress(5)
   tzx=stress(6)
   t2=sx**2
   t6=sy**2
   t10=sz**2
   t14=txy**2
   t16=tyz**2
   t18=tzx**2
   t21=SQRT(two*t2-two*sx*sy+two*t6-two*sy*sz+two*t10-                    &
     two*sz*sx+d6*t14+d6*t16+d6*t18)
   t23=one/sq2/t21
   t25=two*sy
   t26=two*sz
   t30=two*sx
   dfds(1)=t23*(d4*sx-t25-t26)/two
   dfds(2)=t23*(-t30+d4*sy-t26)/two
   dfds(3)=t23*(-t25+d4*sz-t30)/two
   dfds(4)=t23*d6*txy
   dfds(5)=t23*d6*tyz
   dfds(6)=t23*d6*tzx
 END SELECT
 ddfds=MATMUL(dee,dfds)
 denom=DOT_PRODUCT(ddfds,dfds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddfds(i)*ddfds(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,ddfds)
RETURN
END SUBROUTINE vmdpl

SUBROUTINE mcdpl(phi,psi,dee,stress,pl)
!
! This subroutine forms the plastic stress/strain matrix
! for a Mohr-Coulomb material (phi,psi in degrees).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::stress(:),dee(:,:),phi,psi
 REAL(iwp),INTENT(OUT)::pl(:,:)
 REAL(iwp),ALLOCATABLE::dfds(:),dqds(:),ddqds(:),dfdsd(:)
 REAL(iwp)::t1,t2,t3,t4,t5,t6,t8,t10,t12,t13,t14,t15,t16,t17,t18,t19,t20, &
   t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,   &
   t38,t39,t40,t41,t42,t43,t44,t45,t46,t48,t49,t50,t51,t53,t54,t55,t56,   &
   t60,t61,t63,t64,t68,t69,t70,t71,t73,t74,t77,t79,t80,t82,t83,t84,t85,   &
   t86,t89,t92,t93,t94,t97,t98,t101,t103,t106,t110,t111,t113,t122,t129,   &
   t133,t140,t145,t152,t166,t186,t206,pm,pi,phir,snph,snth,sq3,sx,sy,sz,  &
   txy,tyz,tzx,zero=0.0_iwp,pt49=0.49_iwp,one=1.0_iwp,two=2.0_iwp,        &
   d3=3.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d180=180.0_iwp,snps,psir
 REAL(iwp)::denom
 INTEGER::i,j,ih
 ih=SIZE(stress)
 ALLOCATE(dfds(ih),dqds(ih),ddqds(ih),dfdsd(ih))
 pi=ACOS(-one)
 phir=phi*pi/d180
 snph=SIN(phir)
 psir=psi*pi/d180
 snps=SIN(psir)
 sq3=SQRT(d3)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(3)
   sz=stress(4)
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t8=sz-t4
   t10=txy**2
   t16=(sx-sy)**2
   t18=(sy-sz)**2
   t20=(sz-sx)**2
   t25=SQRT((t16+t18+t20)/d6+t10)
   t26=t25**2
   t30=d3*sq3*((sx-t4)*(sy-t4)*t8-t8*t10)/two/t26/t25
   IF(t30>one)t30=one
   IF(t30<-one)t30=-one
   t31=ASIN(t30)
   t33=SIN(t31*t3)
   snth=-t33
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     if(snth.LT.zero)pm=one
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t14=SQRT((t4+t6+t8)*t10+t12)
     t17=one/t14/sq3
     t18=one/two
     t19=t17*t18
     t21=d3+pm*snph
     t23=two*sy
     t24=two*sz
     t31=two*sx
     dfds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dfds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dfds(3)=t17*t18*t21*txy
     dfds(4)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
     t2=snps/d3
     t21=d3+pm*snps
     dqds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dqds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dqds(3)=t17*t18*t21*txy
     dqds(4)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=(t4+t6+t8)*t10+t12
     t14=SQRT(t13)
     t16=d3*sq3
     t18=(sx+sy+sz)*t1
     t19=sx-t18
     t20=sy-t18
     t21=t19*t20
     t22=sz-t18
     t25=t21*t22-t22*t12
     t26=one/two
     t28=t14**2
     t30=one/t28/t14
     t31=t16*t25*t26*t30
     IF(t31>one)t31=one
     IF(t31<-one)t31=-one
     t33=ASIN(t31)
     t34=t33*t1
     t35=COS(t34)
     t36=SIN(t34)
     t38=one/sq3
     t41=one/t14*(t35+t36*snph*t38)
     t43=two*sy
     t44=two*sz
     t46=(d4*sx-t43-t44)*t10
     t49=one-t1
     t53=t19*t1*t22
     t54=t21*t1
     t55=t1*t12
     t60=t16*t25
     t61=t28**2
     t64=t26/t61/t14
     t68=t16*(t49*t20*t22-t53-t54+t55)*t26*t30-d3/two*t60*t64*t46
     t70=d3**2
     t71=sq3**2
     t73=t25**2
     t74=two**2
     t77=t13**2
     t83=SQRT(one-t70*t71*t73/t74/t77/t13)
     t84=one/t83
     t85=t84*t1
     t89=t2*t38
     t94=two*sx
     t97=(-t94+d4*sy-t44)*t10
     t101=t1*t20*t22
     t111=t16*(-t101+t19*t49*t22-t54+t55)*t26*t30-d3/two*t60*t64*t97
     t129=-two*t16*t22*txy*t26*t30-d3*t60*t64*txy
     t140=(-t43+d4*sz-t94)*t10
     t152=t16*(-t101-t53+t21*t49-t49*t12)*t26*t30-d3/two*t60*t64*t140
     dfds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dfds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dfds(3)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dfds(4)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
     t2=snps*t1
     t41=one/t14*(t35+t36*snps*t38)
     t89=t2*t38
     dqds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dqds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dqds(3)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dqds(4)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
   END IF
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)
   txy=stress(4)
   tyz=stress(5)
   tzx=stress(6)
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t5=sx-t4
   t6=sy-t4
   t8=sz-t4
   t10=tyz**2
   t12=tzx**2
   t14=txy**2
   t23=(sx-sy)**2
   t25=(sy-sz)**2
   t27=(sz-sx)**2
   t32=SQRT((t23+t25+t27)/d6+t14+t10+t12)
   t33=t32**2
   t37=d3*sq3*(t5*t6*t8-t5*t10-t6*t12-t8*t14+two*txy*tyz*tzx)/two/t33/t32
   IF(t37>one)t37=one
   IF(t37<-one)t37=-one
   t38=ASIN(t37)
   t40=SIN(t38*t3)
   snth=-t40
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     IF(snth.LT.zero)pm=one
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t16=SQRT((t4+t6+t8)*t10+t12+t13+t14)
     t19=one/t16/sq3
     t20=one/two
     t21=t19*t20
     t23=d3+pm*snph
     t25=two*sy
     t26=two*sz
     t33=two*sx
     t48=t20*t23
     dfds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dfds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dfds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dfds(4)=t19*t48*txy
     dfds(5)=t19*t48*tyz
     dfds(6)=t19*t48*tzx
     t2=snps/d3
     t23=d3+pm*snps
     t48=t20*t23
     dqds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dqds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dqds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dqds(4)=t19*t48*txy
     dqds(5)=t19*t48*tyz
     dqds(6)=t19*t48*tzx
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t15=(t4+t6+t8)*t10+t12+t13+t14
     t16=SQRT(t15)
     t18=d3*sq3
     t20=(sx+sy+sz)*t1
     t21=sx-t20
     t22=sy-t20
     t23=t21*t22
     t24=sz-t20
     t29=two*txy
     t32=t23*t24-t21*t13-t22*t14-t24*t12+t29*tyz*tzx
     t33=one/two
     t35=t16**2
     t37=one/t35/t16
     t39=t18*t32*t33*t37
     IF(t39>one)t39=one
     IF(t39<-one)t39=-one
     t40=ASIN(t39)
     t41=t40*t1
     t42=COS(t41)
     t43=SIN(t41)
     t45=one/sq3
     t48=one/t16*(t42+t43*snph*t45)
     t50=two*sy
     t51=two*sz
     t53=(d4*sx-t50-t51)*t10
     t56=one-t1
     t60=t21*t1*t24
     t61=t23*t1
     t63=t1*t14
     t64=t1*t12
     t69=t18*t32
     t70=t35**2
     t73=t33/t70/t16
     t77=t18*(t56*t22*t24-t60-t61-t56*t13+t63+t64)*t33*t37-               &
       d3/two*t69*t73*t53
     t79=d3**2
     t80=sq3**2
     t82=t32**2
     t83=two**2
     t86=t15**2
     t92=SQRT(one-t79*t80*t82/t83/t86/t15)
     t93=one/t92
     t94=t93*t1
     t98=t2*t45
     t103=two*sx
     t106=(-t103+d4*sy-t51)*t10
     t110=t1*t22*t24
     t113=t1*t13
     t122=t18*(-t110+t21*t56*t24-t61+t113-t56*t14+t64)*t33*t37-           &
       d3/two*t69*t73*t106
     t133=(-t50+d4*sz-t103)*t10
     t145=t18*(-t110-t60+t23*t56+t113+t63-t56*t12)*t33*t37-               &
       d3/two*t69*t73*t133
     t166=t18*(-two*t24*txy+two*tyz*tzx)*t33*t37-d3*t69*t73*txy
     t186=t18*(-two*t21*tyz+t29*tzx)*t33*t37-d3*t69*t73*tyz
     t206=t18*(-two*t22*tzx+t29*tyz)*t33*t37-d3*t69*t73*tzx
     dfds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dfds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dfds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dfds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dfds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dfds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
     t2=snps*t1
     t48=one/t16*(t42+t43*snps*t45)
     t98=t2*t45
     dqds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dqds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dqds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dqds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dqds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dqds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
   END IF
 END SELECT
 ddqds=MATMUL(dee,dqds)
 dfdsd=MATMUL(dee,dfds)
 denom=DOT_PRODUCT(dfdsd,dqds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddqds(i)*dfdsd(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,dqds,ddqds,dfdsd)
RETURN
END SUBROUTINE mcdpl

SUBROUTINE vmflow(stress,dsbar,vmfl)
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
END SUBROUTINE vmflow

   SUBROUTINE fmacat(vmfl,acat)
   !
   ! This subroutine sets up an intermediate matrix acat.
   !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::vmfl(:)
    REAL(iwp),INTENT(OUT)::acat(:,:)
    REAL(iwp)::temp(4,4),zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,three=3.0_iwp
    INTEGER::i,j
    temp=zero
    temp(1,1)=one
    temp(1,2)=-pt5
    temp(1,4)=-pt5
    temp(2,1)=-pt5
    temp(2,2)=one
    temp(2,4)=-pt5
    temp(3,3)=three
    temp(4,1)=-pt5
    temp(4,2)=-pt5
    temp(4,4)=one
    DO i=1,4
      DO j=1,4
        acat(i,j)=vmfl(i)*vmfl(j)
      END DO
    END DO
    acat=temp-acat
   RETURN
   END SUBROUTINE fmacat

      SUBROUTINE plasticmod(k,Maxplastic,SModulus,h)
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::k,Maxplastic,SModulus
    REAL(iwp),INTENT(OUT)::h
    REAL(iwp)::top,bot

    IF(k>=Maxplastic)h=0.0
    IF(k<Maxplastic)h=SModulus

   RETURN
   END SUBROUTINE plasticmod

       SUBROUTINE datafile(Maxplastic,mpcr,cohesion,smethod,E,V,tol,k0,nsrf,dtim, &
                   damping,frictfact)
      !--This subroutine creates a file with the information used in the analysis
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::Maxplastic,mpcr,cohesion,E,V,tol,k0,dtim,damping,frictfact
      INTEGER,INTENT(IN)::smethod,nsrf

      OPEN(1,FILE="Paraview/Information.dat")

      write(1,*) '#'
      write(1,*) '# This file contains the data used for the analysis'
      write(1,*) '#'
      write(1,*) '# Analysis used:1=MPM, 2=CMPM, 3=GIMP',smethod
      write(1,*) '#'
      write(1,*) '# Maximum deviatoric strain:',Maxplastic
      write(1,*) '# Cohesion:',cohesion
      write(1,*) '# Residual cohesion:',mpcr
      write(1,*) '# Youngs modulus:',E
      write(1,*) '# Poisson ration:',V
      write(1,*) '# Tolerance:',tol
      write(1,*) '# k0:',k0
      write(1,*) '# dtim:',dtim
      write(1,*) '# damping:',damping
      write(1,*) '# Friction factor between soil and boundarie:',frictfact
      write(1,*) '# nsrf:',nsrf
      CLOSE(1)

                   END SUBROUTINE datafile

    SUBROUTINE Sloan_KT(lode_theta,Theta_t,phi,Ktheta,DTDT)

     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

     REAL(iwp),INTENT(IN)::lode_theta,Theta_t,phi
     REAL(iwp),INTENT(OUT)::Ktheta,DTDT
     REAL(iwp)::Aval,Bval,three=3.0_iwp,lode_deg,one=1.0_iwp,phi_rad,sign_t,theta_trad

     lode_deg=lode_theta*180.0_iwp/3.14159265_iwp
     phi_rad=phi*3.14159265_iwp/180.0_iwp
     theta_trad=Theta_t*3.14159265_iwp/180.0_iwp


     IF(ABS(lode_deg)>Theta_t)THEN
        IF(lode_deg>=0.0)sign_t=1.0
        IF(lode_deg<0.0)sign_t=-1.0

        Aval=(one/three)*cos(theta_trad)*(three+tan(theta_trad)*tan(three*theta_trad)+ &
                (one/sqrt(three)*sign_t*(tan(three*theta_trad)-three*tan(theta_trad))*sin(phi_rad)))

        Bval=one/(three*cos(three*theta_trad))*(sign_t*sin(theta_trad)+(one/sqrt(three)*sin(phi_rad)*cos(lode_theta)))

        Ktheta=Aval-Bval*sin(three*lode_theta)
        DTDT=-three*Bval*cos(three*lode_theta)

     ELSE
        DTDT=-sin(lode_theta)-(one/sqrt(three))*sin(phi_rad)*cos(lode_theta)
        Ktheta=(cos(lode_theta)-(one/sqrt(three))*sin(phi_rad)*sin(lode_theta))
     END IF


    RETURN
    END SUBROUTINE Sloan_KT

    SUBROUTINE mocouq_Sloan(coh,psi,dsbar,theta,dq1,dq2,dq3,TLoad,a_slo,DTDT)
     !
     ! This subroutine forms the derivatives of a Mohr-Coulomb potential
     ! function with respect to the three invariants (psi in degrees).
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::psi,dsbar,theta,TLoad,a_slo,coh,DTDT
      REAL(iwp),INTENT(OUT)::dq1,dq2,dq3
      REAL(iwp)::psi_rad,three=3.0_iwp,two=2.0_iwp,alpha

      psi_rad=psi*3.14159265_iwp/180.0_iwp

      alpha=dsbar*TLoad/sqrt((dsbar**2)*(TLoad**2)+(a_slo**2)*(sin(psi_rad)**2))
      dq1=sin(psi_rad)
      dq2=(TLoad-tan(three*theta)*DTDT)*alpha
      dq3=(-(sqrt(three)/(two*cos(three*theta)*dsbar**2))*DTDT)*alpha



     RETURN
     END SUBROUTINE mocouq_Sloan

 	SUBROUTINE form_Ai_Mohr(softval,dQ,h,a_slo,dsbar,sigm,coh,phif,psi,KT,theta,TLoad,Ai,dFdc)

	IMPLICIT NONE
	INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::softval
	REAL(iwp),INTENT(IN)::dQ,h,a_slo,coh,phif,psi,KT,dsbar,theta,sigm,TLoad
    REAL(iwp),INTENT(OUT)::Ai,dFdc
    REAL(iwp)::TOP_f,BOT_f,d180=180.0_iwp,pi,psirad,two=2.0_iwp,cotpsi,one=1.0_iwp,   &
    			phirad,sinps,three=3.0_iwp,cospsi,afact,C1,dA,dB,costht,thantht, &
                than3tht,cos3tht,C2,C3,LOADTRad,signth,sin3th,sinth

	pi=ACOS(-one)
    LOADTRad=TLoad*pi/d180
	phirad=phif*pi/d180
    psirad=psi*pi/d180
    sinth=SIN(theta)
    sin3th=SIN(three*theta)
    sinps=SIN(psirad)
    cospsi=COS(psirad)
    cotpsi=one/TAN(psirad)
    afact=a_slo*coh/TAN(psirad)
    costht=COS(LOADTRad)
    thantht=TAN(LOADTRad)
    than3tht=TAN(three*LOADTRad)
    cos3tht=COS(three*LOADTRad)

    C1=0.0
    C2=0.0
    C3=0.0

    IF(theta>0.0)THEN
    signth=one
    ELSE
      signth=-one
    END IF


    IF(softval==1)THEN !Softening reducing cohesino
	TOP_f=sinps**2*two*(a_slo*coh*cotpsi)*(a_slo*cotpsi)
    BOT_f=two*sqrt(dsbar**2*KT**2+(a_slo*coh*cotpsi)**2*sinps**2)
    Ai=-(two/three)*sqrt(dQ)*h*(TOP_f/BOT_f-cospsi)
	dFdc=TOP_f/BOT_f-cospsi
    ELSE   !Softening reducing psi
     IF(ABS(theta)<=LOADTRad)THEN
       C1=2.0*dsbar**2*(-one/sqrt(three)*sinps*sinth)*(-one/sqrt(three)*cospsi*sinth)*(one/sqrt(three))
     ELSE
       dA=(one/three)*costht*((one/SQRT(three))*signth*(than3tht-three*thantht)*cospsi)
       dB=one/(three*cos3tht)*(one/SQRT(three))*cospsi*costht*sin3th
       C1=(dA-dB)*(one/sqrt(three))
     END IF
     	C2=afact**2*two*cospsi*sinps+sinps**2*two*afact*(-a_slo*coh/sinps**2)
        C3=-coh*sinps
        Ai=-(two/three)*sqrt(dQ)*h*(sigm*cospsi+(one/(two*sqrt(dsbar**2*KT**2+afact**2*sinps**2))*(C1+C2))-C3)
        dFdc=(sigm*cospsi+(one/(two*sqrt(dsbar**2*KT**2+afact**2*sinps**2))*(C1+C2))-C3)
    END IF

    END SUBROUTINE form_Ai_Mohr

     SUBROUTINE Small_stiffness(Young,Poiss,Gtan,Gmin,Gamma_r,eps_acum_rev,Unloading,reloading,eps_acum,w,alpha,beta)

        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::Gamma_r,Gtan,Poiss,Gmin,alpha,beta
        REAL(iwp),INTENT(IN)::eps_acum,eps_acum_rev
        LOGICAL,INTENT(INOUT)::Unloading,reloading
        REAL(iwp):: Gval,fac,fac2,zero=0.0_iwp,one=1.0_iwp,two=2.0_iwp
        REAL(iwp),INTENT(INOUT)::Young
        INTEGER,INTENT(IN)::w

       !-Ramber-Osgood
       IF(w==1)THEN
            Young=Gtan*two*(one+Poiss)
       ELSE IF(w>1)THEN
           IF(Unloading)THEN
                   fac=abs((eps_acum-eps_acum_rev)/(two*Gamma_r))
                   fac2=(eps_acum-eps_acum_rev)/(two)
                    !IF(eps_acum_rev<1.0e-10)fac=zero
                    !IF(eps_acum-eps_acum_rev<1.0e-10)fac2=zero
                    IF(abs(fac2)<1.0e-20)THEN
                        Gval=Gtan!(Gtan*fac2)/((one+fac)*fac2)
                    ELSE
                        Gval=(Gtan*fac2)/((one+beta*abs(fac)**(alpha-one))*fac2)
                    END IF

                    IF(Gval<Gmin)Gval=Gmin
                    Young=Gval*two*(one+Poiss)
           ELSE
                    fac=abs(eps_acum/Gamma_r)
                    IF(abs(eps_acum)<1.0e-10)fac=zero
                    Gval=(Gtan)/(one+beta*abs(fac)**(alpha-one))
                    IF(Gval<Gmin)Gval=Gmin
                    Young=Gval*two*(one+Poiss)
           END IF
       END IF

        !-Prevost


    END SUBROUTINE Small_stiffness

    SUBROUTINE Reversal(eps_acum_rev,Unloading,eps_acum,eps_acum_prev,stiffness_rev,reloading,w)

        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::eps_acum_prev,eps_acum
        REAL(iwp),INTENT(INOUT)::eps_acum_rev
        LOGICAL,INTENT(INOUT)::Unloading,stiffness_rev,reloading
        REAL(iwp):: Gval,fac,zero=0.0_iwp,one=1.0_iwp,two=2.0_iwp
        INTEGER,INTENT(IN)::w

        IF(abs(eps_acum_prev)>abs(eps_acum_rev))eps_acum_rev=eps_acum_prev
        IF(abs(eps_acum)>abs(eps_acum_prev).and.abs(eps_acum)>=(abs(eps_acum_rev)+1.0e-5))THEN
            Unloading=.false.
            IF(reloading)stiffness_rev=.true.
            Reloading=.false.
        ELSE IF(abs(eps_acum)<abs(eps_acum_prev).and.(abs(eps_acum)<=abs(eps_acum_rev).or.abs(eps_acum_rev)<1.0e-20))THEN!(abs(eps_acum)<=abs(eps_acum_rev).or.(abs(eps_acum)<abs(eps_acum_prev).and.abs(eps_acum_rev)<1.0e-20)))THEN

            IF(abs(eps_acum)<abs(eps_acum_prev).and.abs(eps_acum_prev)==abs(eps_acum_rev).and. .not.stiffness_rev)THEN
                stiffness_rev=.true.
            END IF
            Unloading=.true.

         ELSE IF(abs(eps_acum)>abs(eps_acum_prev).and.abs(eps_acum)<abs(eps_acum_rev))THEN
            Unloading=.true.
            reloading=.true.
         END IF


    END SUBROUTINE Reversal

   SUBROUTINE from3Dto2DD(deeumat,dee)
 !-This soubrutine creates the 2D stiffness matrix from the 2D stiffnes matrix

  IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(INOUT)::dee(:,:)
 REAL(iwp),INTENT(INOUT)::deeumat(:,:)
 REAL(iwp)::zero=0.0_iwp

 dee=zero
 IF(size(dee,1)==4)THEN
 dee(1:2,1:2)=deeumat(1:2,1:2)
 dee(1,4)=deeumat(1,3)
 dee(2,4)=deeumat(2,3)
 dee(4,1)=deeumat(3,1)
 dee(4,2)=deeumat(3,2)
 dee(4,4)=deeumat(3,3)
 dee(3,3)=deeumat(6,6)
 ELSE
 dee(1:2,1:2)=deeumat(1:2,1:2)
 dee(3,3)=deeumat(6,6)
 END IF


 RETURN
 END SUBROUTINE from3Dto2DD


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

SUBROUTINE formtb2(pb,km,g)
!
! This subroutine assembles an unsymmetrical band matrix pb from
! element constituent matrices km.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::km(:,:)
 INTEGER,INTENT(IN)::g(:)
 REAL(iwp),INTENT(OUT)::pb(:,:)
 INTEGER::i,j,idof,icd,kg,ig,jg,k
 idof=SIZE(km,1)
 k=(SIZE(pb,1)-1)/3
 !n = SIZE(pb,2)
 DO i=1,idof
   ig = g(i)
   IF(ig/=0)THEN
     DO j=1,idof
       jg = g(j)
       !if (ig<=min(m,jg+k)
       IF(jg/=0)THEN
         kg = k+k+1+ig-jg
         pb(kg,jg) = pb(kg,jg) + km(i,j)
       END IF
     END DO
   END IF
 END DO
RETURN
END SUBROUTINE formtb2

SUBROUTINE gauss_band(pb,work)
!
! This subroutine performs gaussian reduction of an unsymmetric
! banded matrix pb. Array work used as working space.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::pb(:,:),work(:,:)
 REAL(iwp)::s,zero=0.0_iwp,small=1.e-10_iwp
 INTEGER::n,iwp1,iq,iqp,iwp11,i,j,k,l,ip,k1
 n=UBOUND(pb,1)
 iwp1=(UBOUND(pb,2)-1)/2+1
 iq=2*iwp1-1
 iqp=iwp1
 iwp11=iwp1-1
 DO i=1,iwp11
   DO j=1,iq
     IF(j>=iwp1+i)THEN
       pb(i,j)=zero
       pb(n-i+1,j)=zero
     ELSE
       pb(i,j)=pb(i,j+iwp1-i)
     END IF
   END DO
 END DO
 DO k=1,n
   l=k+iwp1-1
   IF(l>n)l=n
   ip=0
   s=small
   DO i=k,l
     IF(ABS(pb(i,1))<=s)CYCLE
     s=ABS(pb(i,1))
     ip=i
   END DO
   IF(ip==0)THEN
     WRITE(6,'("singular")')
     EXIT
   END IF
   IF(k==n)EXIT
   work(iwp1,k)=ip
   iqp=iqp-1
   j=iwp1+ip-k
   IF(iqp<j)iqp=j
   IF(j/=iwp1)THEN
     DO j=1,iqp
       s=pb(k,j)
       pb(k,j)=pb(ip,j)
       pb(ip,j)=s
     END DO
   END IF
   k1=k+1
   DO i=k1,l
     s=pb(i,1)/pb(k,1)
     DO j=2,iq
       IF(j>iqp)THEN
         pb(i,j-1)=pb(i,j)
       ELSE
         pb(i,j-1)=pb(i,j)-s*pb(k,j)
       END IF
     END DO
     pb(i,iq)=zero
     work(i-k,k)=s
   END DO
 END DO
RETURN
END SUBROUTINE gauss_band

SUBROUTINE solve_band(pb,work,loads)
!
! This subroutine performs Gaussian forward and back-substitution
! on the reduced unsymmetric band matrix pb.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::pb(:,:),work(:,:)
 REAL(iwp),INTENT(OUT)::loads(0:)
 INTEGER::iwp1,n,n1,i,iv,l,iq,iv1
 REAL(iwp)::s,pt5=0.5_iwp
 iwp1=(UBOUND(pb,2)-1)/2+1
 n=UBOUND(pb,1)
 iq=2*iwp1-1
 n1=n-1
 DO iv=1,n1
   i=INT(work(iwp1,iv)+pt5)
   IF(i/=iv)THEN
     s=loads(iv)
     loads(iv)=loads(i)
     loads(i)=s
   END IF
   l=iv+iwp1-1
   IF(l>n)l=n
   iv1=iv+1
   DO i=iv1,l
     loads(i)=loads(i)-work(i-iv,iv)*loads(iv)
   END DO
 END DO
 loads(n)=loads(n)/pb(n,1)
 iv=n-1
 DO WHILE(iv/=0)
   s=loads(iv)
   l=iq
   IF(iv+l-1>n)l=n-iv+1
   DO i=2,l
     s=s-pb(iv,i)*loads(iv+i-1)
     loads(iv)=s/pb(iv,1)
   END DO
 iv=iv-1
 END DO
RETURN
END SUBROUTINE solve_band

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


SUBROUTINE cross_product(b,c,a)
!
! This subroutine forms the cross product of two vectors, a = b x c
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::b(:),c(:)
 REAL(iwp),INTENT(OUT)::a(:,:)
 INTEGER::ib,ic,i,j
 ib=SIZE(b)
 ic=SIZE(c)
 DO i=1,ib
   DO j=1,ic
     a(i,j)=b(i)*c(j)
   END DO
 END DO
RETURN
END SUBROUTINE cross_product

SUBROUTINE sparin_gauss(kv,kdiag)
!
! This subroutine performs Gaussian factorisation of a skyline matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::kdiag(:)
 REAL(iwp),INTENT(OUT)::kv(:)
 REAL(iwp)::num,den,fac,zero=0.0_iwp
 INTEGER::n,ii,i,j,k,l,kk,l1,l2,l3
 n=UBOUND(kdiag,1)
 DO j=1,n-1
   den=kv(kdiag(j))
   ii=0
   DO i=j+1,n
     ii=ii+1
     l=kdiag(i)-ii
     IF(l-kdiag(i-1)>zero)THEN
       num=kv(l)
       fac=num/den
       kk=-1
       DO k=i,n
         kk=kk+1
         l1=kdiag(i+kk)-kk
         l2=l1-ii
         l3=kdiag(i+kk-1)
         IF(l2-l3>zero)kv(l1)=kv(l1)-fac*kv(l2)
       END DO
     END IF
   END DO
 END DO
RETURN
END SUBROUTINE sparin_gauss

SUBROUTINE spabac_gauss(kv,loads,kdiag)
!
! This subroutine performs Gaussian forwrad and back-substitution on a
! skyline matrix.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::kv(:)
 REAL(iwp),INTENT(IN OUT)::loads(0:)
 INTEGER,INTENT(IN)::kdiag(:)
 REAL(iwp)::num,den,fac,asum,zero=0.0_iwp
 INTEGER::i,j,l,n,ii,jj,l1,l2
 n=UBOUND(kdiag,1)
 DO j=1,n-1
   den=kv(kdiag(j))
   ii=0
   DO i=j+1,n
     ii=ii+1
     l=kdiag(i)-ii
     IF(l-kdiag(i-1)>zero)THEN
       num=kv(l)
       fac=num/den
       loads(i)=loads(i)-fac*loads(j)
     END IF
   END DO
 END DO
 loads(n)=loads(n)/kv(kdiag(n))
 DO i=n-1,1,-1
   jj=0
   asum=zero
   DO j=i+1,n
     jj=jj+1
     l1=kdiag(i+jj)-jj
     l2=kdiag(i+jj-1)
     IF(l1-l2>zero)asum=asum+kv(l1)*loads(j)
   END DO
   loads(i)=(loads(i)-asum)/kv(kdiag(i))
 END DO
RETURN
END SUBROUTINE spabac_gauss






SUBROUTINE TransM(Tmatrix,alpha)

IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::alpha
 REAL(iwp),INTENT(IN OUT)::Tmatrix(:,:)

 INTEGER::i,j,m

 j=SIZE(Tmatrix,1)
 m=1
 DO i=1,j/2
     Tmatrix(m,m)=COS(alpha)
     Tmatrix(m+1,m+1)=COS(alpha)
     Tmatrix(m+1,m)=SIN(alpha)
     Tmatrix(m,m+1)=-SIN(alpha)
     m=m+2
 END DO

END SUBROUTINE TransM










end module main
