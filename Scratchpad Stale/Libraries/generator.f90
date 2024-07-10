MODULE GEOMETRY_GENERATOR
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RECTANGULAR_2D(coord,num,nn,nels,nx,ny,w1,h1,node,offsetx,offsety)
    !
    ! This subroutine forms the nodal coordinates and numbering for a 2-d
    ! of n-noded rectangles. Nodes numbering in the y-direction,
    ! elements numbered in the x-direction.
    !
    ! GLOBAL SYSTEM
    !                  facx
    !              |----------|
    !
    !          🎯(0,h1)                🎯(w1,h1)    
    !              1----------4----------7       −.
    !              |          |          |        |
    !              |  3(1,2)  |  4(2,2)  |        | 
    !              |          |          |        |
    !        .-    2----------5----------8        |  'ny' elements
    !        |     |          |          |        |
    !   facy |     |  1(1,1)  |  2(2,1)  |        | 
    !        |     |  (ip,iq) |          |        |
    !        !-    3----------6----------9       -!
    !          🎯(0,0)                 🎯(w1,0)    
    !        
    !              |---------------------|
    !                   'nx' elements
    !
    ! LOCAL SYSTEM
    !
    !    1----------5----------4   −.
    !    |          .          |    |
    !    |          .          |    | 
    !    |          .          |    |
    !    6- - - - - 9 - - - - -8    |  facy
    !    |          .          |    |
    !    |          .          |    | 
    !    |          .          |    |
    !    2----------7----------3   -!
    !
    !    |---------------------|
    !             facx
    !
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    ! Dummy variable
    INTEGER,INTENT(IN)::nx,ny,node
    INTEGER,INTENT(IN),OPTIONAL::offsetx,offsety
    REAL(iwp),INTENT(IN)::w1,h1
    REAL(iwp),ALLOCATABLE,INTENT(OUT)::coord(:,:)
    INTEGER,ALLOCATABLE,INTENT(OUT)::num(:,:)
    INTEGER,INTENT(OUT)::nn,nels
    ! Local variable
    REAL(iwp)::facx,facy,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp
    INTEGER::i,j,iel,ip,iq,start,end,ndim=2
    ! Calculate number of element and number of nodes
    nels=nx*ny
    SELECT CASE(node)
    CASE(4)
      nn=(nx+1)*(ny+1)
    CASE(8)
      nn=(nx*2+1)*(ny*2+1)-nels
    CASE(9)
      nn=(nx*2+1)*(ny*2+1)
    END SELECT  
    ! Calculate mesh size facx and facy
    facx=w1/nx/2
    facy=h1/ny/2
    ! Determine element node number (ip,iq)
    ALLOCATE(num(node,nels))
    DO iel=1,nels
      iq=(iel-1)/nx+1
      ip=iel-(iq-1)*nx
      SELECT CASE(node)
      CASE(4)
        num(1,iel)=ip*(ny+1)-iq
        num(4,iel)=(ip+1)*(ny+1)-iq
        num(2,iel)=num(1,iel)+1
        num(3,iel)=num(4,iel)+1
      CASE(8)
        num(1,iel)=ip*(ny*2+1)+(ip-1)*(ny+1) - iq*2
        num(7,iel)=ip*(ny*2+1)+ip*(ny+1) - iq
        num(8,iel)=(ip+1)*(ny*2+1)+ip*(ny+1) - iq*2
        num(2,iel)=num(1,iel)+1
        num(3,iel)=num(1,iel)+2
        num(4,iel)=num(8,iel)+1
        num(5,iel)=num(7,iel)+2
        num(6,iel)=num(7,iel)+1
      CASE(9)
        num(1,iel)=(ip*2-1)*(ny*2+1) - iq*2
        num(7,iel)=(ip*2)*(ny*2+1) - iq*2
        num(8,iel)=(ip*2+1)*(ny*2+1) - iq*2
        num(2,iel)=num(1,iel)+1
        num(3,iel)=num(1,iel)+2
        num(4,iel)=num(8,iel)+2
        num(5,iel)=num(7,iel)+2
        num(6,iel)=num(7,iel)+1
        num(9,iel)=num(8,iel)+1
      END SELECT
    END DO
    
    ! Calculate node coordinates
    ALLOCATE(coord(ndim,nn))
    SELECT CASE(node)
    CASE(4)
      DO i=1,nx+1
        start=i*(ny+1)-ny
        end=start+ny
        coord(1,start:end)=(i-1)*facx*2
      END DO
      DO i=1,ny+1
        coord(2,i::ny+1)=(ny+1-i)*facy*2
      END DO
    CASE(8)
    CASE(9)
      DO i=1,nx*2+1
        start=i*(ny*2+1)-ny*2
        end=start+ny
        coord(1,start:end)=i*facy
      END DO
      DO i=1,ny*2+1
        coord(2,i::ny*2+1)=(ny*2*facy)-(i-1)*facy
      END DO
    END SELECT
    IF(present(offsetx)) coord(1,:)=coord(1,:)+facx*2*offsetx
    IF(present(offsety)) coord(2,:)=coord(2,:)+facy*2*offsety
  END SUBROUTINE RECTANGULAR_2D

END MODULE GEOMETRY_GENERATOR