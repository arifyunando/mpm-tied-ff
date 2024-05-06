module geom

SAVE

CONTAINS

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


SUBROUTINE emb_2d_geom2(iel,bod,nx1,ny1,s1,newel,slopeel,row,column,slope1,w1,h1,coord,num,A,B,dist_x,dist_y,slopeopt)
!
! This subroutine forms the nodal coordinates and numbering for a 2-d
! slope of 4-node quadrilaterals. Nodes numbering in the y-direction,
! elements numbered in the x-direction. Considering only slope without
! foundation soil
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::w1,h1,s1
 INTEGER,INTENT(IN)::iel,nx1,ny1,newel,slopeel,slopeopt,bod,dist_x,dist_y
 INTEGER,INTENT(INOUT)::row,column,slope1,A,B
 REAL(iwp),INTENT(OUT)::coord(:,:) !!coord(4,2)
 INTEGER,INTENT(OUT)::num(:)
 REAL(iwp)::dimy,dimx,facs,frh,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,   &
            slopelong,dimx2
 INTEGER::nxe,nye,nc,nt,ip,iq,ntotal
 nxe=nx1  !numero de elementos en X
 nye=ny1   !numero de elementos en Y
 nt=nx1*ny1   !numero de elementos en la parte superior sin la pendiente
 ntotal=nt+slopeel
 nc=(nye+1)*nx1+ny1   !ECUACION MODIFICADA PARA 4 NODOS CORRECTA
 dimy=-h1/ny1          !element dimension in y direction
 dimx=w1/nx1     !element dimension in x direction
 IF(newel>0.0)dimx2=s1/newel
 IF(newel==0.0)dimx2=0.0
 facs=zero
 frh=zero

 IF(column==1 .and. row==1 .and. iel==1)THEN
   num(1)=2
   num(2)=1
   num(3)=num(1)+ny1
   num(4)=num(3)+1
   A=num(3)
   B=num(4)

   !This are the GIMP coordinates, including the boundary elements   
   !x coordenates
    coord(1,1)=dimx*dist_x
    coord(2,1)=dimx*dist_x
    coord(3,1)=dimx*(dist_x+1)   
    coord(4,1)=dimx*(dist_x+1)
   !y coordenates
    coord(1,2)=(row+dist_y)*dimy
    coord(2,2)=dist_y*dimy
    coord(3,2)=dist_y*dimy
    coord(4,2)=(row+dist_y)*dimy

 END IF

 IF(column==1.and.row>1)THEN
   num(1)=row+1
   num(2)=num(1)-1
   num(3)=num(1)+ny1
   num(4)=num(3)+1
   A=num(3)
   B=num(4)
   !This are the GIMP coordinates, including the boundary elements   
   !x coordenates
    coord(1,1)=dimx*(column+dist_x-1)
    coord(2,1)=dimx*(column+dist_x-1)
    coord(3,1)=dimx*(column+dist_x)  
    coord(4,1)=dimx*(column+dist_x)
   !y coordenates
    coord(1,2)=(row+dist_y)*dimy
    coord(2,2)=((row-1)+dist_y)*dimy
    coord(3,2)=((row-1)+dist_y)*dimy
    coord(4,2)=(row+dist_y)*dimy
  END IF  

 IF(column>1.and.row>1.and.column<=nx1+1)THEN
   num(1)=B
   num(2)=A
   num(3)=num(1)+ny1
   num(4)=num(3)+1
   A=num(3)
   B=num(4)
   !This are the GIMP coordinates, including the boundary elements   
   !x coordenates
    coord(1,1)=dimx*(column+dist_x-1)
    coord(2,1)=dimx*(column+dist_x-1)
    coord(3,1)=dimx*(column+dist_x)  
    coord(4,1)=dimx*(column+dist_x)
   !y coordenates
    coord(1,2)=(row+dist_y)*dimy
    coord(2,2)=((row-1)+dist_y)*dimy
    coord(3,2)=((row-1)+dist_y)*dimy
    coord(4,2)=(row+dist_y)*dimy
    
 END IF

 IF(column<=nx1+1 .and. column>1 .and. row==1)THEN
  num(1)=B
   num(2)=A
   num(3)=num(1)+ny1
   num(4)=num(3)+1
   A=num(3)
   B=num(4)
   !x coordenates
    coord(1,1)=dimx*(column+dist_x-1)
    coord(2,1)=dimx*(column+dist_x-1)
    coord(3,1)=dimx*(column+dist_x)  
    coord(4,1)=dimx*(column+dist_x)
    
   IF(column==slope1.and.slopeopt==1)coord(3,1)=coord(2,1)+dimx2/2.0_iwp
   IF(column==slope1.and.slopeopt==1)coord(4,1)=coord(1,1)+dimx2

   !y coordenates
    coord(1,2)=(row+dist_y)*dimy
    coord(2,2)=dist_y*dimy
    coord(3,2)=dist_y*dimy
    coord(4,2)=(row+dist_y)*dimy

   IF(column==slope1.and.slopeopt==1)coord(3,2)=coord(2,2)+(coord(1,2)-coord(2,2))/2.0_iwp
 END IF

  IF(column>nx1+1 .and. iel>1)THEN
   num(1)=B
   num(2)=A
   num(3)=num(1)+ny1-(column-nx1-1)
   num(4)=num(3)+1
   A=num(3)
   B=num(4)
   
   !x coordenates
    coord(1,1)=dimx*(column+dist_x-1)
    coord(2,1)=dimx*(column+dist_x-1)
    coord(3,1)=dimx*(column+dist_x)  
    coord(4,1)=dimx*(column+dist_x)
    
   IF(column==slope1.and.slopeopt==1)coord(3,1)=coord(2,1)+dimx2/2.0_iwp
   IF(column==slope1.and.slopeopt==1)coord(4,1)=coord(1,1)+dimx2

   !y coordenates
    coord(1,2)=(row+dist_y)*dimy
   coord(2,2)=(row*(dist_y-1))*dimy
    coord(3,2)=(row*(dist_y-1))*dimy
    coord(4,2)=(row+dist_y)*dimy
    
   IF(column==slope1.and.slopeopt==1)coord(3,2)=coord(2,2)+(coord(1,2)-coord(2,2))/2.0_iwp
 END IF

   IF(column==slope1)THEN
     row=row+1
     slope1=slope1+1
     IF(slopeopt==2)slope1=slope1-1
     column=0
   END IF
   column=column+1

RETURN
END SUBROUTINE emb_2d_geom2

SUBROUTINE emb_2d_bc(nx1,nx2,ny1,ny2,nf)
!
! This subroutine generates the nf array for a 2-d slope geometry.
!
 IMPLICIT NONE
 INTEGER,INTENT(IN)::nx1,nx2,ny1,ny2
 INTEGER,INTENT(OUT)::nf(:,:)
 INTEGER::nm,ic,i,j,nye
 nye=ny1+ny2
 nm=0
 ic=0
 !DO i=1,2*nye
 DO i=1,nye
   nm=nm+1
   nf(1,nm)=0
   ic=ic+1
   nf(2,nm)=ic
 END DO
 nm=nm+1
 nf(1,nm)=0
 nf(2,nm)=0
 DO j=1,nx1-1
   !DO i=1,nye
    ! nm=nm+1
     !ic=ic+1
     !nf(1,nm)=ic
    ! ic=ic+1
    ! nf(2,nm)=ic
   !END DO
   !nm=nm+1
   !nf(1,nm)=0
   !nf(2,nm)=0
   !DO i=1,2*nye
   DO i=1,nye
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
 END DO
   DO i=1,nye
     nm=nm+1
     ic=ic+1
     nf(1,nm)=ic
     ic=ic+1
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
   !DO i=1,2*ny1
   !DO i=1,ny1
   !  nm=nm+1
   !  ic=ic+1
   !  nf(1,nm)=ic
   !  ic=ic+1
   !  nf(2,nm)=ic
   !END DO
   !IF(nx2==0)THEN
   !DO i=1,2*ny2
   ! DO i=1,ny2
   !  nm=nm+1
   !  nf(1,nm)=0
   !  ic=ic+1
   !  nf(2,nm)=ic
   !END DO
   !nm=nm+1
   !nf(1,nm)=0
   !nf(2,nm)=0
! ELSE
   !DO i=1,2*ny2

   DO j=1,nx2
      DO i=1,ny2
     nm=nm+1
     ic=ic+1
     IF(j==nx2)then
       ic=0
     END IF
     nf(1,nm)=ic
     ic=ic+1
    IF(j==nx2)then
       ic=0
     END IF
     nf(2,nm)=ic
   END DO
   nm=nm+1
   nf(1,nm)=0
   nf(2,nm)=0
  END DO
   !DO j=1,nx2
    ! DO i=1,ny2
     !  nm=nm+1
      ! ic=ic+1
      ! nf(1,nm)=ic
      ! ic=ic+1
      ! nf(2,nm)=ic
     !END DO
     !nm=nm+1
     !nf(1,nm)=0
     !nf(2,nm)=0
     !DO i=1,2*ny2
     !DO i=1,ny2
       !nm=nm+1
      ! IF(j==nx2)THEN
       !  nf(1,nm)=0
       !ELSE
        ! ic=ic+1
         !nf(1,nm)=ic
       !END IF
       !ic=ic+1
       !nf(2,nm)=ic
     !END DO
     !nm=nm+1
     !nf(1,nm)=0
     !nf(2,nm)=0
   !END DO
 !END IF
RETURN
END SUBROUTINE emb_2d_bc




SUBROUTINE emb_2d_geom4(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2,coord,num)
!
! This subroutine forms the nodal coordinates and numbering for a 2-d
! slope of 8-node quadrilaterals. Nodes numbering in the y-direction,
! elements numbered in the x-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::w1,s1,w2,h1,h2
 INTEGER,INTENT(IN)::iel,nx1,nx2,ny1,ny2
 REAL(iwp),INTENT(OUT)::coord(:,:) !!coord(4,2)
 INTEGER,INTENT(OUT)::num(:)
 REAL(iwp)::facx,facy,facb,facs,frh,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp
 INTEGER::nxe,nye,nc,nt,ip,iq
 nxe=nx1+nx2  !numero de elementos en X
 nye=ny1+ny2   !numero de elementos en Y
 nt=nx1*ny1   !numero de elementos sin considerar la cimentacion
 nc=(nye+1)*nx1+ny1   !ECUACION MODIFICADA PARA 4 NODOS CORRECTA
 facx=s1/ny1  ! this is zero
 facy=h1/ny1  !element high
 facb=(w1+s1)/nx1  !element widht
 facs=zero
 IF(ny2/=0)facs=h2/ny2   !This is zero
 frh=zero
 IF(nx2/=0)frh=w2/nx2    !This is zero
 IF(iel<=nt)THEN
   iq=(iel-1)/nx1+1
   ip=iel-(iq-1)*nx1
 ELSE
   iq=(iel-nt-1)/nxe+ny1+1
   ip=iel-nt-(iq-ny1-1)*nxe
 END IF
 IF(ip<=nx1)THEN
   num(1)=(ip-1)*(nye+1)+1*iq+1
   num(2)=num(1)-1
   num(3)=ip*(nye+1)+1*iq
   num(4)=num(3)+1
   IF(iq<=ny1)THEN
     coord(1,1)=(ip-one)*(w1+iq*facx)/nx1
     coord(2,1)=(ip-one)*(w1+(iq-1)*facx)/nx1
     coord(3,1)=ip*(w1+(iq-1)*facx)/nx1
     coord(4,1)=ip*(w1+iq*facx)/nx1
     coord(1,2)=-iq*facy
     coord(2,2)=-(iq-1)*facy
     coord(3,2)=-(iq-1)*facy   !This are original
     coord(4,2)=-iq*facy
!       coord(1,2)=-(iq-1)*facy
!     coord(2,2)=-(iq-2)*facy  !modiyed for the GIMP
!     coord(3,2)=-(iq-2)*facy
!     coord(4,2)=-(iq-1)*facy
   ELSE
     coord(1,1)=(ip-one)*facb
     coord(2,1)=(ip-one)*facb
     coord(3,1)=ip*facb
     coord(4,1)=ip*facb
     coord(1,2)=-h1-(iq-ny1)*facs
     coord(2,2)=-h1-(iq-ny1-1)*facs
     coord(3,2)=-h1-(iq-ny1-1)*facs
     coord(4,2)=-h1-(iq-ny1)*facs
   END IF
 ELSE
   num(1)=nc+(ip-nx1-1)*(ny2+1)+1*(iq-ny1)+1
   num(2)=num(1)-1
   num(3)=nc+(ip-nx1)*(ny2+1)+1*(iq-ny1)
   num(4)=num(3)+1

   coord(1,1)=w1+s1+(ip-nx1-1)*frh
   coord(2,1)=coord(1,1)
   coord(3,1)=w1+s1+(ip-nx1)*frh
   coord(4,1)=coord(3,1)
   coord(1,2)=-h1-(iq-ny1)*facs
   coord(2,2)=-h1-(iq-ny1-1)*facs
   coord(3,2)=-h1-(iq-ny1-1)*facs
   coord(4,2)=-h1-(iq-ny1)*facs
 END IF

RETURN
END SUBROUTINE emb_2d_geom4

SUBROUTINE emb_2d_geom(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2,coord,num)
!
! This subroutine forms the nodal coordinates and numbering for a 2-d
! slope of 8-node quadrilaterals. Nodes numbering in the y-direction,
! elements numbered in the x-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::w1,s1,w2,h1,h2
 INTEGER,INTENT(IN)::iel,nx1,nx2,ny1,ny2
 REAL(iwp),INTENT(OUT)::coord(:,:)
 INTEGER,INTENT(OUT)::num(:)
 REAL(iwp)::facx,facy,facb,facs,frh,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp
 INTEGER::nxe,nye,nc,nt,ip,iq
 nxe=nx1+nx2
 nye=ny1+ny2
 nt=nx1*ny1
 nc=(3*nye+2)*nx1+2*ny1
 facx=s1/ny1
 facy=h1/ny1
 facb=(w1+s1)/nx1
 facs=zero
 IF(ny2/=0)facs=h2/ny2
 frh=zero
 IF(nx2/=0)frh=w2/nx2
 IF(iel<=nt)THEN
   iq=(iel-1)/nx1+1
   ip=iel-(iq-1)*nx1
 ELSE
   iq=(iel-nt-1)/nxe+ny1+1
   ip=iel-nt-(iq-ny1-1)*nxe
 END IF
 IF(ip<=nx1)THEN
   num(1)=(ip-1)*(3*nye+2)+2*iq+1
   num(2)=num(1)-1
   num(3)=num(1)-2
   num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
   num(5)=ip*(3*nye+2)+2*iq-1
   num(6)=num(5)+1
   num(7)=num(5)+2
   num(8)=num(4)+1
   IF(iq<=ny1)THEN
     coord(1,1)=(ip-one)*(w1+iq*facx)/nx1
     coord(3,1)=(ip-one)*(w1+(iq-1)*facx)/nx1
     coord(5,1)=ip*(w1+(iq-1)*facx)/nx1
     coord(7,1)=ip*(w1+iq*facx)/nx1
     coord(1,2)=-iq*facy
     coord(3,2)=-(iq-1)*facy
     coord(5,2)=-(iq-1)*facy
     coord(7,2)=-iq*facy
   ELSE
     coord(1,1)=(ip-one)*facb
     coord(3,1)=(ip-one)*facb
     coord(5,1)=ip*facb
     coord(7,1)=ip*facb
     coord(1,2)=-h1-(iq-ny1)*facs
     coord(3,2)=-h1-(iq-ny1-1)*facs
     coord(5,2)=-h1-(iq-ny1-1)*facs
     coord(7,2)=-h1-(iq-ny1)*facs
   END IF
 ELSE
   num(1)=nc+(ip-nx1-1)*(3*ny2+2)+2*(iq-ny1)+1
   num(2)=num(1)-1
   num(3)=num(1)-2
   num(4)=nc+(ip-nx1-1)*(3*ny2+2)+2*ny2+iq-ny1+1
   num(5)=nc+(ip-nx1)*(3*ny2+2)+2*(iq-ny1)-1
   num(6)=num(5)+1
   num(7)=num(5)+2
   num(8)=num(4)+1
   coord(1,1)=w1+s1+(ip-nx1-1)*frh
   coord(3,1)=coord(1,1)
   coord(5,1)=w1+s1+(ip-nx1)*frh
   coord(7,1)=coord(5,1)
   coord(1,2)=-h1-(iq-ny1)*facs
   coord(3,2)=-h1-(iq-ny1-1)*facs
   coord(5,2)=-h1-(iq-ny1-1)*facs
   coord(7,2)=-h1-(iq-ny1)*facs
 END IF
 coord(2:6:2,:)=pt5*(coord(1:5:2,:)+coord(3:7:2,:))
 coord(8,:)=pt5*(coord(7,:)+coord(1,:))
RETURN
END SUBROUTINE emb_2d_geom

SUBROUTINE emb_2d_geom_9n(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2,coord,num,rowy,colx)
!
! This subroutine forms the nodal coordinates and numbering for a 2-d
! slope of 8-node quadrilaterals. Nodes numbering in the y-direction,
! elements numbered in the x-direction.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::w1,s1,w2,h1,h2
 INTEGER,INTENT(IN)::iel,nx1,nx2,ny1,ny2,rowy,colx
 REAL(iwp),INTENT(INOUT)::coord(:,:)
 INTEGER,INTENT(INOUT)::num(:)
 REAL(iwp)::facx,facy,facb,facs,frh,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,slope,midllex,midlley
 INTEGER::nxe,nye,nc,nt,ip,iq
 nxe=nx1+nx2
 nye=ny1+ny2
 nt=nx1*nye !Elements in the slope and the base of the slope (not the front foundation soil)

 facx=w1/nx1 !Cellsize x direction in the slope section (top of the slope)
 facy=h1/ny1 !Cellsize y direction in the slope section
 facb=(w1+s1)/nx1 !Cellsize x direction (base of the slope)
 facs=zero
 IF(ny2/=0)facs=h2/ny2
 frh=zero
 IF(nx2/=0)frh=w2/nx2
 
 IF(facx>zero)THEN
    midllex=facx/2.0_iwp
 ELSE
    midllex=zero  
 END IF
     
 IF(facy>zero)THEN
   midlley=facy/2.0_iwp
 ELSE
    midlley=zero  
 END IF


     ip=(iel-1)/nye+1
     iq=iel-(ip-1)*nye
     IF(colx==1.and.rowy==1)THEN
     num=0
     num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
     num(5)=ip*2*(2*nye+1)+2*iq-1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     coord=0.0_iwp
     coord(1:3,1)=0.0_iwp

     
     coord(4,1)=coord(3,1)+midllex
     coord(5,1)=coord(3,1)+facx
     slope=(colx*facb-colx*facx)/h1
     coord(6,1)=coord(5,1)+slope*midlley
     coord(7,1)=coord(6,1)+slope*midlley
     coord(8,1)=(coord(7,1)-coord(1,1))/2.0_iwp+coord(1,1)
     coord(9,1)=(coord(6,1)-coord(2,1))/2.0_iwp+coord(2,1)
     
     coord(3:5,2)=0.0_iwp+((rowy-1)*facy)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(7,2)
     coord(8,2)=coord(7,2)
     coord(2,2)=coord(6,2)
     coord(9,2)=coord(6,2)
     ELSE IF(colx==1.and.rowy>1.and.rowy<=ny1)THEN
         
    
     num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
     num(5)=ip*2*(2*nye+1)+2*iq-1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     coord(1:3,1)=0.0_iwp
     coord(4,1)=coord(8,1)
     coord(5,1)=coord(7,1)
     slope=(colx*facb-colx*facx)/h1
     coord(6,1)=coord(5,1)+slope*midlley
     coord(7,1)=coord(6,1)+slope*midlley
     coord(8,1)=(coord(7,1)-coord(1,1))/2.0_iwp+coord(1,1)
     coord(9,1)=(coord(6,1)-coord(2,1))/2.0_iwp+coord(2,1)
     
     coord(3:5,2)=0.0_iwp-((rowy-1)*facy)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
     ELSE IF(colx==1.and.rowy>ny1)THEN
         
     num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
     num(5)=ip*2*(2*nye+1)+2*iq-1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     coord(1:3,1)=0.0_iwp
     coord(4,1)=coord(8,1)
     coord(5,1)=coord(7,1)
     slope=(colx*facb-colx*facx)/h1
     coord(6,1)=coord(5,1)
     coord(7,1)=coord(6,1)
     coord(8,1)=coord(1,1)+facb/2.0_iwp
     coord(9,1)=coord(4,1)
     facy=h2/ny2
     
       IF(facy>zero)THEN
        midlley=facy/2.0_iwp
      ELSE
        midlley=zero  
      END IF
      
     !coord(3:5,2)=0.0_iwp-((rowy-1)*h1/ny1)
     coord(3,2)=coord(1,2)
     coord(4,2)=coord(8,2)
     coord(5,2)=coord(7,2)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
     ELSE IF(colx>1.and.rowy==1)THEN
     
     num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
     num(5)=ip*2*(2*nye+1)+2*iq-1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     coord(3,1)=(colx-1)*facx
     coord(4,1)=coord(3,1)+midllex
     coord(5,1)=coord(3,1)+facx
     slope=(colx*facb-colx*facx)/h1
     coord(6,1)=coord(5,1)+slope*midlley
     coord(7,1)=coord(6,1)+slope*midlley
     slope=((colx-1)*facb-(colx-1)*facx)/h1
     coord(2,1)=coord(3,1)+slope*midlley
     coord(1,1)=coord(3,1)+slope*facy
     coord(8,1)=(coord(7,1)-coord(1,1))/2.0_iwp+coord(1,1)
     coord(9,1)=(coord(6,1)-coord(2,1))/2.0_iwp+coord(2,1)
     
     coord(3:5,2)=0.0_iwp-((rowy-1)*facy)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
    ELSE IF(colx>1.and.rowy>1.and.rowy<=ny1)THEN
     
     num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
     num(5)=ip*2*(2*nye+1)+2*iq-1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     coord(3,1)=coord(1,1)
     coord(4,1)=coord(8,1)
     coord(5,1)=coord(7,1)
     slope=(colx*facb-colx*facx)/h1
     coord(6,1)=coord(5,1)+slope*midlley
     coord(7,1)=coord(6,1)+slope*midlley
     slope=((colx-1)*facb-(colx-1)*facx)/h1
     coord(2,1)=coord(3,1)+slope*midlley
     coord(1,1)=coord(3,1)+slope*facy
     coord(8,1)=(coord(7,1)-coord(1,1))/2.0_iwp+coord(1,1)
     coord(9,1)=(coord(6,1)-coord(2,1))/2.0_iwp+coord(2,1)
     
     coord(3:5,2)=0.0_iwp-((rowy-1)*facy)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
  
    ELSE IF(colx>1.and.colx<nx1+1.and.rowy>ny1)THEN
        
     num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
     num(5)=ip*2*(2*nye+1)+2*iq-1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     coord(3,1)=coord(1,1)
     coord(4,1)=coord(8,1)
     coord(5,1)=coord(7,1)
     slope=(colx*facb-colx*facx)/h1
     coord(6,1)=coord(5,1)
     coord(7,1)=coord(6,1)
     coord(8,1)=coord(1,1)+facb/2.0_iwp
     coord(9,1)=coord(4,1)
     coord(1,1)=(colx-1)*facb
     coord(2,1)=(colx-1)*facb
     facy=h2/ny2
     
       IF(facy>zero)THEN
        midlley=facy/2.0_iwp
      ELSE
        midlley=zero  
      END IF
     !coord(3:5,2)=0.0_iwp-((rowy-1)*h1/ny1)
      coord(3,2)=coord(1,2)
     coord(4,2)=coord(8,2)
     coord(5,2)=coord(7,2)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
   ELSE IF(colx==nx1+1.and.rowy==ny1+1)THEN  
     
     
     num(1)=(nye*2+1)*(nx1*2+1)-((ny2-1)*2)
     num(2)=num(1)-1
     num(3)=num(1)-2
     num(4)=num(3)+(ny2*2)+1
     num(5)=num(4)+(ny2*2)+1
     num(6)=num(5)+1
     num(7)=num(5)+2
     num(8)=num(4)+2
     num(9)=num(4)+1
     
     facx=w2/nx2
     facy=h2/ny2
     
      IF(facx>zero)THEN
        midllex=facx/2.0_iwp
      ELSE
        midllex=zero  
      END IF
     
      IF(facy>zero)THEN
        midlley=facy/2.0_iwp
      ELSE
        midlley=zero  
      END IF
     
     coord(1:3,1)=(colx-1)*facb
     coord(4,1)=coord(3,1)+midlley
     coord(5,1)=coord(3,1)+facx
     coord(6,1)=coord(5,1)
     coord(7,1)=coord(6,1)
     coord(8,1)=coord(1,1)+midlley
     coord(9,1)=coord(4,1)
     
     coord(3:5,2)=0.0_iwp-((rowy-1)*h1/ny1)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
   ELSE IF(colx==nx1+1.and.rowy>=ny1+1)THEN  
       
       
   
     num(3)=num(1)
     num(4)=num(8)
     num(5)=num(7)
     num(6)=num(5)+1
     num(7)=num(6)+1
     num(8)=num(4)+2
     num(9)=num(4)+1
     num(1)=num(3)+2
     num(2)=num(3)+1
     
      facx=w2/nx2
     facy=h2/ny2
     
      IF(facx>zero)THEN
        midllex=facx/2.0_iwp
      ELSE
        midllex=zero  
      END IF
     
      IF(facy>zero)THEN
        midlley=facy/2.0_iwp
      ELSE
        midlley=zero  
      END IF
    
     coord(1:3,1)=(colx-1)*facb
     coord(4,1)=coord(3,1)+midllex
     coord(5,1)=coord(3,1)+facx
     coord(6,1)=coord(5,1)
     coord(7,1)=coord(6,1)
     coord(8,1)=coord(1,1)+midllex
     coord(9,1)=coord(4,1)
     
     !coord(3:5,2)=0.0_iwp-((rowy-1)*h1/ny1)
      coord(3,2)=coord(1,2)
     coord(4,2)=coord(8,2)
     coord(5,2)=coord(7,2)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
   ELSE IF(colx>=nx1+1.and.rowy==ny1+1)THEN   
       
     num(3)=num(8)+1
     num(4)=num(7)+1
     num(5)=num(4)+ny2*2+1
     num(6)=num(5)+1
     num(7)=num(6)+1
     num(8)=num(4)+2
     num(9)=num(4)+1
     num(1)=num(3)+2
     num(2)=num(3)+1
     
      facx=w2/nx2
     facy=h2/ny2
     
      IF(facx>zero)THEN
        midllex=facx/2.0_iwp
      ELSE
        midllex=zero  
      END IF
     
      IF(facy>zero)THEN
        midlley=facy/2.0_iwp
      ELSE
        midlley=zero  
      END IF
     
         coord(1:3,1)=(nx1)*facb+(colx-1-nx1)*facx
     coord(4,1)=coord(3,1)+midllex
     coord(5,1)=coord(3,1)+facx
     coord(6,1)=coord(5,1)
     coord(7,1)=coord(6,1)
     coord(8,1)=coord(1,1)+midllex
     coord(9,1)=coord(4,1)
   
     coord(3:5,2)=0.0_iwp-((rowy-1)*h1/ny1)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
    ELSE IF(colx>=nx1+1.and.rowy>ny1+1)THEN   
       
     num(3)=num(1)
     num(4)=num(8)
     num(5)=num(7)
     num(6)=num(5)+1
     num(7)=num(6)+1
     num(8)=num(4)+2
     num(9)=num(4)+1
     num(1)=num(3)+2
     num(2)=num(3)+1
     
      facx=w2/nx2
     facy=h2/ny2
     
      IF(facx>zero)THEN
        midllex=facx/2.0_iwp
      ELSE
        midllex=zero  
      END IF
     
      IF(facy>zero)THEN
        midlley=facy/2.0_iwp
      ELSE
        midlley=zero  
      END IF
     
         coord(1:3,1)=(nx1)*facb+(colx-1-nx1)*facx
     coord(4,1)=coord(3,1)+midllex
     coord(5,1)=coord(3,1)+facx
     coord(6,1)=coord(5,1)
     coord(7,1)=coord(6,1)
     coord(8,1)=coord(1,1)+midllex
     coord(9,1)=coord(4,1)
   
     !coord(3:5,2)=0.0_iwp-((rowy-1)*h1/ny1)
       coord(3,2)=coord(1,2)
     coord(4,2)=coord(8,2)
     coord(5,2)=coord(7,2)
     coord(6,2)=coord(5,2)-midlley
     coord(7,2)=coord(6,2)-midlley
     coord(1,2)=coord(3,2)-facy
     coord(2,2)=coord(3,2)-midlley
     coord(8,2)=coord(1,2)
     coord(9,2)=coord(2,2)
     
  END IF
     
 
RETURN
END SUBROUTINE emb_2d_geom_9n


SUBROUTINE mesh_size(element,nod,nels,nn,nxe,nye,nze)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE
 CHARACTER(LEN=15),INTENT(IN)::element
 INTEGER,INTENT(IN)::nod,nxe,nye
 INTEGER,INTENT(IN),OPTIONAL::nze
 INTEGER,INTENT(OUT)::nels,nn
 IF(element=="triangle")THEN
   nels=nxe*nye*2
   IF(nod==3)nn=(nxe+1)*(nye+1)
   IF(nod==6)nn=(2*nxe+1)*(2*nye+1)
   IF(nod==10)nn=(3*nxe+1)*(3*nye+1)
   IF(nod==15)nn=(4*nxe+1)*(4*nye+1)
 ELSE IF(element=="quadrilateral")THEN
   nels=nxe*nye
   IF(nod==4)nn=(nxe+1)*(nye+1)
   IF(nod==5)nn=(nxe+1)*(nye+1)+nxe*nye
   IF(nod==8)nn=(2*nxe+1)*(nye+1)+(nxe+1)*nye
   IF(nod==9)nn=(2*nxe+1)*(2*nye+1)
 ELSE IF(element=="hexahedron")THEN
   nels=nxe*nye*nze
   IF(nod==8)nn=(nxe+1)*(nye+1)*(nze+1)
   IF(nod==14)nn=4*nxe*nye*nze+2*(nxe*nye+nye*nze+nze*nxe)+nxe+nye+nze+1
   IF(nod==20)nn=((2*nxe+1)*(nze+1)+(nxe+1)*nze)*(nye+1)+                 &
     (nxe+1)*(nze+1)*nye
 END IF
RETURN
END SUBROUTINE mesh_size

SUBROUTINE mesh_size_slope9(nx1,nx2,ny1,ny2,nn,nels)
!
!  This subroutine returns the number of elements (nels) and the number
!  of nodes (nn) in a 2-d geometry-created mesh.
!
 IMPLICIT NONE

 INTEGER,INTENT(IN)::nx1,nx2,ny1,ny2
 INTEGER,INTENT(OUT)::nels,nn
 
 nels=nx1*(ny1+ny2)+(nx2*ny2)
 nn=((ny1+ny2)*2+1)*(nx1*2+1)+(ny2*2+1)*(nx2*2)
 
 
 
RETURN
END SUBROUTINE mesh_size_slope9

SUBROUTINE geom_rect(iel,x_coords,y_coords,coord,num,dir)
!
! This subroutine forms the coordinates and connectivity for a
! rectangular mesh of rt. angled triangular elements (3, 6, 10 or 15-node)
! or quadrilateral elements (4, 8 or 9-node) counting in the
! x- or y-dir. 
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN)::x_coords(:),y_coords(:)
 REAL(iwp),INTENT(OUT)::coord(:,:)
 CHARACTER(LEN=15)::element
 CHARACTER(LEN=1),INTENT(IN)::dir
 INTEGER,INTENT(IN)::iel
 INTEGER,INTENT(OUT)::num(:)
 INTEGER::ip,iq,jel,fac1,nod,nxe,nye
 REAL(iwp)::pt5=0.5_iwp,two=2.0_iwp,d3=3.0_iwp 
 nxe=UBOUND(x_coords,1)-1
 nod=UBOUND(num,1)
 element='square'
 IF(element=='triangle')THEN
  CONTINUE
 ELSE
   nye=UBOUND(y_coords,1)-1
   IF(dir=='x'.OR.dir=='r')THEN
     iq=(iel-1)/nxe+1
     ip=iel-(iq-1)*nxe
   ELSE
     ip=(iel-1)/nye+1
     iq=iel-(ip-1)*nye
   END IF
   SELECT CASE(nod)
   CASE(4)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(nxe+1)+ip		        		
       num(2)=(iq-1)*(nxe+1)+ip				
       num(3)=num(2)+1					
       num(4)=num(1)+1					
     ELSE
       num(1)=(ip-1)*(nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(nye+1)+iq
       num(4)=num(3)+1
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
   CASE(5)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(2*nxe+1)+ip
       num(2)=(iq-1)*(2*nxe+1)+ip
       num(3)=num(2)+1
       num(4)=num(1)+1
       num(5)=iq*(2*nxe+1)+ip-nxe
     ELSE
       num(1)=(ip-1)*(2*nye+1)+iq+1
       num(2)=num(1)-1
       num(3)=ip*(2*nye+1)+iq
       num(4)=num(3)+1
       num(5)=ip*(2*nye+1)+iq-nye
     END IF
!
     coord(1:2,1)=x_coords(ip)
     coord(3:4,1)=x_coords(ip+1)
     coord(1,2)=y_coords(iq+1)
     coord(2:3,2)=y_coords(iq)
     coord(4,2)=coord(1,2)
     coord(5,:)=0.25_iwp*(coord(1,:)+coord(2,:)+coord(3,:)+coord(4,:))
   CASE(8)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(3*nxe+2)+2*ip-1                 
       num(2)=iq*(3*nxe+2)+ip-nxe-1		  
       num(3)=(iq-1)*(3*nxe+2)+2*ip-1		   
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+1
       num(7)=num(1)+2
       num(8)=num(1)+1
     ELSE
       num(1)=(ip-1)*(3*nye+2)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*(3*nye+2)+2*nye+iq+1
       num(5)=ip*(3*nye+2)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
   CASE(9)
     IF(dir=='x'.OR.dir=='r')THEN
       num(1)=iq*(4*nxe+2)+2*ip-1
       num(1)=iq*(2*nxe+1)*2+2*ip-1
       num(2)=iq*(4*nxe+2)+2*ip-nxe-4
       num(2)=(2*nxe+1)+(ip-1)*2+1+(iq-1)*(2*nxe+1)*2
       num(3)= (iq-1)*(4*nxe+2)+2*ip-1
       num(4)=num(3)+1
       num(5)=num(4)+1
       num(6)=num(2)+2
       num(7)=num(1)+2
       num(8)=num(1)+1
       num(9)=num(2)+1
     ELSE
       num(1)=(ip-1)*2*(2*nye+1)+2*iq+1
       num(2)=num(1)-1
       num(3)=num(1)-2
       num(4)=(ip-1)*2*(2*nye+1)+2*nye+2*iq
       num(5)=ip*2*(2*nye+1)+2*iq-1
       num(6)=num(5)+1
       num(7)=num(5)+2
       num(8)=num(4)+2
       num(9)=num(4)+1
     END IF
!
     coord(1:3,1)=x_coords(ip)
     coord(5:7,1)=x_coords(ip+1)
     coord(4,1)=pt5*(coord(3,1)+coord(5,1))
     coord(8,1)=pt5*(coord(7,1)+coord(1,1))
     coord(1,2)=y_coords(iq+1)
     coord(7:8,2)=y_coords(iq+1)
     coord(3:5,2)=y_coords(iq)
     coord(2,2)=pt5*(coord(1,2)+coord(3,2))
     coord(6,2)=pt5*(coord(5,2)+coord(7,2))
     coord(9,:)=pt5*(coord(4,:)+coord(8,:))
   CASE DEFAULT
     WRITE(11,'(a)')"Wrong number of nodes for quadrilateral element"
     STOP
   END SELECT
 END IF
RETURN
END SUBROUTINE geom_rect


end module geom
