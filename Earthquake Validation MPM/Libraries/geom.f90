module geom
    save
    contains
    
SUBROUTINE mesh_size_slope9(nx1,nx2,ny1,ny2,nn,nels)
    !
    !  This subroutine returns the number of elements (nels) and the number
    !  of nodes (nn) in a 2-d geometry-created mesh.
    !
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nx1,nx2,ny1,ny2
    INTEGER,INTENT(OUT):: nels,nn
    
    nels = nx1*(ny1+ny2)+(nx2*ny2)
    nn   = ((ny1+ny2)*2+1)*(nx1*2+1)+(ny2*2+1)*(nx2*2)
    RETURN
END SUBROUTINE mesh_size_slope9


SUBROUTINE mesh_size_slope4(nx1,nx2,ny1,ny2,nn,nels)
  !
  !  This subroutine returns the number of elements (nels) and the number
  !  of nodes (nn) in a 2-d geometry-created mesh.
  !
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nx1,nx2,ny1,ny2
  INTEGER,INTENT(OUT):: nels,nn
  
  nels = nx1*(ny1+ny2)+(nx2*ny2)
  nn   = ((ny1+ny2)+1)*(nx1+1)+(ny2+1)*(nx2)
  RETURN
END SUBROUTINE mesh_size_slope4


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


SUBROUTINE emb_2d_geom_4n(iel,nx1,nx2,ny1,ny2,w1,s1,w2,h1,h2,coord,num)
  !
  ! This subroutine forms the nodal coordinates and numbering for a 2-d
  ! slope of 4-node quadrilaterals. Nodes numbering in the y-direction,
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
    ip=(iel-1)/nye+1
    iq=iel-(ip-1)*nye
   ELSE
     iq=(iel-nt-1)/nxe+ny1+1
     ip=iel-nt-(iq-ny1-1)*nxe
   END IF
   IF(ip<=nx1)THEN
     num(1)=(ip-1)*(nye+1)+iq+1
     num(2)=num(1)-1
     num(3)=ip*(nye+1)+iq
     num(4)=num(3)+1

     IF(iq<=ny1)THEN
       coord(1,1)=(ip-one)*(w1+iq*facx)/nx1
       coord(2,1)=(ip-one)*(w1+(iq-1)*facx)/nx1
       coord(3,1)=ip*(w1+(iq-1)*facx)/nx1
       coord(4,1)=ip*(w1+iq*facx)/nx1
       coord(1,2)=-iq*facy
       coord(2,2)=-(iq-1)*facy
       coord(3,2)=-(iq-1)*facy
       coord(4,2)=-iq*facy
     END IF
   END IF
  RETURN
END SUBROUTINE emb_2d_geom_4n

end module geom