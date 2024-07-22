module mpm_geom
    save
    contains

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
            coord(3,2)=-(iq-1)*facy
            coord(4,2)=-iq*facy
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
END SUBROUTINE emb_2d_geom


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
    
      IF(column==slope1.and.slopeopt==1) coord(3,2)=coord(2,2)+(coord(1,2)-coord(2,2))/2.0_iwp
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




end module mpm_geom