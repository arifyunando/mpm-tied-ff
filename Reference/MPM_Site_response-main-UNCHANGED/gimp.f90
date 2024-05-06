
 MODULE GIMP
 SAVE 
 CONTAINS

 !-This library contains all subroutins regarding GIMP method

    SUBROUTINE GIMP_activenode(g_num,nip,g_coord,gm_coord,lp_mp,gimptol,neighb,a_ele,nf,GIMP_node_mp)
     !     
     ! Subroutine to activate al nodes inside the support domain of all material points
     !
    IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::nip,g_num(:,:),neighb(:,:),a_ele(:)
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),lp_mp(:,:),gimptol
   INTEGER,INTENT(INOUT)::nf(:,:),GIMP_node_mp(:,:)
   REAL(iwp),ALLOCATABLE::supdom(:,:) !supdom = particle suport domain
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,maxx,miny,maxy,zero=0.0_iwp,cellsize,tol,lpx,lpy
   INTEGER::i,j,k,neighbour,nn,nmps,s,iel,neig
   !LOGICAL,INTENT(IN)::gimptol

   ALLOCATE(supdom(4,2))
   cellsize=g_coord(2,1)-g_coord(2,2)
   !IF(nip==1)twolp=cellsize/(nip)
   !IF(nip==4)twolp=cellsize/(nip/two)
   !lp=twolp/two
   GIMP_node_mp=0

   nn=UBOUND(g_coord,2)
   nmps=UBOUND(gm_coord,2)
   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   tol=gimptol
   MP_Loop:DO s=1,nmps
     lpx=lp_mp(1,s);lpy=lp_mp(2,s)
     iel=a_ele(s)
     !-Evaluation of the coordinates of the boundaries of the support domain for the GIMP
     supdom(1,1)=gm_coord(1,s)-(cellsize+lpx) !Left boundarie of the support domain
     supdom(1,2)=gm_coord(2,s)+(cellsize+lpy) !Upper boundarie of the support domain
     supdom(2,1)=gm_coord(1,s)+(cellsize+lpx) !Right boundarie of the support domain
     supdom(4,2)=gm_coord(2,s)-(cellsize+lpy) !Lower boundarie of the support domain
     !-Tolerance is included to be sure the nodes are considered inside a support domain
     !-when they are just over the the edge of the boudarie
     minx=supdom(1,1)+tol
     maxx=supdom(2,1)-tol
     miny=supdom(4,2)+tol
     maxy=supdom(1,2)-tol

     nf(:,g_num(:,iel))=1

   DO j=1,8
     neig=neighb(iel,j)
   
     DO i=1,4
      IF(g_coord(1,g_num(i,neig))>=minx.and.g_coord(1,g_num(i,neig))<=maxx)THEN ! x coordenate of the node is inside the domain
        IF(g_coord(2,g_num(i,neig))<=maxy.and.g_coord(2,g_num(i,neig))>=miny)THEN  ! y coordenate of the node is inside the domain      
         nf(:,g_num(i,neig))=1 
         GIMP_node_mp(s,g_num(i,neig))=1
        END IF
      END IF  
     END DO  

   END DO

!$$$$$$    j=j
!$$$$$$   DO j=1,size(nf,2)
!$$$$$$    IF(nf(1,j)==1)GIMP_node_mp(s,j)=1
!$$$$$$   END DO  

  END DO MP_Loop


20 RETURN
   END SUBROUTINE GIMP_activenode

   SUBROUTINE GIMP_nodsup(s,g_num,nip,g_coord,gm_coord,lp_mp,mpoints,valuesg,gimptol,neighb,nny,a_ele,GIMP_nodes)
     !     
     ! Subroutine to save al the nodes (number of the nodes) inside the support domain for each material point with GIMP 
     ! node numbering should be in y direction
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::nip,s,g_num(:,:),neighb(:,:),a_ele(:),nny
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),mpoints(:,:),lp_mp(:,:),gimptol
   INTEGER,INTENT(IN OUT)::GIMP_nodes(:,:),valuesg(:)
   REAL(iwp),ALLOCATABLE::supdom(:,:) !supdom = particle suport domain
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,maxx,miny,maxy,zero=0.0_iwp,cellsize,tol,lpx,lpy
   INTEGER::i,j,k,n,m,neighbour,nn,nmps,nod,neig,iel,mainnod
   !LOGICAL,INTENT(IN)::gimptol

   ALLOCATE(supdom(4,2))
   !GIMP_nodes=zero
   cellsize=g_coord(2,1)-g_coord(2,2)
   !IF(nip==1)twolp=cellsize/(nip)
   !IF(nip==4)twolp=cellsize/(nip/two)
   !lp=twolp/two
   lpx=lp_mp(1,s)!+gimptol
   lpy=lp_mp(2,s)!+gimptol

   nn=UBOUND(g_coord,2)
   nmps=UBOUND(gm_coord,2)
   !tol=1.0e-8

   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   tol=gimptol
   
     !-Evaluation of the coordinates of the boundaries of the support domain for the GIMP
     supdom(1,1)=gm_coord(1,s)-(cellsize+lpx) !Left boundarie of the support domain
     supdom(1,2)=gm_coord(2,s)+(cellsize+lpy) !Upper boundarie of the support domain
     supdom(2,1)=gm_coord(1,s)+(cellsize+lpx) !Right boundarie of the support domain
     supdom(4,2)=gm_coord(2,s)-(cellsize+lpy) !Lower boundarie of the support domain
     !-Tolerance is included to be sure the nodes are considered inside a support domain
     !-when they are just over the the edge of the boudarie
     minx=supdom(1,1)+tol  ! The tolerance is making the supoprt domain smaller
     maxx=supdom(2,1)-tol
     miny=supdom(4,2)+tol
     maxy=supdom(1,2)-tol
     
     k=0
     n=0
     m=0
     j=1
     valuesg(s)=0
     iel=a_ele(s)
     neig=neighb(iel,1)
     nod=g_num(2,neig)
     mainnod=g_num(2,neig)
     valuesg(s)=0

    DO WHILE(m<16) !-Maximum number of nodes inside the support domain of a particle
      n=n+1;m=m+1
      IF(g_coord(1,nod)>=minx.and.g_coord(1,nod)<=maxx)THEN ! x coordenate of the node is inside the domain
        IF(g_coord(2,nod)<=maxy.and.g_coord(2,nod)>=miny)THEN  ! y coordenate of the node is inside the domain
          k=k+1 
          valuesg(s)=valuesg(s)+1
          GIMP_nodes(k,s)=nod  
        END IF
      END IF 
          IF(n==4)THEN
             n=0
             nod=mainnod+nny*j
             j=j+1
          ELSE
             nod=nod+1
          END IF
    END DO 

!$$$$$$       IF(m==9)GOTO 20 
     
20 RETURN
   END SUBROUTINE GIMP_nodsup

   SUBROUTINE GIMP_nodsup2(nmps,g_num,nip,g_coord,gm_coord,mpoints,valuesg,gimptol,neighb,nny,a_ele,GIMP_nodes,  &
   							GIMP_node_mp_aux,nf)
     !     
     ! Subroutine to save al the nodes (number of the nodes) inside the support domain for each material point with GIMP 
     ! node numbering should be in y direction
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::nip,nmps,g_num(:,:),neighb(:,:),a_ele(:),nny
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),mpoints(:,:),gimptol
   INTEGER,INTENT(IN OUT)::GIMP_nodes(:,:),valuesg(:),GIMP_node_mp_aux(:,:),nf(:,:)
   REAL(iwp),ALLOCATABLE::supdom(:,:) !supdom = particle suport domain
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,maxx,miny,maxy,zero=0.0_iwp,cellsize,tol
   INTEGER::i,j,k,n,m,neighbour,nn,nod,neig,iel,mainnod,s
   !LOGICAL,INTENT(IN)::gimptol

   ALLOCATE(supdom(4,2))
   !GIMP_nodes=zero
   cellsize=g_coord(2,1)-g_coord(2,2)
   IF(nip==1)twolp=cellsize/(nip)
   IF(nip==4)twolp=cellsize/(nip/two)
   lp=twolp/two

   nn=UBOUND(g_coord,2)
!$$$$$$    nmps=UBOUND(gm_coord,2)
   tol=1.0e-8

   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol

   
  DO s=1,nmps 
     !-Evaluation of the coordinates of the boundaries of the support domain for the GIMP
     supdom(1,1)=gm_coord(1,s)-(cellsize+lp) !Left boundarie of the support domain
     supdom(1,2)=gm_coord(2,s)+(cellsize+lp) !Upper boundarie of the support domain
     supdom(2,1)=gm_coord(1,s)+(cellsize+lp) !Right boundarie of the support domain
     supdom(4,2)=gm_coord(2,s)-(cellsize+lp) !Lower boundarie of the support domain
     !-Tolerance is included to be sure the nodes are considered inside a support domain
     !-when they are just over the the edge of the boudarie
     minx=supdom(1,1)+tol  ! The tolerance is making the supoprt domain smaller
     maxx=supdom(2,1)-tol
     miny=supdom(4,2)+tol
     maxy=supdom(1,2)-tol
     
     k=0
     n=0
     m=0
     j=1
     valuesg(s)=0
     iel=a_ele(s)
     neig=neighb(iel,1)
     nod=g_num(2,neig)
     mainnod=g_num(2,neig)
     valuesg(s)=0
 
    DO WHILE(m<16) !-Maximum number of nodes inside the support domain of a particle
      n=n+1;m=m+1
      IF(g_coord(1,nod)>=minx.and.g_coord(1,nod)<=maxx)THEN ! x coordenate of the node is inside the domain
        IF(g_coord(2,nod)<=maxy.and.g_coord(2,nod)>=miny)THEN  ! y coordenate of the node is inside the domain
          k=k+1 
          valuesg(s)=valuesg(s)+1
          GIMP_nodes(k,s)=nod  
          GIMP_node_mp_aux(s,nod)=1
          nf(:,nod)=1
        END IF
      END IF 
          IF(n==4)THEN
             n=0
             nod=mainnod+nny*j
             j=j+1
          ELSE
             nod=nod+1
          END IF
    END DO 

   END DO
!$$$$$$       IF(m==9)GOTO 20 
     
20 RETURN
   END SUBROUTINE GIMP_nodsup2

 
  SUBROUTINE GIMP_values(GIMP_nodes,s,bound,values)
      !     
      ! Subroutine to compute the amount of nodes inside a material point support domain 
      !
     IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::GIMP_nodes(:,:),s,bound
    INTEGER,INTENT(OUT)::values
    INTEGER::i

    values=0

      DO i=1,9
        IF(GIMP_nodes(i,s)>0)values=values+1
      END DO  
 
  END SUBROUTINE GIMP_values

  SUBROUTINE GIMP_funder(s,nip,g_coord,cellsize,gm_coord,lp_mp,GIMP_nodes,gimptol,derGIMP,funGIMP)
     !     
     ! Subroutine to compute GIMP shape functions and GIMP derivatives of the shape functions 
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::GIMP_nodes(:,:),s,nip
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),cellsize,lp_mp(:,:),gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,miny,tol
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,Sr,Sz,dSr,dSz,fact,lpx,lpy
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   INTEGER::i,j,k,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol


   derGIMP=0.0
   funGIMP=0.0
   Sr=0.0
   Sz=0.0
   dSr=0.0
   dSz=0.0
   !tol=1.0e-8
   nodes=UBOUND(GIMP_nodes,1)
   !IF(nip==1)twolp=cellsize/(nip)
   !IF(nip==4)twolp=cellsize/(nip/two)
   !lp=twolp/two
   lpx=lp_mp(1,s);lpy=lp_mp(2,s)
   i=1
   k=1

   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol

   DO WHILE(i<=nodes)   ! 9 is the maximum of nodes inside a material point suport domain 
    IF(GIMP_nodes(i,s)>0)THEN
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(i,s)     !--node 'i' in the suport domain to interpolate value to the current material point 's' 
     minx=g_coord(1,nod)     !--minx coordinate of node 'i'
     xpoint=gm_coord(1,s)    !--xpoint coordinate of current material point 's'
     xdist=xpoint-minx       !--xdist distance between node 'i' and material point 's' (if material point is at the right side of the node, distance is positive)
     miny=g_coord(2,nod)     !--miny coordinate of node 'i'
     ypoint=gm_coord(2,s)    !--ypoint coordinate of material point 's'
     ydist=(-ypoint)-(-miny) !--ydist distance between node 'i' and material point 's' (if material point is below the node, distance is positive)
     
     !--GIMP shape functions and derivatives in the x direction
     IF(g_coord(1,nod)<=gm_coord(1,s))THEN
        fact=xdist
        lp=lpx
        IF(fact<=lpx)THEN    
         Sr=one-((fact**2+lpx**2)/(two*cellsize*lpx))
         dSr=-(fact)/(cellsize*lp)
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sr=one-(fact/cellsize)
         dSr=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<=(cellsize+lp-tol))THEN
         Sr=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSr=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(1,nod)>=gm_coord(1,s))THEN  
        fact=xdist 
         lp=lpx
        IF(fact>=(-lp))THEN    
          Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sr=one+(fact/cellsize)
          dSr=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>=(-(cellsize+lp-tol)))THEN
         Sr=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSr=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF     
     
      !--GIMP shape functions and derivatives in the y direction
     IF(g_coord(2,nod)>=gm_coord(2,s))THEN
        fact=ydist  !--ydist positive
         lp=lpy
        IF(fact<=lp)THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sz=one-(fact/cellsize)
         dSz=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<=(cellsize+lp-tol))THEN
         Sz=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSz=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(2,nod)<=gm_coord(2,s))THEN  
        fact=ydist !--ydist negative
         lp=lpy
        IF(fact>=(-lp))THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sz=one+(fact/cellsize)
         dSz=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>=(-(cellsize+lp-tol)))THEN
         Sz=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSz=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF  

     funGIMP(k)=Sr*Sz
     derGIMP(1,k)=dSr*Sz
     derGIMP(2,k)=-dSz*Sr

     IF(Sr<=0.0.or.Sz<=0)THEN
       derGIMP(1,k)=0.0
       derGIMP(2,k)=0.0
     END IF

     k=k+1
    END IF
    i=i+1
   END DO 
   RETURN
  END SUBROUTINE GIMP_funder

    SUBROUTINE GIMP_funder3(s,nip,g_coord,cellsize,gm_coord,GIMP_nodes,gimptol,derGIMP,funGIMP)
     !     
     ! Subroutine to compute GIMP shape functions and GIMP derivatives of the shape functions 
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::GIMP_nodes(:,:),s,nip
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),cellsize,gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,miny,tol
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,Sr,Sz,dSr,dSz,fact
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   INTEGER::i,j,k,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol

   funGIMP=0.0  
   derGIMP=0.0
   Sr=0.0
   Sz=0.0
   dSr=0.0
   dSz=0.0
   !tol=1.0e-8
   nodes=UBOUND(GIMP_nodes,1)
   IF(nip==1)twolp=cellsize/(nip)
   IF(nip==4)twolp=cellsize/(nip/two)
   lp=twolp/two
   i=1
   k=1

   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol

   DO WHILE(i<=nodes)   ! 9 is the maximum of nodes inside a material point suport domain 
    IF(GIMP_nodes(i,s)>0)THEN
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(i,s)     !--node 'i' in the suport domain to interpolate value to the current material point 's' 
     minx=g_coord(1,nod)     !--minx coordinate of node 'i'
     xpoint=gm_coord(1,s)    !--xpoint coordinate of current material point 's'
     xdist=xpoint-minx       !--xdist distance between node 'i' and material point 's' (if material point is at the right side of the node, distance is positive)
     miny=g_coord(2,nod)     !--miny coordinate of node 'i'
     ypoint=gm_coord(2,s)    !--ypoint coordinate of material point 's'
     ydist=(-ypoint)-(-miny) !--ydist distance between node 'i' and material point 's' (if material point is below the node, distance is positive)
     
     !--GIMP shape functions and derivatives in the x direction
     IF(g_coord(1,nod)<=gm_coord(1,s))THEN
        fact=xdist
        IF(fact<=lp)THEN    
         Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact)/(cellsize*lp)
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sr=one-(fact/cellsize)
         dSr=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<(cellsize+lp-tol))THEN
         Sr=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSr=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(1,nod)>=gm_coord(1,s))THEN  
        fact=xdist 
        IF(fact>=(-lp))THEN    
          Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sr=one+(fact/cellsize)
          dSr=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>(-(cellsize+lp-tol)))THEN
         Sr=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSr=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF     
     
      !--GIMP shape functions and derivatives in the y direction
     IF(g_coord(2,nod)>=gm_coord(2,s))THEN
        fact=ydist  !--ydist positive
        IF(fact<=lp)THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sz=one-(fact/cellsize)
         dSz=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<(cellsize+lp-tol))THEN
         Sz=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSz=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(2,nod)<=gm_coord(2,s))THEN
        fact=ydist !--ydist negative
        IF(fact>=(-lp))THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sz=one+(fact/cellsize)
         dSz=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>(-(cellsize+lp-tol)))THEN
         Sz=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSz=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF  

     funGIMP(k)=Sr*Sz
     derGIMP(1,k)=dSr*Sz
     derGIMP(2,k)=-dSz*Sr

     IF(Sr<=0.0.or.Sz<=0)THEN
       derGIMP(1,k)=0.0
       derGIMP(2,k)=0.0
     END IF

!$$$$$$     IF((derGIMP(1,k)>0.0.or.derGIMP(2,k)>0.0).and.i==1)THEN
!$$$$$$       derGIMP(1,k)=0.0
!$$$$$$       derGIMP(2,k)=0.0
!$$$$$$       funGIMP(k)=0.0
!$$$$$$     ELSE IF((derGIMP(1,k)>0.0.or.derGIMP(2,k)<0.0).and.i==2)THEN
!$$$$$$     derGIMP(1,k)=0.0
!$$$$$$     derGIMP(2,k)=0.0  
!$$$$$$     funGIMP(k)=0.0
!$$$$$$     ELSE IF((derGIMP(1,k)<0.0.or.derGIMP(2,k)<0.0).and.i==3)THEN
!$$$$$$     derGIMP(1,k)=0.0
!$$$$$$     derGIMP(2,k)=0.0
!$$$$$$     funGIMP(k)=0.0
!$$$$$$     ELSE IF((derGIMP(1,k)<0.0.or.derGIMP(2,k)>0.0).and.i==4)THEN
!$$$$$$     derGIMP(1,k)=0.0
!$$$$$$     derGIMP(2,k)=0.0
!$$$$$$     funGIMP(k)=0.0
!$$$$$$     END IF
     
     k=k+1
    END IF
    i=i+1
   END DO 
   RETURN
  END SUBROUTINE GIMP_funder3

  SUBROUTINE GIMP_funder4(s,nip,g_coord,cellsize,gm_coord,GIMP_nodes,gimptol,derGIMP,funGIMP)
     !     
     ! Subroutine to compute GIMP shape functions and GIMP derivatives of the shape functions 
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::GIMP_nodes(:,:),s,nip
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),cellsize,gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,miny,tol
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,Sr,Sz,dSr,dSz,fact
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   INTEGER::i,j,k,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol

   funGIMP=0.0  
   derGIMP=0.0
   Sr=0.0
   Sz=0.0
   dSr=0.0
   dSz=0.0
   !tol=1.0e-8
   nodes=UBOUND(GIMP_nodes,1)
   IF(nip==1)twolp=cellsize/(nip)
   IF(nip==4)twolp=cellsize/(nip/two)
   lp=twolp/two
   i=1
   k=1
   !
   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol

   DO WHILE(i<=nodes)   ! 9 is the maximum of nodes inside a material point suport domain 
    IF(GIMP_nodes(i,s)>0)THEN
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(i,s)     !--node 'i' in the suport domain to interpolate value to the current material point 's' 
     minx=g_coord(1,nod)     !--minx coordinate of node 'i'
     xpoint=gm_coord(1,s)    !--xpoint coordinate of current material point 's'
     xdist=xpoint-minx       !--xdist distance between node 'i' and material point 's' (if material point is at the right side of the node, distance is positive)
     miny=g_coord(2,nod)     !--miny coordinate of node 'i'
     ypoint=gm_coord(2,s)    !--ypoint coordinate of material point 's'
     ydist=(-ypoint)-(-miny) !--ydist distance between node 'i' and material point 's' (if material point is below the node, distance is positive)
     
     !--GIMP shape functions and derivatives in the x direction
     IF(g_coord(1,nod)<=gm_coord(1,s))THEN
        fact=xdist
        IF(fact<=lp)THEN    
         Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact)/(cellsize*lp)
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sr=one-(fact/cellsize)
         dSr=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<(cellsize+lp-tol))THEN
         Sr=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSr=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(1,nod)>=gm_coord(1,s))THEN  
        fact=xdist 
        IF(fact>=(-lp))THEN    
          Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sr=one+(fact/cellsize)
          dSr=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>(-(cellsize+lp-tol)))THEN
         Sr=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSr=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF     
     
      !--GIMP shape functions and derivatives in the y direction
     IF(g_coord(2,nod)>=gm_coord(2,s))THEN
        fact=ydist  !--ydist positive
        IF(fact<=lp)THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sz=one-(fact/cellsize)
         dSz=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<(cellsize+lp-tol))THEN
         Sz=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSz=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(2,nod)<=gm_coord(2,s))THEN
        fact=ydist !--ydist negative
        IF(fact>=(-lp))THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sz=one+(fact/cellsize)
         dSz=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>(-(cellsize+lp-tol)))THEN
         Sz=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSz=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF  

     funGIMP(k)=Sr*Sz
     derGIMP(1,k)=dSr*Sz
     derGIMP(2,k)=-dSz*Sr

     IF(Sr<=0.0.or.Sz<=0)THEN
       derGIMP(1,k)=0.0
       derGIMP(2,k)=0.0
     END IF

     IF((derGIMP(1,k)>0.0.or.derGIMP(2,k)>0.0).and.i==1)THEN
       derGIMP(1,k)=0.0
       derGIMP(2,k)=0.0
       funGIMP(k)=0.0
     ELSE IF((derGIMP(1,k)>0.0.or.derGIMP(2,k)<0.0).and.i==2)THEN
     derGIMP(1,k)=0.0
     derGIMP(2,k)=0.0  
     funGIMP(k)=0.0
     ELSE IF((derGIMP(1,k)<0.0.or.derGIMP(2,k)<0.0).and.i==3)THEN
     derGIMP(1,k)=0.0
     derGIMP(2,k)=0.0
     funGIMP(k)=0.0
     ELSE IF((derGIMP(1,k)<0.0.or.derGIMP(2,k)>0.0).and.i==4)THEN
     derGIMP(1,k)=0.0
     derGIMP(2,k)=0.0
     funGIMP(k)=0.0
     END IF
     
     k=k+1
    END IF
    i=i+1
   END DO 
   RETURN
  END SUBROUTINE GIMP_funder4

  SUBROUTINE GIMP_funder2(s,nip,g_coord,cellsize,gm_coord,lp_mp,GIMP_nodes,gimptol,derGIMP,funGIMP)
                        
     !     
     ! Subroutine to compute GIMP shape functions and GIMP derivatives of the shape functions 
     ! to send information to the nodes (this subroutine is diferent from GIMP_funder only
     ! because the values of the functions and derivatives are stored in a different way)
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::GIMP_nodes(:,:),s,nip
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),cellsize,lp_mp(:,:),gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:,:)
   REAL(iwp)::twolp,two=2.0_iwp,minx,miny,Sr,Sz,dSr,dSz,lp,fact,lpx,lpy
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,tol
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   INTEGER::i,j,k,m,nod,side
   !LOGICAL,INTENT(IN)::gimptol
   
   funGIMP=0.0
   derGIMP=0.0
   Sr=0.0
   Sz=0.0
   dSr=0.0
   dSz=0.0
   !tol=1.0e-8  !-tolerance values is to avoid the effect of nodes over the support domain boundarie of particles, that can give errors
   !IF(nip==1)twolp=cellsize/(nip)
   !IF(nip==4)twolp=cellsize/(nip/two)
   !lp=twolp/two
   lpx=lp_mp(1,s)!-gimptol
   lpy=lp_mp(2,s)!-gimptol
   i=1
   k=1
   m=1

   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol
   

   DO WHILE(i<=9)   ! 9 is the maximum of nodes inside a material point suport domain 
    IF(GIMP_nodes(i,s)>0)THEN
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(i,s)     !--node 'i' in the suport domain to interpolate value to the current material point 's' 
     minx=g_coord(1,nod)     !--x coordinate of node 'i'
     xpoint=gm_coord(1,s)    !--x coordinate of current material point 's'
     xdist=xpoint-minx       !--x distance between node 'i' and material point 's' (if material point is at the right side of the node, distance is positive)
     miny=g_coord(2,nod)     !--y coordinate of node 'i'
     ypoint=gm_coord(2,s)    !--y coordinate of material point 's'
     ydist=(-ypoint)-(-miny) !--y distance between node 'i' and material point 's' (if material point is below the node, distance is positive)
     
     !--GIMP shape functions and derivatives in the x direction
     IF(g_coord(1,nod)<=gm_coord(1,s))THEN
        fact=xdist !--particle at te right side
        IF(fact<=lpx)THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sr=one-((fact**2+lpx**2)/(two*cellsize*lpx))
         dSr=-(fact)/(cellsize*lpx)
        ELSE IF(fact>=lpx.and.fact<=(cellsize-lpx))THEN
         Sr=one-(fact/cellsize)
         dSr=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lpx).and.fact<=(cellsize+lpx-tol))THEN
         Sr=((cellsize+(lpx-tol)-fact)**2)/(four*cellsize*(lpx-tol))
         dSr=-((cellsize+(lpx-tol)-fact)/(two*(lpx-tol)*cellsize))
        END IF     
     ELSE IF(g_coord(1,nod)>=gm_coord(1,s))THEN  
        fact=xdist !--particle at the left side
        IF(fact>=(-lpx))THEN    
          Sr=one-((fact**2+lpx**2)/(two*cellsize*lpx))
         dSr=-(fact/(cellsize*lpx))
        ELSE IF(fact<=(-lpx).and.fact>=(-(cellsize-lpx)))THEN
         Sr=one+(fact/cellsize)
          dSr=one/cellsize
        ELSE IF(fact<=(-(cellsize-lpx)).and.fact>(-(cellsize+lpx-tol)))THEN
         Sr=((cellsize+(lpx-tol)+fact)**2)/(four*cellsize*(lpx-tol))
         dSr=((cellsize+(lpx-tol)+fact)/(two*(lpx-tol)*cellsize))
        END IF  
     END IF     
     
      !--GIMP shape functions and derivatives in the y direction
     IF(g_coord(2,nod)>=gm_coord(2,s))THEN
        fact=ydist  !--ydist positive (particle below the node)
        IF(fact<=lpy)THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
          Sz=one-((fact**2+lpy**2)/(two*cellsize*lpy))
          dSz=-(fact/(cellsize*lpy))
         ELSE IF(fact>=lpy.and.fact<=(cellsize-lpy))THEN
          Sz=one-(fact/cellsize)
          dSz=-(one/cellsize)
         ELSE IF(fact>=(cellsize-lpy).and.fact<=(cellsize+lpy-tol))THEN
          Sz=((cellsize+(lpy-tol)-fact)**2)/(four*cellsize*(lpy-tol))
          dSz=-((cellsize+(lpy-tol)-fact)/(two*(lpy-tol)*cellsize)) 
        END IF   
     ELSE IF(g_coord(2,nod)<=gm_coord(2,s))THEN
        fact=ydist !--ydist negative (particle over the node)
        IF(fact>=(-lpy))THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lpy**2)/(two*cellsize*lpy))
         dSz=-(fact/(cellsize*lpy))
        ELSE IF(fact<=(-lpy).and.fact>=(-(cellsize-lpy)))THEN
         Sz=one+(fact/cellsize)
         dSz=one/cellsize
        ELSE IF(fact<=(-(cellsize-lpy)).and.fact>(-(cellsize+lpy-tol)))THEN
         Sz=((cellsize+(lpy-tol)+fact)**2)/(four*cellsize*(lpy-tol))
         dSz=((cellsize+(lpy-tol)+fact)/(two*(lpy-tol)*cellsize))
        END IF  
     END IF  

     funGIMP(m,1)=Sr*Sz
     funGIMP(m+1,2)=Sz*Sr
     IF(Sr<=0.0 .or. Sz<=0.0)THEN
      derGIMP(1,k)=0.0
      derGIMP(2,k)=0.0
     ELSE  
     derGIMP(1,k)=dSr*SZ
     derGIMP(2,k)=-dSz*Sr
     END IF
     k=k+1
     m=m+2
    END IF
    i=i+1
   END DO 
   RETURN
  END SUBROUTINE GIMP_funder2

  SUBROUTINE iGIMP_funder(s,mpoints,elemmpoints,nip,coord,cellsize,gm_coord,gimptol,derGIMP,funGIMP)
     !     
     ! Subroutine to compute implicit GIMP shape functions and shape functions derivatives 
     ! to send information to the nodes (this subroutine is diferent from GIMP_funder only
     ! because the values of the functions and derivatives are stored in a different way)
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::s,nip,elemmpoints(:,:)
   REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:),cellsize,mpoints(:,:),gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::twolp,two=2.0_iwp,minx,miny,Sr,Sz,dSr,dSz,lp,factx,facty,maxx,maxy
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,tol
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   REAL(iwp)::xi1,xi2,ni1,ni2,p5=-0.5_iwp
   INTEGER::i,j,k,m,nod,side
   !LOGICAL,INTENT(IN)::gimptol
   
   derGIMP=0.0
   funGIMP=0.0
   Sr=0.0
   Sz=0.0
   dSr=0.0
   dSz=0.0
   IF(nip==1)lp=two/(nip)
   IF(nip==4)lp=two/(nip/two)
   i=1
   k=1
   m=1
   !tol=1.0e-8

   !IF(gimptol)THEN
   !  tol=1.0e-8  !-tolerance values is to avoid the effect of nodes over the support domain boundarie of particles, that can give errors
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol
   
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=UBOUND(coord,1)
     minx=coord(1,1)       !--x coordinate of node 'i'
     maxx=coord(3,1)
     miny=coord(1,2)       !--y coordinate of node 'i'
     maxy=coord(2,2)
     xpoint=gm_coord(1,s)    !--x coordinate of current material point 's'
     ypoint=gm_coord(2,s)    !--y coordinate of material point 's'

     !--Evaluate the local coordinate of the material point inside and outside the element affected
     IF(xpoint>=minx.and.xpoint<=maxx)THEN !--Material point inside the element
      factx=mpoints(s,1)
     ELSE     !--Material pooint outside the element
      IF(xpoint<minx)THEN !--point outside the element in the x direction
        xdist=coord(1,1)-xpoint
        factx=((cellsize/two+xdist)/(cellsize/two))*(-one)   !--local coordinate point at the left side of the element
      ELSE IF(xpoint>maxx)THEN
        xdist=xpoint-coord(3,1)
        factx=((cellsize/two+xdist)/(cellsize/two))          !--local coordinate point at the right side of the element
      END IF
     END IF

     IF(ypoint>=miny.and.ypoint<=maxy)THEN
      facty=mpoints(s,2)
     ELSE 
      IF(ypoint>maxy)THEN !--point outside the element in the y direction
        ydist=(coord(2,2)-ypoint)*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))          !--local coordinate point over the element
      ELSE IF(ypoint<miny)THEN
        ydist=(ypoint-coord(1,2))*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))*(-one)   !--local coordinate point below the element
      END IF  
     END IF
     !--End of the evaluation of the local coordinate of the material point
     

     !--Evaluate integration limits Xi1 and Xi2 in x direction
     IF(factx-(lp/2.0)<=(-one))THEN
          xi1=-one 
       ELSE
          xi1=factx-(lp/2.0)
     END IF      
 
     IF(factx+(lp/2.0)>=(one))THEN
          xi2=one 
       ELSE
          xi2=factx+(lp/2.0)
     END IF 

     !--Evaluate integration limits Xi1 and Xi2 in x direction

     IF(facty-(lp/2.0)<=(-one))THEN
          ni1=-one 
       ELSE
          ni1=facty-(lp/2.0)
     END IF
     
     IF(facty+(lp/2.0)>=(one))THEN
          ni2=one 
       ELSE
          ni2=facty+(lp/2.0)
     END IF 
   
    !--End of local coordinates and integration limits


    DO i=1,nod  !Shape functions and shape functions derivatives for an element
        
      IF(i==1)THEN
       Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
       IF(factx>=p5)THEN
         dSr=(Xi1-Xi2)/(two*lp)
       ELSE  
         dSr=2*p5-factx
       END IF
       Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
       IF(facty>=p5)THEN
         dSz=(ni1-ni2)/(two*lp)
       ELSE
         dSz=2*p5-facty
       END IF  

      ELSE IF(i==2)THEN
       Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
         IF(factx>=p5)THEN
           dSr=(Xi1-Xi2)/(two*lp)
         ELSE  
           dSr=2*p5-factx
         END IF
       Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
         IF(facty<=-p5)THEN
           dSz=(ni2-ni1)/(two*lp)
         ELSE
           dSz=2*(-p5)-facty
         END IF  

      ELSE IF(i==3)THEN
       Sr=(one/(four*lp))*(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
         IF(factx<=-p5)THEN
          dSr=(Xi2-Xi1)/(two*lp)
         ELSE
          dSr=-factx+2*(-p5)
         END IF
       Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
         IF(facty<=-p5)THEN
           dSz=(ni2-ni1)/(two*lp)
         ELSE
           dSz=2*(-p5)-facty
         END IF 

      ELSE IF(i==4)THEN 
       Sr=(one/(four*lp))*(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
         IF(factx<=-p5)THEN
          dSr=(Xi2-Xi1)/(two*lp)
         ELSE
          dSr=-factx+2*(-p5)
         END IF
       Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
         IF(facty>=p5)THEN
           dSz=(ni1-ni2)/(two*lp)
         ELSE
           dSz=2*p5-facty
         END IF
         
      END IF  
      
     funGIMP(m)=Sr*Sz
     IF(Sr<=0.0.or.Sz<=0.0)THEN
      derGIMP(1,k)=0.0
      derGIMP(2,k)=0.0
     ELSE  
     derGIMP(1,k)=dSr*Sz
     derGIMP(2,k)=dSz*Sr
     END IF
     k=k+1
     m=m+1

    END DO
    

   RETURN
  END SUBROUTINE iGIMP_funder

  SUBROUTINE Sup_coord(s,GIMP_nodes,g_coord,jac_coord)
     !     
     ! Subroutine to save the coordenates of the activated nodes in jac_coord
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::GIMP_nodes(:,:),s
   REAL(iwp),INTENT(IN)::g_coord(:,:)
   REAL(iwp),INTENT(OUT)::jac_coord(:,:)
   INTEGER::i,j

   j=1
   DO i=1,9
     IF(GIMP_nodes(i,s)>0)THEN
       jac_coord(j,:)=g_coord(:,GIMP_nodes(i,s))
       j=j+1
     END IF  
   END DO  
   
  RETURN
  END SUBROUTINE Sup_coord
 
   SUBROUTINE eldformgimp(s,eld,loads,nf,GIMP_nodes,values,g_g)
    !
    ! Subroutine to create the eld steering vector
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     INTEGER,INTENT(IN)::GIMP_nodes(:,:),values,g_g(:,:),nf(:,:),s
     REAL(iwp),INTENT(IN)::loads(:)
     REAL(iwp),INTENT(OUT)::eld(:)
     INTEGER,ALLOCATABLE::b(:)
     INTEGER::i,n
     
     ALLOCATE(b(2))
     eld=0.0_iwp
     n=1
        !loads(b+1) is due because when sending the load vector to the subroutine, 
        !there is no loads(0) value anymore, every value is now loads(0+1)
       DO i=1,values
           b=nf(:,GIMP_nodes(i,s))
           eld(n:n+1)=loads(b+1)
           n=n+2
        END DO
  
   RETURN 
  END SUBROUTINE eldformgimp

     SUBROUTINE eldformgimp2(s,eld,nf,GIMP_nodes,values,g_g)
    !
    ! Subroutine to create the eld steering vector
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     INTEGER,INTENT(IN)::GIMP_nodes(:,:),values,g_g(:,:),nf(:,:),s
     INTEGER,INTENT(OUT)::eld(:)
     INTEGER::i,n,m
     
     eld=0.0_iwp
     n=1
     m=1
        !loads(b+1) is due because when sending the load vector to the subroutine, 
        !there is no loads(0) value anymore, every value is now loads(0+1)
       DO i=1,values
           eld(m)=nf(1,GIMP_nodes(i,s))
           eld(m+1)=nf(2,GIMP_nodes(i,s))
           m=m+2
          
        END DO
  
   RETURN 
  END SUBROUTINE eldformgimp2

     SUBROUTINE gimpfunform(s,iel,eldddylds,nf,GIMP_nodes,values,g_g,mvval)
    !
    ! Subroutine to create the shape function vector to interpolate values from particles to nodes
    !
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     INTEGER,INTENT(IN)::iel,GIMP_nodes(:,:),g_g(:,:),nf(:,:),s,values,mvval
     INTEGER,INTENT(OUT)::eldddylds(:)
     INTEGER::i,n
     
     eldddylds=0
     n=1
        !loads(b+1) is due because when sending the load vector to the subroutine, 
        !there is no loads(0) value anymore, every value is now loads(0+1)
       DO i=1,values
           !eldddylds(n:n+1)=loads(b+1)
           !IF(num(s)==8.or.num(s)==23.or.num(s)==38.or.num(s)==53.or.num(s)==68.or.num(s)==83.or.num(s)==98)THEN
           !If((GIMP_nodes(i,s)==8.or.GIMP_nodes(i,s)==23.or.GIMP_nodes(i,s)==38.or.GIMP_nodes(i,s)==53.or. &
           !GIMP_nodes(i,s)==68.or.GIMP_nodes(i,s)==83.or.GIMP_nodes(i,s)==98).and.mvval==2)THEN
           !    eldddylds(n:n+1)=0
           !ELSE   
           eldddylds(n:n+1)=nf(:,GIMP_nodes(i,s))
           !END IF    
           n=n+2
        END DO
  
   RETURN 
  END SUBROUTINE gimpfunform

      SUBROUTINE beematgimp(s,bee,deriv,values)
     !
     ! This subroutine forms the bee matrix for 16 and 9 nodes.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::values,s
      REAL(iwp),INTENT(IN)::deriv(:,:)
      REAL(iwp),INTENT(OUT)::bee(:,:)
      INTEGER::l,m,n,i,nod,k
      REAL(iwp)::x,y

      bee=0.0_iwp
      
      DO i=1,values
       k=2*i
       l=k-1
       x=deriv(1,i)
       y=deriv(2,i)
       bee(1,l)=x
       bee(3,k)=x
       bee(2,k)=y
       bee(3,l)=y
      END DO
     
     RETURN
     END SUBROUTINE beematgimp

   SUBROUTINE elemposition(elemmpoints,k,iel,nx1,position)
     !
     ! This subroutine save the position of the material point in the elemets suport domain
      !
     !     -------------------------
     !     -           -           -
     !     -     1     -     2     -
     !     -           -           -
     !     -           -      x    -        position = 2
     !     -------------------------
     !     -           -           -
     !     -           -           -
     !     -     3     -     4     -
     !     -           -           -
     !     -------------------------
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::k,elemmpoints(:,:),iel,nx1
      INTEGER,INTENT(OUT)::position

      IF(iel==elemmpoints(k,1))position=1
      IF(iel==elemmpoints(k,1)+1)position=2
      IF(iel==elemmpoints(k,1)+nx1)position=3
      IF(iel==elemmpoints(k,1)+nx1+1)position=4


   END SUBROUTINE elemposition

    SUBROUTINE Elem_suport(c_ele,nf,g_num,cellsize,nx1,nip,nels,g_coord,gm_coord,lp_mp,gimptol,smethod,elemmpoints)
     !
     ! This subroutine save all the elements inside the support domain of each material point
     !
     !     -------------------------
     !     -           -           -
     !     -     1     -     2     -
     !     -           -           -
     !     -           -           -
     !     -------------------------
     !     -           -           -
     !     -           -           -
     !     -     3     -     4     -
     !     -           -           -
     !     -------------------------
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::nels,g_num(:,:),nip,nf(:,:),c_ele(:),nx1,smethod
      REAL(iwp),INTENT(IN)::cellsize,gimptol
      REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),lp_mp(:,:)
      INTEGER,INTENT(OUT)::elemmpoints(:,:)
      INTEGER::l,m,n,i,nod,k,nmps,cont,colx,rowy,ielloc
      REAL(iwp)::x,y,two=2.0_iwp,twolp,lp,tol,neigbour,lpx,lpy
      !LOGICAL,INTENT(IN)::gimptol

      elemmpoints=0
      nmps=UBOUND(gm_coord,2)
      !twolp=cellsize/(nip/two)
      !lp=twolp/two
      
      !tol=1.0e-8

      
   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
      
      tol=gimptol

   !The neigbour option will allow neigbour elements to be affected by points outside
   !If(smethod==3)neigbour=1
   !If(smethod>3)neigbour=2 !IF neigbour=1 then the element affected is only the one containing the material point
              !IF neigbour=2 then the element can be affected by points in neigbour elements
   neigbour=2
   IF(neigbour==1)THEN
     cont=1
     DO k=1,nmps
      colx=gm_coord(1,k)/cellsize+1.0
      rowy=ABS(gm_coord(2,k))/cellsize+2.0
      ielloc=(rowy-1.0)*nx1+colx
      elemmpoints(k,cont)=ielloc
     END DO   
   ELSE
     DO k=1,nmps
       lpx=lp_mp(1,k);lpy=lp_mp(2,k)
       cont=1
       Ele:DO i=1,nels
         ActiveEl:IF(((nf(1,g_num(1,i))>0.or.nf(2,g_num(1,i))>0).and.(nf(1,g_num(2,i))>0.or.nf(2,g_num(2,i))>0).and.  &
            (nf(1,g_num(3,i))>0.or.nf(2,g_num(3,i))>0).and.(nf(1,g_num(4,i))>0.or.nf(2,g_num(4,i))>0)).or.   &
            c_ele(i)>0)THEN ! Check if the element is active (4 nodes of the element free at least in one direction or with a material point inside)
            m=g_num(3,i)
            n=g_num(2,i)
         
           IF((gm_coord(1,k)<=g_coord(1,m)+lpx-tol).and.(gm_coord(1,k)>=g_coord(1,n)-lpx+tol))THEN !Chek if the material point is inside the x suport doman of one element
             m=g_num(2,i);n=g_num(1,i)
           IF((gm_coord(2,k)<=g_coord(2,m)+lpy-tol).and.(gm_coord(2,k)>=g_coord(2,n)-lpy+tol))THEN !If is inside the x support domain, then chek if the material point is inside the y suport doman of one element
              elemmpoints(k,cont)=i
              cont=cont+1
              IF(cont==5)GOTO 20
           END IF
           END IF
        
         END IF ActiveEl
       END DO Ele
 20  CONTINUE  
      END DO
     END IF  
     
     RETURN
    END SUBROUTINE Elem_suport

   SUBROUTINE Point_support(g_num,cellsize,nip,nels,g_coord,gm_coord,gimptol,pointsupport,nodeweight)
     !
     ! This subroutine save all the elements inside the support domain of each material point
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::nels,g_num(:,:),nip
      REAL(iwp),INTENT(IN)::cellsize
      REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),gimptol
      INTEGER,INTENT(OUT)::pointsupport(:,:),nodeweight(:,:)
      REAL(iwp),ALLOCATABLE::supdom(:,:)
      INTEGER::s,nod,k,j,nmps,cont,colx,rowy,ielloc,der,izq,arr,aba
      REAL(iwp)::x,y,two=2.0_iwp,twolp,lp,tol,minx,maxx,miny,maxy
      !LOGICAL,INTENT(IN)::gimptol

      pointsupport=0
      nodeweight=0
      nmps=UBOUND(gm_coord,2)
      twolp=cellsize/(nip/two)
      lp=twolp/two
      !tol=1.0e-8 
      ALLOCATE(supdom(4,2))

      
   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
      
      tol=gimptol

  DO j=1,nels 
    der=1;izq=1;arr=1;aba=1
    cont=1
   DO s=1,nmps
     !-Evaluation of the coordinates of the boundaries of the support domain for the GIMP
     supdom(1,1)=gm_coord(1,s)-(cellsize+lp) !Left boundarie of the support domain
     supdom(1,2)=gm_coord(2,s)+(cellsize+lp) !Upper boundarie of the support domain
     supdom(2,1)=gm_coord(1,s)+(cellsize+lp) !Right boundarie of the support domain
     supdom(4,2)=gm_coord(2,s)-(cellsize+lp) !Lower boundarie of the support domain
     !-Tolerance is included to be sure the nodes are considered inside a support domain
     !-when they are just over the the edge of the boudarie
     minx=supdom(1,1)+tol
     maxx=supdom(2,1)-tol
     miny=supdom(4,2)+tol
     maxy=supdom(1,2)-tol
     
     
    
    IF(g_coord(1,g_num(1,j))>=minx.and.g_coord(1,g_num(4,j))<=maxx)THEN ! x coordenate of the node is inside the domain 
     IF(g_coord(2,g_num(2,j))<=maxy.and.g_coord(2,g_num(1,j))>=miny)THEN  ! y coordenate of the node is inside the domain

      IF(g_coord(1,g_num(1,j))>gm_coord(1,s).and.g_coord(2,g_num(2,j))<gm_coord(2,s))THEN
        nodeweight(j,4)=nodeweight(j,4)+1

      ELSE IF(g_coord(1,g_num(1,j))>gm_coord(1,s).and.g_coord(2,g_num(2,j))>gm_coord(2,s)     &
                       .and.g_coord(2,g_num(1,j))<gm_coord(2,s))THEN
        nodeweight(j,3)=nodeweight(j,3)+1
        nodeweight(j,4)=nodeweight(j,4)+1
      ELSE IF(g_coord(1,g_num(1,j))>gm_coord(1,s).and.g_coord(2,g_num(1,j))>gm_coord(2,s))THEN  
        nodeweight(j,3)=nodeweight(j,3)+1
      END IF 


      IF(g_coord(1,g_num(1,j))<gm_coord(1,s).and.g_coord(1,g_num(3,j))>gm_coord(1,s).and.      &
                          g_coord(2,g_num(2,j))<gm_coord(2,s))THEN
        nodeweight(j,1)=nodeweight(j,1)+1
        nodeweight(j,4)=nodeweight(j,4)+1
      ELSE IF(g_coord(1,g_num(1,j))<gm_coord(1,s).and.g_coord(1,g_num(3,j))>gm_coord(1,s).and.       &
            g_coord(2,g_num(2,j))>gm_coord(2,s).and.g_coord(2,g_num(1,j))<gm_coord(2,s))THEN
        nodeweight(j,1)=nodeweight(j,1)+1
        nodeweight(j,2)=nodeweight(j,2)+1
        nodeweight(j,3)=nodeweight(j,3)+1
        nodeweight(j,4)=nodeweight(j,4)+1
      ELSE IF(g_coord(1,g_num(1,j))<gm_coord(1,s).and.g_coord(1,g_num(3,j))>gm_coord(1,s).and.       &
                        g_coord(2,g_num(1,j))>gm_coord(2,s))THEN
        nodeweight(j,2)=nodeweight(j,2)+1
        nodeweight(j,3)=nodeweight(j,3)+1
      END IF 

      IF(g_coord(1,g_num(3,j))<gm_coord(1,s).and.g_coord(2,g_num(2,j))<gm_coord(2,s))THEN
        nodeweight(j,1)=nodeweight(j,1)+1
      ELSE IF(g_coord(1,g_num(3,j))<gm_coord(1,s).and.g_coord(2,g_num(2,j))>gm_coord(2,s).and.  &
             g_coord(2,g_num(1,j))<gm_coord(2,s))THEN
        nodeweight(j,1)=nodeweight(j,1)+1
        nodeweight(j,2)=nodeweight(j,2)+1
      ELSE IF(g_coord(1,g_num(3,j))<gm_coord(1,s).and.g_coord(2,g_num(1,j))>gm_coord(2,s))THEN  
        nodeweight(j,2)=nodeweight(j,2)+1
      END IF 

         pointsupport(j,cont)=s
         cont=cont+1 
         IF(cont==10)GOTO 20
        END IF
      END IF  
    END DO  
20  CONTINUE    
   END DO 
     
   RETURN
  END SUBROUTINE Point_support

    SUBROUTINE weightvalue(a_ele,k,iel,mpoints,coord,cellsize,nip,gm_coord,gimptol,waverage)           
     !
     ! This subroutine computes the area of the particle support domain inside the element
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::nip,iel,k,a_ele(:)
      REAL(iwp),INTENT(IN)::cellsize,mpoints(:,:),gimptol
      REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:)
      REAL(iwp),INTENT(OUT)::waverage
      INTEGER::l,m,n,i,nod,nmps,cont
      REAL(iwp)::x,y,two=2.0_iwp,twolp,lp,centerx,centery,xdif,ydif
      !LOGICAL,INTENT(IN)::gimptol

      twolp=cellsize/(nip/two)
      lp=twolp/two
      waverage=0.0

      !-Center of the element in the domain
      centerx=((coord(3,1)-coord(2,1))/2.0)+coord(2,1)
      centery=((coord(4,2)-coord(3,2))/2.0)+coord(3,2)

      IF(gm_coord(1,k)<=centerx)THEN !-The material point is in the left side of the element
        
         IF(gm_coord(2,k)>=centery)THEN !-The material point is over the center of the element
           IF(mpoints(k,1)>=-0.5_iwp.and.mpoints(k,2)<=0.5_iwp.and.a_ele(k)==iel)THEN
             waverage=twolp**2  !- The whole support domain is inside the element
           ELSE
             xdif=gm_coord(1,k)+lp-(centerx-twolp);IF(xdif>twolp)xdif=twolp
             ydif=(gm_coord(2,k)-lp-(centery+twolp))*(-1.0_iwp);IF(ydif>twolp)ydif=twolp
             waverage=xdif*ydif
           END IF    

         ELSE !-The mp is below the center of the element
           IF(mpoints(k,1)>=-0.5_iwp.and.mpoints(k,2)>=-0.5_iwp.and.a_ele(k)==iel)THEN
             waverage=twolp**2  !- The whole support domain is inside the element
           ELSE
             xdif=gm_coord(1,k)+lp-(centerx-twolp);IF(xdif>twolp)xdif=twolp
             ydif=((centery-twolp)-(gm_coord(2,k)+lp))*(-1.0_iwp);IF(ydif>twolp)ydif=twolp
             waverage=xdif*ydif
           END IF 
         END IF

      ELSE !-The material point is at the right side of the element

        IF(gm_coord(2,k)>=centery)THEN !-The material point is over the center of the element
           IF(mpoints(k,1)<=0.5_iwp.and.mpoints(k,2)<=0.5_iwp.and.a_ele(k)==iel)THEN
             waverage=twolp**2  !- The whole support domain is inside the element
           ELSE
             xdif=(centerx+twolp)-(gm_coord(1,k)-lp);IF(xdif>twolp)xdif=twolp
             ydif=(gm_coord(2,k)-lp-(centery+twolp))*(-1.0_iwp);IF(ydif>twolp)ydif=twolp
             waverage=xdif*ydif
           END IF    

         ELSE !-The mp is below the center of the element
           IF(mpoints(k,1)<=0.5_iwp.and.mpoints(k,2)>=-0.5_iwp.and.a_ele(k)==iel)THEN
             waverage=twolp**2  !- The whole support domain is inside the element
           ELSE
             xdif=(centerx+twolp)-(gm_coord(1,k)-lp);IF(xdif>twolp)xdif=twolp
             ydif=((centery-twolp)-(gm_coord(2,k)+lp))*(-1.0_iwp);IF(ydif>twolp)ydif=twolp
             waverage=xdif*ydif
           END IF 
         END IF

      END IF  

    END SUBROUTINE weightvalue

    SUBROUTINE km_weight(mpoints,k,km_w)
     !
     ! This subroutine computes the weightning for each node before add each 
     ! stiffnes matrix to the global stiffnes matrix
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::k
      REAL(iwp),INTENT(IN)::mpoints(:,:)
      REAL(iwp),INTENT(OUT)::km_w(:)
      REAL(iwp)::zero=0.0_iwp,mp5=-0.5_iwp,pp5=0.5_iwp,one=1.0_iwp,p25=0.25_iwp

      km_w=zero

      IF(mpoints(k,1)<mp5.and.mpoints(k,2)<mp5)THEN !Section 1
         km_w(1)=one;km_w(2)=pp5;km_w(3)=one;km_w(4)=pp5;km_w(5)=p25
         km_w(6)=pp5;km_w(7)=one;km_w(8)=pp5;km_w(9)=one
      ELSE IF(mpoints(k,1)<mp5.and.mpoints(k,2)>=mp5.and.mpoints(k,2)<=pp5)THEN !Section 2
         km_w(1)=one;km_w(2)=one;km_w(3)=pp5;km_w(4)=pp5;km_w(5)=one
         km_w(6)=one
      ELSE IF(mpoints(k,1)<mp5.and.mpoints(k,2)>pp5)THEN  !Section 3
         km_w(1)=one;km_w(2)=pp5;km_w(3)=one;km_w(4)=pp5;km_w(5)=p25
         km_w(6)=pp5;km_w(7)=one;km_w(8)=pp5;km_w(9)=one
      ELSE IF(mpoints(k,1)>=mp5.and.mpoints(k,1)<=pp5.and.mpoints(k,2)>pp5)THEN   !Section 4
         km_w(1)=one;km_w(2)=pp5;km_w(3)=one;km_w(4)=one;km_w(5)=pp5
         km_w(6)=one
      ELSE IF(mpoints(k,1)>=pp5.and.mpoints(k,2)>pp5)THEN   !Section 5
         km_w(1)=one;km_w(2)=pp5;km_w(3)=one;km_w(4)=pp5;km_w(5)=p25
         km_w(6)=pp5;km_w(7)=one;km_w(8)=pp5;km_w(9)=one
      ELSE IF(mpoints(k,1)>=pp5.and.mpoints(k,2)<=pp5.and.mpoints(k,2)>=mp5)THEN   !Section 6
         km_w(1)=one;km_w(2)=one;km_w(3)=pp5;km_w(4)=pp5;km_w(5)=one
         km_w(6)=one
      ELSE IF(mpoints(k,1)>=pp5.and.mpoints(k,2)<mp5)THEN  !Section 7
         km_w(1)=one;km_w(2)=pp5;km_w(3)=one;km_w(4)=pp5;km_w(5)=p25
         km_w(6)=pp5;km_w(7)=one;km_w(8)=pp5;km_w(9)=one
      ELSE IF(mpoints(k,1)>=mp5.and.mpoints(k,1)<=pp5.and.mpoints(k,2)<mp5)THEN !Section 8
         km_w(1)=one;km_w(2)=pp5;km_w(3)=one;km_w(4)=one;km_w(5)=pp5
         km_w(6)=one
      ELSE    !Section 9
         km_w(1)=one;km_w(2)=one;km_w(3)=one;km_w(4)=one
      END IF     

    END SUBROUTINE km_weight

  SUBROUTINE iGIMP_funder2(s,mpoints,nip,coord,cellsize,gm_coord,gimptol,GIMP_nodes,  &
                             g_coord,a_ele,c_ele,iel,derGIMP,funGIMP)
     !     
     ! Subroutine to compute implicit GIMP shape functions and shape functions derivatives 
     ! to send information to the nodes (this subroutine is diferent from GIMP_funder only
     ! because the values of the functions and derivatives are stored in a different way)
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::s,nip,GIMP_nodes(:,:),a_ele(:),iel,c_ele(:)
   REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:),cellsize,mpoints(:,:),g_coord(:,:),gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::twolp,two=2.0_iwp,minx,miny,Sr,Sz,dSr,dSz,lp,factx,facty,maxx,maxy
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,tol
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   REAL(iwp)::xi1,xi2,ni1,ni2,p5=-0.5_iwp,coordx,coordy
   INTEGER::i,j,k,m,n,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol
   
   derGIMP=0.0
   funGIMP=0.0
   IF(nip==1)lp=two/(nip)
   IF(nip==4)lp=two/(nip/two)
   i=1
   k=1
   m=1
   !tol=1.0e-8
   !IF(gimptol)THEN
   !  tol=1.0e-8  !-tolerance values is to avoid the effect of nodes over the support domain boundarie of particles, that can give errors
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol
   
     xpoint=gm_coord(1,s)    !--x coordinate of current material point 's'
     ypoint=gm_coord(2,s)    !--y coordinate of material point 's'
     !nod=UBOUND(GIMP_nodes,1)
     nodes=UBOUND(derGIMP,2)

    Fun_Deriv:DO n=1,nodes
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(n,s)
     coordx=g_coord(1,nod)
     coordy=g_coord(2,nod)
     !coord=g_coord(:,GIMP_nodes(j,s))
     minx=coordx               !--x coordinate of node 'i'
     maxx=coordx+cellsize
     miny=coordy-cellsize      
     maxy=coordy               !--y coordinate of node 'i'

     !--Evaluate the local coordinate of the material point inside and outside the element affected
     IF(ABS(xpoint-coordx)<=cellsize)THEN !--Material point inside the element
      factx=mpoints(s,1)
     ELSE     !--Material pooint outside the element
      IF(xpoint<coordx)THEN !--point outside the element in the x direction
        xdist=(coordx-cellsize)-xpoint
        factx=((cellsize/two+xdist)/(cellsize/two))*(-one)   !--local coordinate point at the left side of the element
      ELSE IF(xpoint>coordx)THEN
        xdist=xpoint-(coordx+cellsize)
        factx=((cellsize/two+xdist)/(cellsize/two))          !--local coordinate point at the right side of the element
      END IF
     END IF

     IF(ABS(ypoint-coordy)<=cellsize)THEN
      facty=mpoints(s,2)
     ELSE 
      IF(ypoint>coordy)THEN !--point outside the element in the y direction
        ydist=((coordy+cellsize)-ypoint)*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))          !--local coordinate point over the element
      ELSE IF(ypoint<coordy)THEN
        ydist=(ypoint-(coordy-cellsize))*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))*(-one)   !--local coordinate point below the element
      END IF  
     END IF
     !--End of the evaluation of the local coordinate of the material point
     

     !--Evaluate integration limits Xi1 and Xi2 in x direction
     IF(factx-(lp/2.0)<=(-one).and.ABS(factx)<=1.0)THEN
          xi1=-one 
       ELSE IF(factx-(lp/2.0)>(-one))THEN
          xi1=factx-(lp/2.0)
       ELSE IF(factx-(lp/2.0)<=(-one).and.iel/=a_ele(s).and.ABS(factx)>=1.0)THEN   !Material point in another element and computing derivatives for empty elements
          xi1=factx-(lp/2.0)
     END IF      
 
     IF(factx+(lp/2.0)>=(one).and.ABS(factx)<=1.0)THEN
          xi2=one 
       ELSE IF(factx+(lp/2.0)<(one))THEN
          xi2=factx+(lp/2.0)
       ELSE IF(factx+(lp/2.0)>=(one).and.iel/=a_ele(s).and.ABS(factx)>=1.0)THEN   
          xi2=factx+(lp/2.0)
     END IF 

     !--Evaluate integration limits Xi1 and Xi2 in x direction

     IF(facty-(lp/2.0)<=(-one).and.ABS(facty)<=1.0)THEN
          ni1=-one 
       ELSE IF(facty-(lp/2.0)>(-one))THEN
          ni1=facty-(lp/2.0)
       ELSE IF(facty-(lp/2.0)<=(-one).and.iel/=a_ele(s).and.ABS(facty)>=1.0)THEN
          ni1=facty-(lp/2.0)
     END IF
     
     IF(facty+(lp/2.0)>=(one).and.ABS(facty)<=1.0)THEN
          ni2=one 
       ELSE IF(facty+(lp/2.0)<(one))THEN
          ni2=facty+(lp/2.0)
       ELSE IF(facty+(lp/2.0)>=(one).and.iel/=a_ele(s).and.ABS(facty)>=1.0)THEN
         ni2=facty+(lp/2.0)
     END IF 
   
    !--End of local coordinates and integration limits


    !Deriv_Fun:DO i=1,nod  !Shape functions and shape functions derivatives for an element
        
      IF(xpoint>=coordx.and.ypoint>=coordy)THEN  !Comparing with node 1
       Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
       IF(ABS(factx)>=1.0)Sr=(one/(four*lp))*(two*one-one**2-two*Xi1+Xi1**2) !To recover the original shape function without modify xi1 and xi2
       IF(factx>=p5)THEN
         dSr=(Xi1-Xi2)/(two*lp)
       ELSE  
         dSr=2*p5-factx
       END IF
       Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
       IF(ABS(facty)>=1.0)Sz=(one/(four*lp))*(two*one-one**2-two*ni1+ni1**2)
       IF(facty>=p5)THEN
         dSz=(ni1-ni2)/(two*lp)
       ELSE
         dSz=2*p5-facty
       END IF  

       

      ELSE IF(xpoint>=coordx.and.ypoint<=coordy)THEN
       Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
       IF(ABS(factx)>=1.0)Sr=(one/(four*lp))*(two*one-one**2-two*Xi1+Xi1**2) !To recover the original shape function without modify xi1 and xi2
         IF(factx>=p5)THEN
           dSr=(Xi1-Xi2)/(two*lp)
         ELSE  
           dSr=2*p5-factx
         END IF
       Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
       IF(ABS(facty)>=1.0)Sz=(one/(four*lp))*(two*ni2+ni2**2-two*(-one)-(-one)**2)
         IF(facty<=-p5)THEN
           dSz=(ni2-ni1)/(two*lp)
         ELSE
           dSz=2*(-p5)-facty
         END IF  

      ELSE IF(xpoint<=coordx.and.ypoint<=coordy)THEN
         Sr=(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
         Sr=Sr*(one/(four*lp))
         IF(ABS(factx)>=1.0)Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*(-one)+(-one)**2) !To recover the original shape function without modify xi1 and xi2
         IF(factx<=-p5)THEN
          dSr=(Xi2-Xi1)/(two*lp)
         ELSE
          dSr=-factx+2*(-p5)
         END IF
       Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
       IF(ABS(facty)>=1.0)Sz=(one/(four*lp))*(two*ni2+ni2**2-two*(-one)-(-one)**2)
         IF(facty<=-p5)THEN
           dSz=(ni2-ni1)/(two*lp)
         ELSE
           dSz=2*(-p5)-facty
         END IF 

      ELSE IF(xpoint<=coordx.and.ypoint>=coordy)THEN 
       Sr=(one/(four*lp))*(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
       IF(ABS(factx)>=1.0)Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*(-one)+(-one)**2)  !To recover the original shape function without modify xi1 and xi2
         IF(factx<=-p5)THEN
          dSr=(Xi2-Xi1)/(two*lp)
         ELSE
          dSr=-factx+2*(-p5)
         END IF
       Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
       IF(ABS(facty)>=1.0)Sz=(one/(four*lp))*(two*one-one**2-two*ni1+ni1**2)
         IF(facty>=p5)THEN
           dSz=(ni1-ni2)/(two*lp)
         ELSE
           dSz=2*p5-facty
         END IF
         
      END IF  
      
     funGIMP(m)=Sr*Sz
     IF((Sr<=0.0.or.Sz<=0.0))THEN
      derGIMP(1,m)=0.0
      derGIMP(2,m)=0.0
     ELSE  
     derGIMP(1,m)=dSr*Sz
     derGIMP(2,m)=dSz*Sr
     IF((derGIMP(1,m)>0.0.or.derGIMP(2,m)>0.0).and.n==1)THEN
       derGIMP(1,m)=0.0
       derGIMP(2,m)=0.0
       funGIMP(m)=0.0
     ELSE IF((derGIMP(1,m)>0.0.or.derGIMP(2,m)<0.0).and.n==2)THEN
     derGIMP(1,m)=0.0
     derGIMP(2,m)=0.0  
     funGIMP(m)=0.0
     ELSE IF((derGIMP(1,m)<0.0.or.derGIMP(2,m)<0.0).and.n==3)THEN
     derGIMP(1,m)=0.0
     derGIMP(2,m)=0.0
     funGIMP(m)=0.0
     ELSE IF((derGIMP(1,m)<0.0.or.derGIMP(2,m)>0.0).and.n==4)THEN
     derGIMP(1,m)=0.0
     derGIMP(2,m)=0.0
     funGIMP(m)=0.0
     END IF
     
     END IF
    
     m=m+1

    !END DO Deriv_Fun

   END DO Fun_Deriv
    

   RETURN
  END SUBROUTINE iGIMP_funder2

  SUBROUTINE iGIMP_funder3(s,mpoints,lp_mp,nip,coord,cellsize,gm_coord,gimptol,GIMP_nodes,  &
                             g_coord,a_ele,c_ele,iel,derGIMP,funGIMP)
     !     
     ! Subroutine to compute implicit GIMP shape functions and shape functions derivatives. 
     ! iGIMP_funder3 uses a constant gradient at the center of the nodes and it reduces as the
     ! support domain ends
     
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::s,nip,GIMP_nodes(:,:),a_ele(:),iel,c_ele(:)
   REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:),cellsize,mpoints(:,:),g_coord(:,:),gimptol,lp_mp(:,:)
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::twolp,two=2.0_iwp,minx,miny,Sr,Sz,dSr,dSz,lp,factx,facty,maxx,maxy
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,tol
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   REAL(iwp)::xi1,xi2,ni1,ni2,mp5=-0.5_iwp,coordx,coordy,lpx,lpy
   INTEGER::i,j,k,m,n,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol
   
   
   derGIMP=0.0
   funGIMP=0.0
   !IF(nip==1)lp=two/(nip)
   !IF(nip==4)lp=two/(nip/two)
   !IF(nip==9)lp=two/3.0
   i=1
   k=1
   m=1
	!tol=1.0e-8
 !  IF(gimptol)THEN
 !    tol=1.0e-8  !-tolerance values is to avoid the effect of nodes over the support domain boundarie of particles, that can give errors
 !    CONTINUE
 !   ELSE 
 !     tol=0.0
 !  END IF
   tol=gimptol
   
     xpoint=gm_coord(1,s)    !--x coordinate of current material point 's'
     ypoint=gm_coord(2,s)    !--y coordinate of material point 's'
     !nod=UBOUND(GIMP_nodes,1)
     nodes=UBOUND(derGIMP,2)
     lpx=lp_mp(1,s)*(two/cellsize)
     lpy=lp_mp(2,s)*(two/cellsize)

    Fun_Deriv:DO n=1,nodes
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(n,1)
     coordx=g_coord(1,nod)
     coordy=g_coord(2,nod)
     !coord=g_coord(:,GIMP_nodes(j,s))

     !--Evaluate the local coordinate of the material point inside and outside the element affected
     IF(a_ele(s)==iel)THEN !--Material point inside the element
      factx=mpoints(s,1)
     ELSE     !--Material pooint outside the element
      IF(xpoint<=coordx)THEN !--point outside the element in the x direction (left side)
        IF(n==1.or.n==2)xdist=coordx-xpoint
        IF(n==3.or.n==4)xdist=(coordx-cellsize)-xpoint
        factx=((cellsize/two+xdist)/(cellsize/two))*(-one)   !--local coordinate point at the left side of the element
      ELSE IF(xpoint>=coordx)THEN  !--point outside the element in the x direction (right side)
        IF(n==1.or.n==2)xdist=xpoint-(coordx+cellsize)
          IF(n==3.or.n==4)xdist=xpoint-(coordx)
        factx=((cellsize/two+xdist)/(cellsize/two))          !--local coordinate point at the right side of the element
      END IF
     END IF

     IF(a_ele(s)==iel)THEN
      facty=mpoints(s,2)
     ELSE 
      IF(ypoint>=coordy)THEN !--point outside the element in the y direction (above the element)
        IF(n==1.or.n==4)ydist=((coordy+cellsize)-ypoint)*(-one)
        IF(n==2.or.n==3)ydist=((coordy)-ypoint)*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))          !--local coordinate point over the element
      ELSE IF(ypoint<=coordy)THEN  !--point outside the element in the y direction (below the element)
        IF(n==2.or.n==3)ydist=(ypoint-(coordy-cellsize))*(-one)
        IF(n==1.or.n==4)ydist=(ypoint-(coordy))*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))*(-one)   !--local coordinate point below the element
      END IF  
     END IF
     !--End of the evaluation of the local coordinate of the material point
     

     !Integration limits
     !--Evaluate integration limits Xi1
     !IF(factx-(lp/2.0)<=(-one))THEN
     !     xi1=-one 
     !  ELSE IF(factx-(lp/2.0)>(-one))THEN
     !     xi1=factx-(lp/2.0)                                                  
     !END IF                                                                       
     !
     !!--Evaluate integration limits Xi2
     !IF(factx+(lp/2.0)>=(one))THEN
     !     xi2=one 
     !  ELSE IF(factx+(lp/2.0)<(one))THEN
     !     xi2=factx+(lp/2.0)
     !END IF 
     !
     !!--Evaluate integration limits ni1
     !IF(facty-(lp/2.0)<=(-one))THEN
     !     ni1=-one 
     !  ELSE IF(facty-(lp/2.0)>(-one))THEN
     !     ni1=facty-(lp/2.0)
     !END IF
     !
     ! !--Evaluate integration limits ni1
     !IF(facty+(lp/2.0)>=(one))THEN
     !     ni2=one 
     !  ELSE IF(facty+(lp/2.0)<(one))THEN
     !     ni2=facty+(lp/2.0)
     !END IF 
   
    !--End of integration limits

    !Shape functions
    !IF(n==1)THEN  !Comparing with node 1
    !  Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
    !  Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
    !
    ! ELSE IF(n==2)THEN    !Comparing with node 2
    !  Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
    !  Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
    ! 
    ! ELSE IF(n==3)THEN   !Comparing with node 3
    !   Sr=(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
    !   Sr=Sr*(one/(four*lp))
    !   Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
    ! 
    ! ELSE IF(n==4)THEN      !Comparing with node 4
    !  Sr=(one/(four*lp))*(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
    !  Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
    !END IF  

    !Shape functions
    IF(n==1)THEN  !Comparing with node 1
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

     ELSE IF(n==2)THEN    !Comparing with node 2
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
     
     ELSE IF(n==4)THEN      !Comparing with node 4
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

     IF(n==1)THEN  !Comparing with node 1
        
          dSr=(Xi1-Xi2)/(two*lp)
          !IF(factx<(-one))dSr=-dSr

          dSz=(ni1-ni2)/(two*lp)
          !IF(facty<(-one))dSz=-dSz


      ELSE IF(n==2)THEN    !Comparing with node 2
       
           dSr=(Xi1-Xi2)/(two*lp)
           !IF(factx<(-one))dSr=-dSr

           dSz=(ni2-ni1)/(two*lp)
           !IF(facty>(one))dSz=-dSz
 
      ELSE IF(n==3)THEN   !Comparing with node 3
         
          dSr=(Xi2-Xi1)/(two*lp)
          !IF(factx>(one))dSr=-dSr

          dSz=(ni2-ni1)/(two*lp)
          !IF(facty>(one))dSz=-dSz
 
      ELSE IF(n==4)THEN      !Comparing with node 4
       
          dSr=(Xi2-Xi1)/(two*lp)
          !IF(factx>(one))dSr=-dSr

          dSz=(ni1-ni2)/(two*lp)
          !IF(facty<(-one))dSz=-dSz
         
      END IF  

      
      
      
     funGIMP(m)=Sr*Sz
     IF(abs(funGIMP(m))>2)THEN
        Sr=Sr
         PAUSE
     END IF    
     IF((Sr<=0.0.or.Sz<=0.0))THEN
      derGIMP(1,m)=0.0
      derGIMP(2,m)=0.0
     ELSE  
     derGIMP(1,m)=dSr*Sz
     derGIMP(2,m)=dSz*Sr
!$$$$$$      IF((derGIMP(1,m)>0.0.or.derGIMP(2,m)>0.0).and.n==1)THEN
!$$$$$$        derGIMP(1,m)=0.0
!$$$$$$        derGIMP(2,m)=0.0
!$$$$$$        funGIMP(m)=0.0
!$$$$$$      ELSE IF((derGIMP(1,m)>0.0.or.derGIMP(2,m)<0.0).and.n==2)THEN
!$$$$$$      derGIMP(1,m)=0.0
!$$$$$$      derGIMP(2,m)=0.0  
!$$$$$$      funGIMP(m)=0.0
!$$$$$$      ELSE IF((derGIMP(1,m)<0.0.or.derGIMP(2,m)<0.0).and.n==3)THEN
!$$$$$$      derGIMP(1,m)=0.0
!$$$$$$      derGIMP(2,m)=0.0
!$$$$$$      funGIMP(m)=0.0
!$$$$$$      ELSE IF((derGIMP(1,m)<0.0.or.derGIMP(2,m)>0.0).and.n==4)THEN
!$$$$$$      derGIMP(1,m)=0.0
!$$$$$$      derGIMP(2,m)=0.0
!$$$$$$      funGIMP(m)=0.0
!$$$$$$      END IF
     
     END IF
    
     m=m+1

    !END DO Deriv_Fun

   END DO Fun_Deriv

!$$$$$$    IF(a_ele(s)/=iel.and.ABS(factx)>one.and.ABS(facty)>one)THEN
!$$$$$$      derGIMP=derGIMP*10
!$$$$$$    END IF  
    

   RETURN
  END SUBROUTINE iGIMP_funder3

  SUBROUTINE iGIMP_funder4(s,mpoints,nip,coord,cellsize,gm_coord,gimptol,GIMP_nodes,  &
                             g_coord,a_ele,c_ele,iel,derGIMP,funGIMP)
     !     
     ! Subroutine to compute implicit GIMP shape functions and shape functions derivatives. 
     ! iGIMP_funder3 uses a constant gradient at the center of the nodes and it reduces as the
     ! support domain ends
     
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::s,nip,GIMP_nodes(:,:),a_ele(:),iel,c_ele(:)
   REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:),cellsize,mpoints(:,:),g_coord(:,:),gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),funGIMP(:)
   REAL(iwp)::two=2.0_iwp,minx,miny,Sr,Sz,dSr,dSz,lp,factx,facty,maxx,maxy
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,fact,tol,lpGIMP
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   REAL(iwp)::xi1,xi2,ni1,ni2,mp5=-0.5_iwp,coordx,coordy
   INTEGER::i,j,k,m,n,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol
   
   derGIMP=0.0
   funGIMP=0.0
   IF(nip==1)lp=two/(nip)
   IF(nip==4)lp=two/(nip/two)
   IF(nip==4)lpGIMP=cellsize/(nip/two);lpGIMP=lpGIMP/two
   i=1
   k=1
   m=1
   !tol=1.0e-8
   !
   !IF(gimptol)THEN
   !  tol=1.0e-8 !-tolerance values is to avoid the effect of nodes over the support domain boundarie of particles, that can give errors
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol
   
     xpoint=gm_coord(1,s)    !--x coordinate of current material point 's'
     ypoint=gm_coord(2,s)    !--y coordinate of material point 's'
     !nod=UBOUND(GIMP_nodes,1)
     nodes=UBOUND(derGIMP,2)

    Fun_Deriv:DO n=1,nodes
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=GIMP_nodes(n,s)
     coordx=g_coord(1,nod)
     coordy=g_coord(2,nod)
     !coord=g_coord(:,GIMP_nodes(j,s))
     minx=coordx               !--x coordinate of node 'i'
     maxx=coordx+cellsize
     miny=coordy-cellsize      
     maxy=coordy               !--y coordinate of node 'i'

     !--Evaluate the local coordinate of the material point inside and outside the element affected
     IF(a_ele(s)==iel)THEN !--Material point inside the element
      factx=mpoints(s,1)
     ELSE     !--Material pooint outside the element
      IF(xpoint<coordx)THEN !--point outside the element in the x direction (left side)
        IF(n==1.or.n==2)xdist=coordx-xpoint
        IF(n==3.or.n==4)xdist=(coordx-cellsize)-xpoint
        factx=((cellsize/two+xdist)/(cellsize/two))*(-one)   !--local coordinate point at the left side of the element
      ELSE IF(xpoint>coordx)THEN  !--point outside the element in the x direction (right side)
        IF(n==1.or.n==2)xdist=xpoint-(coordx+cellsize)
          IF(n==3.or.n==4)xdist=xpoint-(coordx)
        factx=((cellsize/two+xdist)/(cellsize/two))          !--local coordinate point at the right side of the element
      END IF
     END IF

     IF(a_ele(s)==iel)THEN
      facty=mpoints(s,2)
     ELSE 
      IF(ypoint>coordy)THEN !--point outside the element in the y direction (above the element)
        IF(n==1.or.n==4)ydist=((coordy+cellsize)-ypoint)*(-one)
        IF(n==2.or.n==3)ydist=((coordy)-ypoint)*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))          !--local coordinate point over the element
      ELSE IF(ypoint<coordy)THEN  !--point outside the element in the y direction (below the element)
        IF(n==2.or.n==3)ydist=(ypoint-(coordy-cellsize))*(-one)
        IF(n==1.or.n==4)ydist=(ypoint-(coordy))*(-one)
        facty=((cellsize/two+ydist)/(cellsize/two))*(-one)   !--local coordinate point below the element
      END IF  
     END IF
     !--End of the evaluation of the local coordinate of the material point
     

     !Integration limits
     !--Evaluate integration limits Xi1
     IF(factx-(lp/2.0)<=(-one))THEN
          xi1=-one 
       ELSE IF(factx-(lp/2.0)>(-one))THEN
          xi1=factx-(lp/2.0)                                                  
     END IF                                                                       
 
     !--Evaluate integration limits Xi2
     IF(factx+(lp/2.0)>=(one))THEN
          xi2=one 
       ELSE IF(factx+(lp/2.0)<(one))THEN
          xi2=factx+(lp/2.0)
     END IF 

     !--Evaluate integration limits ni1
     IF(facty-(lp/2.0)<=(-one))THEN
          ni1=-one 
       ELSE IF(facty-(lp/2.0)>(-one))THEN
          ni1=facty-(lp/2.0)
     END IF
     
      !--Evaluate integration limits ni1
     IF(facty+(lp/2.0)>=(one))THEN
          ni2=one 
       ELSE IF(facty+(lp/2.0)<(one))THEN
          ni2=facty+(lp/2.0)
     END IF 
   
    !--End of integration limits

    xdist=gm_coord(1,s)-g_coord(1,nod)
    ydist=(-gm_coord(2,s) )-(-g_coord(2,nod))
    !--GIMP shape functions and derivatives in the x direction
     IF(g_coord(1,nod)<=gm_coord(1,s))THEN
        fact=xdist !--particle at te right side
        IF(fact<=lpGIMP)THEN    ! fact is the local coordinate of the material point 
         Sr=one-((fact**2+lpGIMP**2)/(two*cellsize*lpGIMP))
        ELSE IF(fact>=lpGIMP.and.fact<=(cellsize-lpGIMP))THEN
         Sr=one-(fact/cellsize)
        ELSE IF(fact>=(cellsize-lpGIMP).and.fact<(cellsize+lpGIMP-tol))THEN
         Sr=((cellsize+lpGIMP-fact)**2)/(four*cellsize*lpGIMP)
        END IF     
     ELSE IF(g_coord(1,nod)>=gm_coord(1,s))THEN  
        fact=xdist !--particle at the left side
        IF(fact>=(-lpGIMP))THEN    
          Sr=one-((fact**2+lpGIMP**2)/(two*cellsize*lpGIMP))
        ELSE IF(fact<=(-lpGIMP).and.fact>=(-(cellsize-lpGIMP)))THEN
         Sr=one+(fact/cellsize)
        ELSE IF(fact<=(-(cellsize-lpGIMP)).and.fact>(-(cellsize+lpGIMP-tol)))THEN
         Sr=((cellsize+lpGIMP+fact)**2)/(four*cellsize*lpGIMP)
        END IF  
     END IF     
     
      !--GIMP shape functions and derivatives in the y direction
     IF(g_coord(2,nod)>=gm_coord(2,s))THEN
        fact=ydist  !--ydist positive (particle below the node)
        IF(fact<=lpGIMP)THEN    ! fact is the local coordinate of the material point
          Sz=one-((fact**2+lpGIMP**2)/(two*cellsize*lpGIMP))
         ELSE IF(fact>=lpGIMP.and.fact<=(cellsize-lpGIMP))THEN
          Sz=one-(fact/cellsize)
         ELSE IF(fact>=(cellsize-lpGIMP).and.fact<(cellsize+lpGIMP-tol))THEN
          Sz=((cellsize+lpGIMP-fact)**2)/(four*cellsize*lpGIMP)
        END IF   
     ELSE IF(g_coord(2,nod)<=gm_coord(2,s))THEN
        fact=ydist !--ydist negative (particle over the node)
        IF(fact>=(-lpGIMP))THEN    ! fact is the local coordinate of the material point 
         Sz=one-((fact**2+lpGIMP**2)/(two*cellsize*lpGIMP))
        ELSE IF(fact<=(-lpGIMP).and.fact>=(-(cellsize-lpGIMP)))THEN
         Sz=one+(fact/cellsize)
        ELSE IF(fact<=(-(cellsize-lpGIMP)).and.fact>(-(cellsize+lpGIMP-tol)))THEN
         Sz=((cellsize+lpGIMP+fact)**2)/(four*cellsize*lpGIMP)
        END IF  
     END IF 
    
!$$$$$$ 
!$$$$$$     !Shape functions
!$$$$$$     IF(n==1)THEN  !Comparing with node 1
!$$$$$$       Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
!$$$$$$       Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
!$$$$$$ 
!$$$$$$      ELSE IF(n==2)THEN    !Comparing with node 2
!$$$$$$       Sr=(one/(four*lp))*(two*Xi2-Xi2**2-two*Xi1+Xi1**2)
!$$$$$$       Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
!$$$$$$      
!$$$$$$      ELSE IF(n==3)THEN   !Comparing with node 3
!$$$$$$        Sr=(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
!$$$$$$        Sr=Sr*(one/(four*lp))
!$$$$$$        Sz=(one/(four*lp))*(two*ni2+ni2**2-two*ni1-ni1**2)
!$$$$$$      
!$$$$$$      ELSE IF(n==4)THEN      !Comparing with node 4
!$$$$$$       Sr=(one/(four*lp))*(two*Xi2+Xi2**2-two*Xi1-Xi1**2)
!$$$$$$       Sz=(one/(four*lp))*(two*ni2-ni2**2-two*ni1+ni1**2)
!$$$$$$     END IF



      
     IF(n==1)THEN  !Comparing with node 1
        
        dSr=(Xi1-Xi2)/(two*lp)
        IF(factx<(-one))dSr=-dSr

        dSz=(ni1-ni2)/(two*lp)
        IF(facty<(-one))dSz=-dSz


      ELSE IF(n==2)THEN    !Comparing with node 2
       
         dSr=(Xi1-Xi2)/(two*lp)
         IF(factx<(-one))dSr=-dSr

         dSz=(ni2-ni1)/(two*lp)
         IF(facty>(one))dSz=-dSz
 
      ELSE IF(n==3)THEN   !Comparing with node 3
         
        dSr=(Xi2-Xi1)/(two*lp)
        IF(factx>(one))dSr=-dSr

        dSz=(ni2-ni1)/(two*lp)
        IF(facty>(one))dSz=-dSz
 
      ELSE IF(n==4)THEN      !Comparing with node 4
       
        dSr=(Xi2-Xi1)/(two*lp)
        IF(factx>(one))dSr=-dSr

        dSz=(ni1-ni2)/(two*lp)
        IF(facty<(-one))dSz=-dSz
         
      END IF  

   IF(a_ele(s)/=iel.and.c_ele(iel)==0.and.Sr<5.0e-2)THEN!.and.(n==3.or.n==4).and.(Sr<4.0e-2))THEN !In case the node is far from the material point and the shape functions are small, a new value of 0.5 will be assigned to the hape functions
         Sr=0.05!;Sz=0.10
   ELSE IF(a_ele(s)/=iel.and.c_ele(iel)==0.and.abs(dSr)<1.0e-2)THEN
        dSr=(dSr/dSr)*0.05
  END IF


     funGIMP(m)=Sr*Sz
     IF((Sr<=0.0.or.Sz<=0.0))THEN
      derGIMP(1,m)=0.0
      derGIMP(2,m)=0.0
     ELSE  
       derGIMP(1,m)=dSr*Sz
       derGIMP(2,m)=dSz*Sr     
     END IF
    
     m=m+1
     
   END DO Fun_Deriv    

   RETURN
  END SUBROUTINE iGIMP_funder4

  SUBROUTINE deriv_select(i,k,iel,numel,elemmpoints,beeextend,funextend,km_w,fun,bee)

  IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::numel,i,elemmpoints(:,:),k,iel
   REAL(iwp),INTENT(IN)::beeextend(:,:),km_w(:),funextend(:)
   REAL(iwp),INTENT(OUT)::bee(:,:),fun(:)
   REAL(iwp)::p5=0.5_iwp,p25=0.25_iwp,zero=0.0_iwp

   bee=zero
   fun=zero

   IF(numel==1)THEN
     
     bee(:,1)=beeextend(:,3)
     bee(:,2)=beeextend(:,4)
     bee(:,3)=beeextend(:,1)
     bee(:,4)=beeextend(:,2)
     bee(:,5)=beeextend(:,5)
     bee(:,6)=beeextend(:,6)
     bee(:,7)=beeextend(:,7)
     bee(:,8)=beeextend(:,8)
     fun(1)=funextend(2)
     fun(2)=funextend(1)
     fun(3)=funextend(3)
     fun(4)=funextend(4)
     
   ELSE IF (numel==2.and.i==1.and.ABS(elemmpoints(k,i)-iel)>3.0)THEN 

     bee(:,1)=beeextend(:,3)*p5
     bee(:,2)=beeextend(:,4)*p5
     bee(:,3)=beeextend(:,1)
     bee(:,4)=beeextend(:,2)
     bee(:,5)=beeextend(:,7)
     bee(:,6)=beeextend(:,8)
     bee(:,7)=beeextend(:,9)*p5
     bee(:,8)=beeextend(:,10)*p5
     fun(1)=funextend(2)*p5
     fun(2)=funextend(1)
     fun(3)=funextend(4)
     fun(4)=funextend(5)*p5


   ELSE IF (numel==2.and.i==2.and.ABS(elemmpoints(k,i)-iel)>3.0)THEN   

     bee(:,1)=beeextend(:,5)
     bee(:,2)=beeextend(:,6)
     bee(:,3)=beeextend(:,3)*p5
     bee(:,4)=beeextend(:,4)*p5
     bee(:,5)=beeextend(:,9)*p5
     bee(:,6)=beeextend(:,10)*p5
     bee(:,7)=beeextend(:,11)
     bee(:,8)=beeextend(:,12)
     fun(1)=funextend(2)
     fun(2)=funextend(1)*p5
     fun(3)=funextend(4)*p5
     fun(4)=funextend(5)

   ELSE IF (numel==2.and.i==1.and.ABS(elemmpoints(k,i)-iel)<3.0)THEN 

     bee(:,1)=beeextend(:,3)
     bee(:,2)=beeextend(:,4)
     bee(:,3)=beeextend(:,1)
     bee(:,4)=beeextend(:,2)
     bee(:,5)=beeextend(:,5)*p5
     bee(:,6)=beeextend(:,6)*p5
     bee(:,7)=beeextend(:,7)*p5
     bee(:,8)=beeextend(:,8)*p5
     fun(1)=funextend(2)
     fun(2)=funextend(1)
     fun(3)=funextend(3)*p5
     fun(4)=funextend(4)*p5

   ELSE IF (numel==2.and.i==2.and.ABS(elemmpoints(k,i)-iel)<3.0)THEN 

     bee(:,1)=beeextend(:,7)*p5
     bee(:,2)=beeextend(:,8)*p5
     bee(:,3)=beeextend(:,5)*p5
     bee(:,4)=beeextend(:,6)*p5
     bee(:,5)=beeextend(:,9)
     bee(:,6)=beeextend(:,10)
     bee(:,7)=beeextend(:,11)
     bee(:,8)=beeextend(:,12)
     fun(1)=funextend(4)*p5
     fun(2)=funextend(3)*p5
     fun(3)=funextend(5)
     fun(4)=funextend(6)

   ELSE IF (numel==4.and.i==1)THEN   

     bee(:,1)=beeextend(:,3)*p5
     bee(:,2)=beeextend(:,4)*p5
     bee(:,3)=beeextend(:,1)
     bee(:,4)=beeextend(:,2)
     bee(:,5)=beeextend(:,7)*p5
     bee(:,6)=beeextend(:,8)*p5
     bee(:,7)=beeextend(:,9)*p25
     bee(:,8)=beeextend(:,10)*p25
     fun(1)=funextend(2)*p5
     fun(2)=funextend(1)
     fun(3)=funextend(4)*p5
     fun(4)=funextend(5)*p25

   ELSE IF (numel==4.and.i==2)THEN   

     bee(:,1)=beeextend(:,9)*p25
     bee(:,2)=beeextend(:,10)*p25
     bee(:,3)=beeextend(:,7)*p5
     bee(:,4)=beeextend(:,8)*p5
     bee(:,5)=beeextend(:,13)
     bee(:,6)=beeextend(:,14)
     bee(:,7)=beeextend(:,15)*p5
     bee(:,8)=beeextend(:,16)*p5
     fun(1)=funextend(5)*p25
     fun(2)=funextend(4)*p5
     fun(3)=funextend(7)
     fun(4)=funextend(8)*p5

   ELSE IF (numel==4.and.i==3)THEN   

     bee(:,1)=beeextend(:,5)
     bee(:,2)=beeextend(:,6)
     bee(:,3)=beeextend(:,3)*p5
     bee(:,4)=beeextend(:,4)*p5
     bee(:,5)=beeextend(:,9)*p25
     bee(:,6)=beeextend(:,10)*p25
     bee(:,7)=beeextend(:,11)*p5
     bee(:,8)=beeextend(:,12)*p5
     fun(1)=funextend(3)
     fun(2)=funextend(2)*p5
     fun(3)=funextend(5)*p25
     fun(4)=funextend(6)*p5

   ELSE IF (numel==4.and.i==4)THEN   

     bee(:,1)=beeextend(:,11)*p5
     bee(:,2)=beeextend(:,12)*p5
     bee(:,3)=beeextend(:,9)*p25
     bee(:,4)=beeextend(:,10)*p25
     bee(:,5)=beeextend(:,15)*p5
     bee(:,6)=beeextend(:,16)*p5
     bee(:,7)=beeextend(:,17)
     bee(:,8)=beeextend(:,18)
     fun(1)=funextend(6)*p5
     fun(2)=funextend(5)*p25
     fun(3)=funextend(8)*p5
     fun(4)=funextend(9)

   END IF  
  
  END SUBROUTINE deriv_select

 SUBROUTINE mpnumber(pointsupport,i,m)
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::pointsupport(:,:),i
   INTEGER,INTENT(OUT)::m
   INTEGER::j
   
   m=0
   DO j=1,9
     IF(pointsupport(i,j)>0)m=m+1
   END DO  
 

 END SUBROUTINE mpnumber

 
  SUBROUTINE GIMP_funel(s,nip,iel,g_coord,cellsize,gm_coord,g_num,gimptol,derGIMP,fun)
     !     
     ! Subroutine to compute GIMP shape functions and GIMP derivatives of the shape functions 
     !
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::g_num(:,:),s,nip,iel
   REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord(:,:),cellsize,gimptol
   REAL(iwp),INTENT(OUT)::derGIMP(:,:),fun(:)
   REAL(iwp)::lp,twolp,two=2.0_iwp,minx,miny,tol
   REAL(iwp)::Dr,rdist,ni,xpoint,ypoint,xdist,ydist,elfact,Sr,Sz,dSr,dSz,fact
   REAL(iwp)::four=4.0_iwp,three=3.0_iwp,twelve=12.0_iwp,six=6.0_iwp,one=1.0_iwp
   INTEGER::i,j,k,nod,side,nodes
   !LOGICAL,INTENT(IN)::gimptol


   derGIMP=0.0
   fun=0.0
   Sr=0.0
   Sz=0.0
   dSr=0.0
   dSz=0.0
   !tol=1.0e-8
   nodes=UBOUND(g_num,1)
   IF(nip==1)twolp=cellsize/(nip)
   IF(nip==4)twolp=cellsize/(nip/two)
   lp=twolp/two
   i=1
   k=1

   !IF(gimptol)THEN
   !  CONTINUE
   ! ELSE 
   !   tol=0.0
   !END IF
   
   tol=gimptol

   DO WHILE(i<=nodes)   ! 9 is the maximum of nodes inside a material point suport domain 
    !IF(GIMP_nodes(i,s)>0)THEN
     Sr=0.0
     Sz=0.0
     dSr=0.0
     dSz=0.0
     nod=g_num(i,iel)     !--node 'i' in the suport domain to interpolate value to the current material point 's' 
     minx=g_coord(1,nod)     !--minx coordinate of node 'i'
     xpoint=gm_coord(1,s)    !--xpoint coordinate of current material point 's'
     xdist=xpoint-minx       !--xdist distance between node 'i' and material point 's' (if material point is at the right side of the node, distance is positive)
     miny=g_coord(2,nod)     !--miny coordinate of node 'i'
     ypoint=gm_coord(2,s)    !--ypoint coordinate of material point 's'
     ydist=(-ypoint)-(-miny) !--ydist distance between node 'i' and material point 's' (if material point is below the node, distance is positive)
     
     !--GIMP shape functions and derivatives in the x direction
     IF(g_coord(1,nod)<=gm_coord(1,s))THEN
        fact=xdist
        IF(fact<=lp)THEN    
         Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact)/(cellsize*lp)
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sr=one-(fact/cellsize)
         dSr=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<(cellsize+lp-tol))THEN
         Sr=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSr=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(1,nod)>=gm_coord(1,s))THEN  
        fact=xdist 
        IF(fact>=(-lp))THEN    
          Sr=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSr=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sr=one+(fact/cellsize)
          dSr=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>(-(cellsize+lp-tol)))THEN
         Sr=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSr=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF     
     
      !--GIMP shape functions and derivatives in the y direction
     IF(g_coord(2,nod)>=gm_coord(2,s))THEN
        fact=ydist  !--ydist positive
        IF(fact<=lp)THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact>=lp.and.fact<=(cellsize-lp))THEN
         Sz=one-(fact/cellsize)
         dSz=-(one/cellsize)
        ELSE IF(fact>=(cellsize-lp).and.fact<(cellsize+lp-tol))THEN
         Sz=((cellsize+lp-fact)**2)/(four*cellsize*lp)
         dSz=-((cellsize+lp-fact)/(two*lp*cellsize))
        END IF   
     ELSE IF(g_coord(2,nod)<=gm_coord(2,s))THEN
        fact=ydist !--ydist negative
        IF(fact>=(-lp))THEN    ! fact is the local coordinate of the material point in AxisymmetricGIMP
         Sz=one-((fact**2+lp**2)/(two*cellsize*lp))
         dSz=-(fact/(cellsize*lp))
        ELSE IF(fact<=(-lp).and.fact>=(-(cellsize-lp)))THEN
         Sz=one+(fact/cellsize)
         dSz=one/cellsize
        ELSE IF(fact<=(-(cellsize-lp)).and.fact>(-(cellsize+lp-tol)))THEN
         Sz=((cellsize+lp+fact)**2)/(four*cellsize*lp)
         dSz=((cellsize+lp+fact)/(two*lp*cellsize))
        END IF  
     END IF  

     fun(k)=Sr*Sz
     derGIMP(1,k)=dSr*Sz
     derGIMP(2,k)=-dSz*Sr

     IF((dSr*Sz>0.0.or.-dSz*Sr>0.0).and.i==1)THEN
      fun(i)=0.0
     ELSE IF((dSr*Sz>0.0.or.-dSz*Sr<0.0).and.i==2)THEN
       fun(i)=0.0
     ELSE IF((dSr*Sz<0.0.or.-dSz*Sr<0.0).and.i==3)THEN  
       fun(i)=0.0
     ELSE IF((dSr*Sz<0.0.or.-dSz*Sr>0.0).and.i==4)THEN  
       fun(i)=0.0
     END IF  

     k=k+1
    !END IF
    i=i+1
   END DO 
   RETURN
  END SUBROUTINE GIMP_funel

  SUBROUTINE nor_dir_GIMP(normal,m,nf,GIMP_node_mp1,GIMP_node_mp2,gm_coord1,gm_coord2,g_coord)

    USE main
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::m,nf(:,:),GIMP_node_mp1(:,:),GIMP_node_mp2(:,:)
    REAL(iwp),INTENT(IN OUT)::normal(:,:)
    REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord1(:,:),gm_coord2(:,:)
    REAL(iwp)::magnorm,zero=0.0_iwp,div,normal_acum1(size(nf,1)),normal_acum2(size(nf,1)),normal_acum(size(nf,1))
    REAL(iwp)::Xval,Yval,cellsize,distX,distY
    INTEGER::i,ind,auxel,contel,ielint,nmps,k,mp,ndim,direction,minpoints
    INTEGER::iel1,iel2,countp,nmp1,nmp2

  !IF(sqrt((gm_coord1(1,i)-gm_coord2(1,j))**2+(gm_coord1(2,i)-gm_coord2(2,j))**2)<=distol)
    normal_acum=zero
    normal_acum1=zero
    normal_acum2=zero
    normal(:,m)=0.0
    cellsize=ABS(g_coord(2,1)-g_coord(2,2))*1.1
    
    Loop_1:DO i=1,SIZE(GIMP_node_mp1,1)
     IF(GIMP_node_mp1(i,m)==1)THEN
          distX=sqrt((g_coord(1,m)-gm_coord1(1,i))**2);distY=sqrt((g_coord(2,m)-gm_coord1(2,i))**2)
          IF(g_coord(1,m)>gm_coord1(1,i).and.distX<=cellsize)Xval=1.0
          IF(g_coord(1,m)<gm_coord1(1,i).and.distX<=cellsize)Xval=-1.0
          IF(g_coord(2,m)>gm_coord1(2,i).and.distY<=cellsize)Yval=1.0
          IF(g_coord(2,m)<gm_coord1(2,i).and.distY<=cellsize)Yval=-1.0
          normal_acum1(1)=normal_acum1(1)+Xval
          normal_acum1(2)=normal_acum1(2)+Yval
!$$$$$$           IF(nf(1,m)==0)normal_acum1(1)=0.0 !-This option erase the normal vector if the node is fix in this direction
!$$$$$$           IF(nf(2,m)==0)normal_acum1(2)=0.0
!$$$$$$           countp=countp+1
!$$$$$$           IF(countp==nmp1)EXIT
      END IF
    END DO Loop_1

    magnorm=SQRT(SUM(normal_acum1(:)**2))
    IF(magnorm==0)THEN
      normal_acum1=zero
    ELSE  
      normal_acum1=normal_acum1/magnorm
    END IF
    countp=0

    Loop_2:DO i=1,SIZE(GIMP_node_mp2,1)
     IF(GIMP_node_mp2(i,m)==1)THEN
          distX=sqrt((g_coord(1,m)-gm_coord2(1,i))**2);distY=sqrt((g_coord(2,m)-gm_coord2(2,i))**2)
          IF(g_coord(1,m)>gm_coord2(1,i).and.distX<=cellsize)Xval=1.0
          IF(g_coord(1,m)<gm_coord2(1,i).and.distX<=cellsize)Xval=-1.0
          IF(g_coord(2,m)>gm_coord2(2,i).and.distY<=cellsize)Yval=1.0
          IF(g_coord(2,m)<gm_coord2(2,i).and.distY<=cellsize)Yval=-1.0
          normal_acum2(1)=normal_acum2(1)+Xval
          normal_acum2(2)=normal_acum2(2)+Yval
!$$$$$$           IF(nf(1,m)==0)normal_acum2(1)=0.0 !-This option erase the normal vector if the node is fix in this direction
!$$$$$$           IF(nf(2,m)==0)normal_acum2(2)=0.0
!$$$$$$           countp=countp+1
!$$$$$$           IF(countp==nmp2)EXIT
      END IF
    END DO Loop_2

   normal_acum2=-normal_acum2
   magnorm=SQRT(SUM(normal_acum2(:)**2))
   IF(magnorm==0)THEN
     normal_acum2=zero
   ELSE  
   normal_acum2=normal_acum2/magnorm
   END IF

    normal_acum=normal_acum1+normal_acum2
     
    !-If magnorm>0.0 then the bodies are in different elemenmts, otherwise the bodies are in the same element
    magnorm=SQRT(SUM(normal_acum(:)**2))
    IF(magnorm>0.0)THEN
      normal(:,m)=normal_acum(:)/magnorm !Normalise normal vector
      normal_acum=normal_acum/magnorm
    END IF
     
  END SUBROUTINE nor_dir_GIMP


  SUBROUTINE Dist_Points(distact,m,GIMP_node_mp1,GIMP_node_mp2,gm_coord1,gm_coord2,g_coord)

    USE main
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::m,GIMP_node_mp1(:,:),GIMP_node_mp2(:,:)
    REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord1(:,:),gm_coord2(:,:)
    REAL(iwp)::cellsize
    INTEGER::i
    LOGICAL,INTENT(OUT)::distact
    LOGICAL::distbod1=.false.,distbod2=.false.


    distact=.false.
    cellsize=ABS(g_coord(2,2)-g_coord(2,1))*0.50_iwp ! The 0.5 is to concider the half of the cellsize for the contact
    
  Loop_dist_1:DO i=1,SIZE(GIMP_node_mp1,1)
   IF(GIMP_node_mp1(i,m)==1)THEN
     IF(ABS(g_coord(1,m)-gm_coord1(1,i))<=cellsize.and.ABS(g_coord(2,m)-gm_coord1(2,i))<=cellsize)distbod1=.true.
   END IF
   IF(distbod1)EXIT
  END DO Loop_dist_1

  

  Loop_dist_2:DO i=1,SIZE(GIMP_node_mp1,1)
   IF(GIMP_node_mp2(i,m)==1)THEN
     IF(ABS(g_coord(1,m)-gm_coord2(1,i))<=cellsize.and.ABS(g_coord(2,m)-gm_coord2(2,i))<=cellsize)distbod2=.true.
   END IF
   IF(distbod2)EXIT
  END DO Loop_dist_2

  IF(distbod1.and.distbod1)distact=.true.
  
     
  END SUBROUTINE Dist_Points


  SUBROUTINE Dist_Points2(distact,closedist,fardistance,distol,m,GIMP_node_mp1,GIMP_node_mp2,gm_coord1,gm_coord2,g_coord,closec,farc)

    USE main
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::m,GIMP_node_mp1(:,:),GIMP_node_mp2(:,:)
    REAL(iwp),INTENT(IN)::g_coord(:,:),gm_coord1(:,:),gm_coord2(:,:),distol,closedist,fardistance
    INTEGER,INTENT(IN OUT)::closec(:),farc(:)
    REAL(iwp)::cellsize,ncoordx,ncoordy
    INTEGER::i,j
    INTEGER,INTENT(IN OUT)::distact
    INTEGER::distbod1,distbod2


    distact=0
    distbod1=0
    distbod2=0
    cellsize=ABS(g_coord(2,2)-g_coord(2,1))
     ! The 0.5 is to concider the half of the cellsize for the contact
    
  Loop_dist_1:DO i=1,SIZE(GIMP_node_mp1,1)
   IF(GIMP_node_mp1(i,m)==1)THEN
     Loop_dist_2:DO j=1,SIZE(GIMP_node_mp2,1)
      IF(GIMP_node_mp2(j,m)==1)THEN
          ncoordx=g_coord(1,m);ncoordy=g_coord(2,m)
        IF(sqrt(((gm_coord1(1,i)-gm_coord2(1,j))**2+(gm_coord1(2,i)-gm_coord2(2,j))**2))<=distol)THEN
           distbod1=1
         IF(sqrt((ncoordx-gm_coord1(1,i))**2+(ncoordy-gm_coord1(2,i))**2)<=closedist.and.   &
           sqrt((ncoordx-gm_coord2(1,j))**2+(ncoordy-gm_coord2(2,j))**2)<=closedist)THEN
             closec(m)=m
           ELSE IF(sqrt((ncoordx-gm_coord1(1,i))**2+(ncoordy-gm_coord1(2,i))**2)>fardistance.and.   &
                   sqrt((ncoordx-gm_coord2(1,j))**2+(ncoordy-gm_coord2(2,j))**2)>fardistance.and.closec(m)==0)THEN
             farc(m)=m
         END IF
        END IF    
      END IF 
      !IF(distbod1==1.and.closec(m)>1)EXIT
     END DO Loop_dist_2
   END IF
   !IF(distbod1==1.and.closec(m)>1)EXIT
  END DO Loop_dist_1

  IF(distbod1==1)distact=1
  
     
  END SUBROUTINE Dist_Points2




 END MODULE GIMP
