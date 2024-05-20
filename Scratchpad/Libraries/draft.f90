mbod(bod)%kv=zero
    mbod(bod)%mv=zero
    mbod(bod)%kp=zero
    KM_MV:DO k=1,mbod(bod)%nmps
      ! This loop goes until 4 since 4 is the maximum number of 
      ! elements affected by a material point


El_stiff:DO i=1,4 
  iel = mbod(bod)%elemmpoints(k,i)
  ! If is true, an element is affected by a material point 
  ! (material point can be outside the element)
  Full_el:IF(mbod(bod)%elemmpoints(k,i)>0)THEN 
    num = g_num(:,iel)
    nod_num(:,1) = num
    coord = TRANSPOSE(g_coord(:,num))
    g = mbod(bod)%g_g(:,iel)

    km_gauss = zero   
    waverage = 4.0/nip

    !-Double mapping technique
    Double_map:IF(smethod>1)THEN 

      CALL iGIMP_funder3(                                 &
        k,mbod(bod)%mpoints,mbod(bod)%lp_mp,nip,coord,    &
        cellsize,mbod(bod)%gm_coord,mbod(bod)%gimptol,    &
        nod_num,g_coord,mbod(bod)%a_ele,mbod(bod)%c_ele,  &
        iel,der,fun                                       &
      )

      CALL sample(element,points,weights)
      DO s=1,nip
        CALL shape_der(der_gauss,points,s)
        CALL shape_fun(fun_gauss,points,s)
        scalar_dee = fun*fun_gauss*waverage
        sum_scalar = SUM(scalar_dee)
        dee_scal = mbod(bod)%dee*sum_scalar
        jac = MATMUL(der_gauss,coord)
        det = determinant(jac)
        CALL invert(jac)
        deriv_gauss = MATMUL(jac,der_gauss)
        CALL beemat(bee_gauss,deriv_gauss)
        IF(mbod(bod)%c_ele(iel)<1) waverage=1.5_iwp
        km_gauss = km_gauss + MATMUL(MATMUL(TRANSPOSE(bee_gauss),dee_scal),bee_gauss)*det*weights(i)
      END DO
      CALL fsparv(mbod(bod)%kv,km_gauss,g,mbod(bod)%kdiag)
    END IF Double_map
  END IF Full_el
END DO El_stiff






!-----------END STIFNESS MATRICES---------------

!---Diagonal mass matrix---  
values=mbod(bod)%valuesg(k)     
ALLOCATE(                                                           &
  derextend(nodof,values),funextend2(values*2,2),                 &
  eldddylds(values*2),mm_gimp(values*2,values*2)                  &
)


CALL GIMP_funder2(                                                  &
k,nip,g_coord,cellsize,mbod(bod)%gm_coord,mbod(bod)%lp_mp,      &
mbod(bod)%GIMP_nodes,mbod(bod)%gimptol,derextend,funextend2     &
)
mm_gimp = zero
mvval   = 2
CALL gimpfunform(                                                   &
  k,iel,eldddylds,mbod(bod)%nf,mbod(bod)%GIMP_nodes,              &
  values,mbod(bod)%g_g,mvval                                      &
) 

mvval=1
a=1
DO j=1,values
    DO q=1,2
        mm_gimp(a,a)=funextend2(a,q)*mbod(bod)%m_mass(k)
        a=a+1
    END DO
END DO
CALL fsparv(mbod(bod)%mv,mm_gimp,eldddylds,mbod(bod)%kdiag)
DEALLOCATE(derextend,funextend2,eldddylds,mm_gimp)
END DO KM_MV


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

  ALLOCATE(supdom(4,2))
  cellsize=g_coord(2,1)-g_coord(2,2)
  lpx=lp_mp(1,s)!+gimptol
  lpy=lp_mp(2,s)!+gimptol

  nn=UBOUND(g_coord,2)
  nmps=UBOUND(gm_coord,2)

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

  RETURN    
END SUBROUTINE GIMP_nodsup

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
  GIMP_node_mp=0

  nn=UBOUND(g_coord,2)
  nmps=UBOUND(gm_coord,2)
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

  END DO MP_Loop


  RETURN
END SUBROUTINE GIMP_activenode