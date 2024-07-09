module mpm_gimp
    
    save
    contains

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


SUBROUTINE eldform(eld,loads,iel,neighb,values,g_g,bound,m)
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    INTEGER,INTENT(IN)::iel,neighb(:,:),values,bound,g_g(:,:),m
    REAL(iwp),INTENT(IN)::loads(:)
    REAL(iwp),INTENT(OUT)::eld(:)
    INTEGER,ALLOCATABLE::b(:)
    INTEGER::i,n
 
    ALLOCATE(b(4))
    eld=0.0_iwp
    n=1
 
    IF((bound==1.or.bound==2.or.bound==4).and.m==1)THEN
        !bound 1=upper left corner; bound 2=upper boundarie; bound 4=left boundarie
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,4))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,5))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,6))
        b(2)=g_g(2,neighb(iel,6))
        b(3)=g_g(7,neighb(iel,6))
        b(4)=g_g(8,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,5))
        b(2)=g_g(2,neighb(iel,5))
        b(3)=g_g(7,neighb(iel,5))
        b(4)=g_g(8,neighb(iel,5))
        eld(n:n+3)=loads(b+1)

    ELSE IF(bound==2.and.m==2)THEN
        !bound 1=upper left corner; bound 2=upper boundarie; bound 4=left boundarie
        b=g_g(3:6,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,7))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,7))
        b(2)=g_g(2,neighb(iel,7))
        b(3)=g_g(7,neighb(iel,7))
        b(4)=g_g(8,neighb(iel,7))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,6))
        b(2)=g_g(2,neighb(iel,6))
        b(3)=g_g(7,neighb(iel,6))
        b(4)=g_g(8,neighb(iel,6))
        eld(n:n+3)=loads(b+1)

    ELSE IF(bound==4.and.m==2)THEN
        !bound 1=upper left corner; bound 2=upper boundarie; bound 4=left boundarie
        b=g_g(3:6,neighb(iel,2))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,3))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,4))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,iel)
        b(2)=g_g(2,iel)
        b(3)=g_g(7,iel)
        b(4)=g_g(8,iel)
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,4))
        b(2)=g_g(2,neighb(iel,4))
        b(3)=g_g(7,neighb(iel,4))
        b(4)=g_g(8,neighb(iel,4))
        eld(n:n+3)=loads(b+1)

    ELSE IF((bound==3.or.bound==6).and.m==1)THEN
        !bound 3=lowest left corner; bound 6=upper boundarie
        b=g_g(3:6,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,7))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,7))
        b(2)=g_g(2,neighb(iel,7))
        b(3)=g_g(7,neighb(iel,7))
        b(4)=g_g(8,neighb(iel,7))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,6))
        b(2)=g_g(2,neighb(iel,6))
        b(3)=g_g(7,neighb(iel,6))
        b(4)=g_g(8,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
    ELSE IF(bound==6.and.m==2)THEN
        !bound 3=lowest left corner; bound 6=upper boundarie
        b=g_g(3:6,neighb(iel,1))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,2))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,8))
        b(2)=g_g(2,neighb(iel,8))
        b(3)=g_g(7,neighb(iel,8))
        b(4)=g_g(8,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,iel)
        b(2)=g_g(2,iel)
        b(3)=g_g(7,iel)
        b(4)=g_g(8,iel)
        eld(n:n+3)=loads(b+1)
    ELSE IF((bound==7.or.bound==8).and.m==1)THEN
        !bound 7=upper right left corner; bound 6=right boundarie
        b=g_g(3:6,neighb(iel,2))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,3))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,4))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,iel)
        b(2)=g_g(2,iel)
        b(3)=g_g(7,iel)
        b(4)=g_g(8,iel)
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,4))
        b(2)=g_g(2,neighb(iel,4))
        b(3)=g_g(7,neighb(iel,4))
        b(4)=g_g(8,neighb(iel,4))
        eld(n:n+3)=loads(b+1)
    ELSE IF(bound==8.and.m==2)THEN
        !bound 7=upper right left corner; bound 6=right boundarie
        b=g_g(3:6,neighb(iel,1))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,2))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,8))
        b(2)=g_g(2,neighb(iel,8))
        b(3)=g_g(7,neighb(iel,8))
        b(4)=g_g(8,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,iel)
        b(2)=g_g(2,iel)
        b(3)=g_g(7,iel)
        b(4)=g_g(8,iel)
        eld(n:n+3)=loads(b+1)
    ELSE IF(bound==9)THEN
        !bound 9= right lowest corner
        b=g_g(3:6,neighb(iel,1))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,2))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,8))
        b(2)=g_g(2,neighb(iel,8))
        b(3)=g_g(7,neighb(iel,8))
        b(4)=g_g(8,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,iel)
        b(2)=g_g(2,iel)
        b(3)=g_g(7,iel)
        b(4)=g_g(8,iel)
        eld(n:n+3)=loads(b+1)
    ELSE IF(bound==5)THEN
        !bound 5= center element
        b=g_g(3:6,neighb(iel,1))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,2))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,3))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,8))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,iel)
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,4))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b=g_g(3:6,neighb(iel,7))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b=g_g(3:6,neighb(iel,5))
        eld(n:n+3)=loads(b+1)
        n=n+4
        b(1)=g_g(1,neighb(iel,7))
        b(2)=g_g(2,neighb(iel,7))
        b(3)=g_g(7,neighb(iel,7))
        b(4)=g_g(8,neighb(iel,7))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,6))
        b(2)=g_g(2,neighb(iel,6))
        b(3)=g_g(7,neighb(iel,6))
        b(4)=g_g(8,neighb(iel,6))
        eld(n:n+3)=loads(b+1)
        n=n+2
        b(1)=g_g(1,neighb(iel,5))
        b(2)=g_g(2,neighb(iel,5))
        b(3)=g_g(7,neighb(iel,5))
        b(4)=g_g(8,neighb(iel,5))
        eld(n:n+3)=loads(b+1)
    END IF
    RETURN
END SUBROUTINE eldform


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
    LOGICAL::conditions, condx, condy
    elemmpoints=0
    nmps=UBOUND(gm_coord,2)
     
    tol=gimptol

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
                conditions=.false.; condx=.false.; condy=.false.
                ! Check if the element is active (4 nodes of the element free at 
                ! least in one direction or with a material point inside)
                conditions =  (                                                 &
                        (nf(1,g_num(1,i))>0.or.nf(2,g_num(1,i))>0).and.         &
                        (nf(1,g_num(2,i))>0.or.nf(2,g_num(2,i))>0).and.         &
                        (nf(1,g_num(3,i))>0.or.nf(2,g_num(3,i))>0).and.         &
                        (nf(1,g_num(4,i))>0.or.nf(2,g_num(4,i))>0)              &
                    ) .or. c_ele(i)>0
                ActiveEl:IF(conditions) THEN 
                    m=g_num(3,i)
                    n=g_num(2,i)
                    
                    !Chek if the material point is inside the x suport doman of one element
                    condx = (gm_coord(1,k)<=g_coord(1,m)+lpx-tol).and.(gm_coord(1,k)>=g_coord(1,n)-lpx+tol)
                    IF(condx)THEN 
                        m=g_num(2,i);n=g_num(1,i)
                        condy = (gm_coord(2,k)<=g_coord(2,m)+lpy-tol).and.(gm_coord(2,k)>=g_coord(2,n)-lpy+tol)
                        IF(condy)THEN !If is inside the x support domain, then chek if the material point is inside the y suport doman of one element
                            elemmpoints(k,cont)=i
                            cont=cont+1
                            IF(cont==5) THEN
                                GOTO 20 ! exit element loop
                            END IF
                        END IF
                    END IF
                END IF ActiveEl
            END DO Ele
            20 CONTINUE
        END DO
    END IF  
    
    RETURN
END SUBROUTINE Elem_suport


SUBROUTINE iGIMP_funder3(s,mpoints,lp_mp,nip,coord,cellsize,gm_coord,gimptol,GIMP_nodes,  &
    g_coord,a_ele,c_ele,iel,derGIMP,funGIMP)
    !     
    ! Subroutine to compute implicit GIMP shape functions and shape functions derivatives. 
    ! iGIMP_funder3 uses a constant gradient at the center of the nodes and it reduces as the
    ! support domain ends
    !
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

    derGIMP=0.0
    funGIMP=0.0
    i=1
    k=1
    m=1
    tol=gimptol

    xpoint=gm_coord(1,s)    !--x coordinate of current material point 's'
    ypoint=gm_coord(2,s)    !--y coordinate of material point 's'
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

        !--Evaluate the local coordinate of the material point inside and outside the element affected
        IF(a_ele(s)==iel)THEN !--Material point inside the element
            factx=mpoints(s,1)
        ELSE     !--Material pooint outside the element
            IF(xpoint<=coordx)THEN          !--point outside the element in the x direction (left side)
                IF(n==1.or.n==2) xdist=coordx-xpoint
                IF(n==3.or.n==4) xdist=(coordx-cellsize)-xpoint
                factx = ((cellsize/two+xdist)/(cellsize/two))*(-one)   !--local coordinate point at the left side of the element
            ELSE IF(xpoint>=coordx)THEN     !--point outside the element in the x direction (right side)
                IF(n==1.or.n==2) xdist=xpoint-(coordx+cellsize)
                IF(n==3.or.n==4) xdist=xpoint-(coordx)
                factx = ((cellsize/two+xdist)/(cellsize/two))       !--local coordinate point at the right side of the element
            END IF
        END IF

        IF(a_ele(s)==iel)THEN
            facty=mpoints(s,2)
        ELSE 
            IF(ypoint>=coordy)THEN          !--point outside the element in the y direction (above the element)
                IF(n==1.or.n==4)ydist=((coordy+cellsize)-ypoint)*(-one)
                IF(n==2.or.n==3)ydist=((coordy)-ypoint)*(-one)
                facty=((cellsize/two+ydist)/(cellsize/two))          !--local coordinate point over the element
            ELSE IF(ypoint<=coordy)THEN     !--point outside the element in the y direction (below the element)
                IF(n==2.or.n==3)ydist=(ypoint-(coordy-cellsize))*(-one)
                IF(n==1.or.n==4)ydist=(ypoint-(coordy))*(-one)
                facty=((cellsize/two+ydist)/(cellsize/two))*(-one)   !--local coordinate point below the element
            END IF  
        END IF
        !--End of the evaluation of the local coordinate of the material point

        !Shape functions
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

        IF(n==1)THEN        !Comparing with node 1
            dSr=(Xi1-Xi2)/(two*lp)
            dSz=(ni1-ni2)/(two*lp)
        ELSE IF(n==2)THEN   !Comparing with node 2
            dSr=(Xi1-Xi2)/(two*lp)
            dSz=(ni2-ni1)/(two*lp)
        ELSE IF(n==3)THEN   !Comparing with node 3
            dSr=(Xi2-Xi1)/(two*lp)
            dSz=(ni2-ni1)/(two*lp)
        ELSE IF(n==4)THEN   !Comparing with node 4
            dSr=(Xi2-Xi1)/(two*lp)
            dSz=(ni1-ni2)/(two*lp)
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
        END IF

        m=m+1

    END DO Fun_Deriv

    RETURN
END SUBROUTINE iGIMP_funder3


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

    derGIMP=0.0
    funGIMP=0.0
    Sr=0.0
    Sz=0.0
    dSr=0.0
    dSz=0.0
    nodes=UBOUND(GIMP_nodes,1)

    lpx=lp_mp(1,s);lpy=lp_mp(2,s)
    i=1
    k=1
    
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
    
    funGIMP=0.0
    derGIMP=0.0
    Sr=0.0
    Sz=0.0
    dSr=0.0
    dSz=0.0
    lpx=lp_mp(1,s)
    lpy=lp_mp(2,s)
    i=1
    k=1
    m=1

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
    iel=a_ele(s)        ! obtain particle residing cell
    neig=neighb(iel,1)  ! obatain cell id of the neighbours
    nod=g_num(2,neig)   ! 
    mainnod=g_num(2,neig) ! 
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


SUBROUTINE paraview2(input,realisation,argv,g_coord,g_num,nf,nels,nod,nn,nlen,  &
    diag,ddylds,d1x1,d2x1,gravlo,loads,normals,fcont,kv,mv,kdiag,vcm,f_fint)

    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::g_coord(:,:),normals(:,:),kv(:),mv(:)
    REAL(iwp),INTENT(IN)::diag(:),ddylds(:),loads(:),gravlo(:),d1x1(:),d2x1(:),fcont(:),vcm(:),f_fint(:)
    INTEGER,INTENT(IN)::input,nels,nod,nn,nlen,g_num(:,:),nf(:,:),realisation,kdiag(:)
    CHARACTER(*),INTENT(IN)::argv
    INTEGER::i,iel,ss
    REAL(iwp):: zero=0.0_iwp,dis_vec(2,nn)
    character (len=8) ::  cnumber,cnumber1
    ss=input!+100
    write(cnumber,'(i8.6)') ss
    write(cnumber1,'(i8.6)') realisation
    ss=15
    !open(ss,FILE = argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
    OPEN(ss,FILE="Output/Paraview2_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
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

        CASE DEFAULT
            WRITE(*,*)"wrong number of nodes input in paraview"
    end select

    WRITE(ss,'(a)')""
    WRITE(ss,'(1A10,1I5)')"POINT_DATA", nn
    WRITE(ss,'(a)')"vectors mass float "
    DO i=1, nn
        IF(nf(1,i)==0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')zero,zero, zero
        IF(nf(1,i)>0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')diag(nf(1,i)+1),zero,zero
        IF(nf(1,i)==0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')zero,diag(nf(2,i)+1),zero
        IF(nf(1,i)>0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')diag(nf(1,i)+1),diag(nf(2,i)+1),zero
    END DO

    WRITE(ss,'(a)')"vectors Ftotal float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')loads(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors Fint float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')ddylds(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors Fext float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')gravlo(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors fkin float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')vcm(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors ffint float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')f_fint(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors normals float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')normals(1,i),normals(2,i), zero
    END DO 

    WRITE(ss,'(a)')"vectors Fcont float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')fcont(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors Vel float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')d1x1(nf(:,i)+1), zero
    END DO

    WRITE(ss,'(a)')"vectors Acc float "
    DO i=1, nn
        WRITE(ss,'(3f15.6)')d2x1(nf(:,i)+1), zero
    END DO 

    WRITE(ss,'(a)')"vectors KM float "
    DO i=1, nn
        IF(nf(1,i)==0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')zero,zero, zero
        IF(nf(1,i)>0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')kv(kdiag(nf(1,i))),zero, zero
        IF(nf(1,i)==0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')zero,kv(kdiag(nf(2,i))), zero
        IF(nf(1,i)>0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')kv(kdiag(nf(1,i))),kv(kdiag(nf(2,i))), zero
    END DO

    WRITE(ss,'(a)')"vectors Mv float "
    DO i=1, nn
        IF(nf(1,i)==0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')zero,zero, zero
        IF(nf(1,i)>0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')mv(kdiag(nf(1,i))),zero, zero
        IF(nf(1,i)==0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')zero,mv(kdiag(nf(2,i))), zero
        IF(nf(1,i)>0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')mv(kdiag(nf(1,i))),mv(kdiag(nf(2,i))), zero
    END DO

    close(ss)
    RETURN
END SUBROUTINE paraview2

end module mpm_gimp