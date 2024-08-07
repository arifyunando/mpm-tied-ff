PROGRAM p613
!-------------------------------------------------------------------------
! Program 6.13 Three-dimensional strain analysis of an elastic-plastic
!              (Mohr-Coulomb) slope using 20-node hexahedra. Viscoplastic
!              strain method. No global stiffness matrix assembly.
!              Diagonally preconditioned conjugate gradient solver.
!-------------------------------------------------------------------------
 USE main
 USE geom
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER::cg_iters,cg_limit,cg_tot,i,iel,ifix,iters,iy,k,limit,ndim=3,    &
   ndof=60,nels,neq,nip=8,nn,nod=20,nodof=3,nprops=6,np_types,nsrf,nst=6, &
   nx1,nx2,ny1,ny2,nze,nlen
 REAL(iwp)::alpha,beta,cf,cg_tol,ddt,det,dq1,dq2,dq3,dsbar,dt=1.0e15_iwp, &
   d1,d4=4.0_iwp,d180=180.0_iwp,e,f,fmax,h1,h2,lode_theta,one=1.0_iwp,phi,&
   phif,pi,psi,psif,sigm,snph,start_dt=1.e15_iwp,s1,tnph,tnps,tol,        &
   two=2.0_iwp,up,v,w1,w2,zero=0.0_iwp
 CHARACTER(LEN=80)::argv,element='hexahedron' 
 LOGICAL::converged,cg_converged
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),num(:)    
 REAL(iwp),ALLOCATABLE::bdylds(:),bee(:,:),bload(:),coord(:,:),d(:),      &
   dee(:,:),der(:,:),deriv(:,:),devp(:),diag_precon(:),eld(:),eload(:),   &
   eps(:),erate(:),evp(:),evpt(:,:,:),flow(:,:),fun(:),gravlo(:),         &
   g_coord(:,:),jac(:,:),km(:,:),loads(:),m1(:,:),m2(:,:),m3(:,:),        &
   oldis(:),p(:),points(:,:),prop(:,:),sigma(:),srf(:),storkm(:,:,:),u(:),&
   weights(:),x(:),xnew(:)
!-----------------------input and initialisation--------------------------
 CALL getname(argv,nlen)
 OPEN(10,FILE=argv(1:nlen)//'.dat')
 OPEN(11,FILE=argv(1:nlen)//'.res')
!---(ifix=1) smooth-smooth(ifix=2) rough-smooth(ifix=3) rough-rough---
 READ(10,*)w1,s1,w2,h1,h2,d1,nx1,nx2,ny1,ny2,nze,ifix,                    &
   cg_tol,cg_limit,np_types
 nels=(nx1*ny1+ny2*(nx1+nx2))*nze
 nn=((3*(ny1+ny2)+2)*nx1+2*(ny1+ny2)+1+(3*ny2+2)*nx2)*(1+nze)+            &
   ((ny1+ny2+1)*(nx1+1)+(ny2+1)*nx2)*nze
 ALLOCATE(nf(nodof,nn),points(nip,ndim),weights(nip),g_coord(ndim,nn),    &
   num(nod),dee(nst,nst),evpt(nst,nip,nels),coord(nod,ndim),fun(nod),     &
   g_g(ndof,nels),jac(ndim,ndim),der(ndim,nod),etype(nels),               &
   deriv(ndim,nod),g_num(nod,nels),bee(nst,ndof),km(ndof,ndof),eld(ndof), &
   eps(nst),sigma(nst),bload(ndof),eload(ndof),erate(nst),evp(nst),       &
   devp(nst),g(ndof),m1(nst,nst),m2(nst,nst),m3(nst,nst),flow(nst,nst),   &
   prop(nprops,np_types),storkm(ndof,ndof,nels))
 READ(10,*)prop
 etype=1
 IF(np_types>1)READ(10,*)etype
 CALL emb_3d_bc(ifix,nx1,nx2,ny1,ny2,nze,nf)
 neq=MAXVAL(nf)
 WRITE(11,'(A,I7,A)')" There are",neq," equations"
 ALLOCATE(loads(0:neq),bdylds(0:neq),oldis(0:neq),gravlo(0:neq),p(0:neq), &
   x(0:neq),xnew(0:neq),u(0:neq),diag_precon(0:neq),d(0:neq))
!-----------------------loop the elements to find global array sizes-----
 elements_1: DO iel=1,nels
   CALL emb_3d_geom(iel,nx1,nx2,ny1,ny2,nze,w1,s1,w2,h1,h2,d1,coord,num)
   g_num(:,iel)=num
   CALL num_to_g(num,nf,g)
   g_coord(:,num)=TRANSPOSE(coord)
   g_g(:,iel)=g
 END DO elements_1
 pi=ACOS(-one)
 oldis=zero
 gravlo=zero
 p=zero
 xnew=zero 
 diag_precon=zero
 CALL sample(element,points,weights)
!----------element stiffness integration, storage and preconditioner------ 
 elements_2: DO iel=1,nels
   CALL deemat(dee,prop(5,etype(iel)),prop(6,etype(iel)))
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num)) 
   g=g_g(:,iel)
   km=zero
   eld=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_fun(fun,points,i)
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord)
     det=determinant(jac)
     CALL invert(jac)
     deriv=MATMUL(jac,der)
     CALL beemat(bee,deriv)           
     km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
     eld(2:ndof-1:3)=eld(2:ndof-1:3)+fun(:)*det*weights(i)
   END DO gauss_pts_1   
   storkm(:,:,iel)=km 
   DO k=1,ndof
     diag_precon(g(k))=diag_precon(g(k))+km(k,k)
   END DO
   gravlo(g)=gravlo(g)-eld*prop(4,etype(iel))
 END DO elements_2                                                            
 diag_precon(1:)=one/diag_precon(1:)
 diag_precon(0)=zero
!-----------------------trial strength reduction factor loop--------------
 READ(10,*)tol,limit,nsrf
 ALLOCATE(srf(nsrf))
 READ(10,*)srf
 WRITE(11,'(/A)')"    srf   max disp  iters     cg iters/plastic iter"
 srf_trials: DO iy=1,nsrf
   dt=start_dt
   DO i=1,np_types
     phi=prop(1,i)
     tnph=TAN(phi*pi/d180)
     phif=ATAN(tnph/srf(iy))
     snph=SIN(phif)
     e=prop(5,i)
     v=prop(6,i)
     ddt=d4*(one+v)*(one-two*v)/(e*(one-two*v+snph**2))
     IF(ddt<dt)dt=ddt
   END DO
   iters=0
   bdylds=zero
   evpt=zero
   cg_tot=0
!-----------------------plastic iteration loop----------------------------
   its: DO    
     fmax=zero
     iters=iters+1
     loads=gravlo+bdylds
     d=diag_precon*loads
     p=d
     x=zero
     cg_iters=0
!-----------------------pcg equation solution-----------------------------    
     pcg: DO
       cg_iters=cg_iters+1
       u=zero
       elements_3: DO iel=1,nels
         CALL deemat(dee,prop(2,etype(iel)),prop(3,etype(iel)))
         g=g_g(:,iel)
         km=storkm(:,:,iel)
         u(g)=u(g)+MATMUL(km,p(g))
       END DO elements_3
       up=DOT_PRODUCT(loads,d)
       alpha=up/DOT_PRODUCT(p,u)
       xnew=x+p*alpha
       loads=loads-u*alpha
       d=diag_precon*loads
       beta=DOT_PRODUCT(loads,d)/up
       p=d+p*beta
       CALL checon(xnew,x,cg_tol,cg_converged)
       IF(cg_converged.OR.cg_iters==cg_limit)EXIT
     END DO pcg
     cg_tot=cg_tot+cg_iters
     loads=xnew
     loads(0)=zero
!-----------------------check plastic convergence-------------------------
     CALL checon(loads,oldis,tol,converged)
     IF(iters==1)converged=.FALSE.
     IF(converged.OR.iters==limit)bdylds=zero
!-----------------------go round the Gauss Points ------------------------
     elements_4: DO iel=1,nels
       phi=prop(1,etype(iel))
       tnph=TAN(phi*pi/d180)
       phif=ATAN(tnph/srf(iy))*d180/pi
       psi=prop(3,etype(iel))
       tnps=TAN(psi*pi/d180)
       psif=ATAN(tnps/srf(iy))*d180/pi
       cf=prop(2,etype(iel))/srf(iy)
       e=prop(5,etype(iel))
       v=prop(6,etype(iel))
       CALL deemat(dee,e,v)
       num=g_num(:,iel)
       coord=TRANSPOSE(g_coord(:,num))
       g=g_g(:,iel)
       eld=loads(g)         
       bload=zero
       gauss_points_2: DO i=1,nip
         CALL shape_der(der,points,i) 
         jac=MATMUL(der,coord)
         det=determinant(jac)
         CALL invert(jac)
         deriv=MATMUL(jac,der)
         CALL beemat(bee,deriv)
         eps=MATMUL(bee,eld)
         eps=eps-evpt(:,i,iel)
         sigma=MATMUL(dee,eps)
         CALL invar(sigma,sigm,dsbar,lode_theta)                             
!-----------------------check whether yield is violated-------------------
         CALL mocouf(phif,cf,sigm,dsbar,lode_theta,f)
         IF(f>fmax)fmax=f
         IF(converged.OR.iters==limit)THEN
           devp=sigma
         ELSE
           IF(f>=zero.OR.(converged.OR.iters==limit))THEN
             CALL mocouq(psif,dsbar,lode_theta,dq1,dq2,dq3)
             CALL formm(sigma,m1,m2,m3)
             flow=f*(m1*dq1+m2*dq2+m3*dq3) 
             erate=MATMUL(flow,sigma)
             evp=erate*dt
             evpt(:,i,iel)=evpt(:,i,iel)+evp
             devp=MATMUL(dee,evp) 
           END IF
         END IF
         IF(f>=zero.OR.(converged.OR.iters==limit))THEN
           eload=MATMUL(devp,bee)
           bload=bload+eload*det*weights(i)
         END IF
       END DO gauss_points_2
!-----------------------compute the total bodyloads vector----------------
       bdylds(g)=bdylds(g)+bload
       bdylds(0)=zero
     END DO elements_4             
     WRITE(*,'(A,F7.2,A,I4,A,F8.3)')                                      &
       "  srf",srf(iy),"  iteration",iters,"  F_max",fmax
     IF(converged.OR.iters==limit)EXIT
   END DO its
   WRITE(11,'(F7.2,E12.4,I5,F17.2)')                                     &
     srf(iy),MAXVAL(ABS(loads)),iters,REAL(cg_tot)/REAL(iters)
   IF(iters==limit)EXIT
 END DO srf_trials
STOP
END PROGRAM p613























