PROGRAM p56         
!-------------------------------------------------------------------------
! Program 5.6 Three-dimensional strain of an elastic solid using
!             8-, 14- or 20-node brick hexahedra. Mesh numbered in x-z
!             planes then in the y-direction. No global stiffness matrix
!             assembly. Diagonally preconditioned conjugate gradient solver.
!-------------------------------------------------------------------------
 USE main 
 USE geom 
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15),npri=1,nstep=1
 INTEGER::cg_iters,cg_limit,fixed_freedoms,i,iel,k,loaded_nodes,ndim=3,   &
   ndof,nels,neq,nip,nlen,nn,nprops=2,np_types,nod,nodof=3,nr,nst=6,nxe,  &
   nye,nze
 REAL(iwp)::alpha,beta,cg_tol,det,dtim,one=1.0_iwp,penalty=1.0e20_iwp,up,      &
   zero=0.0_iwp 
   CHARACTER(LEN=15)::argv,element='hexahedron'
   LOGICAL::cg_converged,solid=.TRUE.
!-----------------------dynamic arrays------------------------------------
 INTEGER,ALLOCATABLE::etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),no(:),    &
   node(:),num(:),sense(:)
 REAL(iwp),ALLOCATABLE::bee(:,:),coord(:,:),d(:),dee(:,:),der(:,:),       &
   deriv(:,:),diag_precon(:),eld(:),fun(:),gc(:),g_coord(:,:),jac(:,:),   &
   km(:,:),loads(:),p(:),points(:,:),prop(:,:),sigma(:),store(:),         &
   storkm(:,:,:),u(:),value(:),weights(:),x(:),xnew(:),x_coords(:),       &
   y_coords(:),z_coords(:)
!-----------------------input and initialisation--------------------------
 CALL getname(argv,nlen)
 OPEN(10,FILE=argv(1:nlen)//'.dat') 
 OPEN(11,FILE=argv(1:nlen)//'.res')
 READ(10,*)nod,nxe,nye,nze,nip,cg_tol,cg_limit,np_types
 CALL mesh_size(element,nod,nels,nn,nxe,nye,nze) 
 ndof=nod*nodof
 ALLOCATE(nf(nodof,nn),points(nip,ndim),dee(nst,nst),coord(nod,ndim),     &
   jac(ndim,ndim),der(ndim,nod),deriv(ndim,nod),fun(nod),gc(ndim),        &
   bee(nst,ndof),km(ndof,ndof),eld(ndof),sigma(nst),g_coord(ndim,nn),     &
   g_num(nod,nels),weights(nip),num(nod),g_g(ndof,nels),x_coords(nxe+1),  &
   g(ndof),y_coords(nye+1),z_coords(nze+1),etype(nels),                   &
   prop(nprops,np_types),storkm(ndof,ndof,nels))
 READ(10,*)prop 
 etype=1 
 IF(np_types>1)READ(10,*)etype
 READ(10,*)x_coords,y_coords,z_coords
 nf=1 
 READ(10,*)nr,(k,nf(:,k),i=1,nr) 
 CALL formnf(nf) 
 neq=MAXVAL(nf)
 WRITE(11,'(A,I10,A)')" There are",neq," equations"
 ALLOCATE(p(0:neq),loads(0:neq),x(0:neq),xnew(0:neq),u(0:neq),            &
   diag_precon(0:neq),d(0:neq))
 CALL sample(element,points,weights) 
 diag_precon=zero
!----------element stiffness integration, storage and preconditioner------ 
 elements_1: DO iel=1,nels
   CALL hexahedron_xz(iel,x_coords,y_coords,z_coords,coord,num)
   CALL num_to_g(num,nf,g) 
   g_num(:,iel)=num
   g_coord(:,num)=TRANSPOSE(coord) 
   g_g(:,iel)=g
   CALL deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)))
   num=g_num(:,iel) 
   g=g_g(:,iel) 
   coord=TRANSPOSE(g_coord(:,num))
   km=zero
   gauss_pts_1: DO i=1,nip
     CALL shape_der(der,points,i)
     jac=MATMUL(der,coord) 
     det=determinant(jac) 
     CALL invert(jac)
     deriv=MATMUL(jac,der) 
     CALL beemat(bee,deriv)
     km=km+MATMUL(MATMUL(TRANSPOSE(bee),dee),bee)*det*weights(i)
   END DO gauss_pts_1
   storkm(:,:,iel)=km
   DO k=1,ndof 
     diag_precon(g(k))=diag_precon(g(k))+km(k,k) 
   END DO
 END DO elements_1  
!-----------------------invert the preconditioner and get starting loads--
 loads=zero 
 READ(10,*)loaded_nodes,(k,loads(nf(:,k)),i=1,loaded_nodes)
 READ(10,*)fixed_freedoms
 IF(fixed_freedoms/=0)THEN
   ALLOCATE(node(fixed_freedoms),sense(fixed_freedoms),                   &
     value(fixed_freedoms),no(fixed_freedoms),store(fixed_freedoms))
   READ(10,*)(node(i),sense(i),value(i),i=1,fixed_freedoms)
   DO  i=1,fixed_freedoms 
     no(i)=nf(sense(i),node(i)) 
   END DO
   diag_precon(no)=diag_precon(no)+penalty 
   loads(no)=diag_precon(no)*value
   store=diag_precon(no)
 END IF
 CALL mesh_ensi(argv,nlen,g_coord,g_num,element,etype,nf,loads(1:),      &
                nstep,npri,dtim,solid)
 diag_precon(1:)=one/diag_precon(1:) 
 diag_precon(0)=zero 
 d=diag_precon*loads 
 p=d 
 x=zero 
 cg_iters=0
!-----------------------pcg equation solution-----------------------------
 pcg: DO 
   cg_iters=cg_iters+1 
   u=zero
   elements_2: DO iel=1,nels
     g=g_g(:,iel) 
     km=storkm(:,:,iel) 
     u(g)=u(g)+MATMUL(km,p(g)) 
   END DO elements_2
   IF(fixed_freedoms/=0)u(no)=p(no)*store 
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
 WRITE(11,'(A,I5)')" Number of cg iterations to convergence was",cg_iters
 WRITE(11,'(/A)')"  Node   x-disp      y-disp      z-disp" 
 loads=xnew
 DO k=1,nn 
   WRITE(11,'(I6,3E12.4)')k,loads(nf(:,k)) 
 END DO
 CALL dismsh_ensi(argv,nlen,nstep,nf,xnew(1:))
!-----------------------recover stresses at nip integrating point---------
 nip=1 
 DEALLOCATE(points,weights) 
 ALLOCATE(points(nip,ndim),weights(nip))
 CALL sample(element,points,weights) 
 loads(0)=zero
 WRITE(11,'(/A,I2,A)')" The integration point (nip=",nip,") stresses are:"
 WRITE(11,'(A,/,A)')"    Element     x-coord     y-coord     z-coord",    &
   "    sig_x       sig_y       sig_z       tau_xy      tau_yz      tau_zx" 
 elements_3: DO iel=1,nels
   CALL deemat(dee,prop(1,etype(iel)),prop(2,etype(iel))) 
   num=g_num(:,iel)
   coord=TRANSPOSE(g_coord(:,num)) 
   g=g_g(:,iel) 
   eld=loads(g)
   gauss_pts_2: DO i=1,nip
     CALL shape_der(der,points,i) 
     CALL shape_fun(fun,points,i)
     gc=MATMUL(fun,coord) 
     jac=MATMUL(der,coord) 
     CALL invert(jac)
     deriv=MATMUL(jac,der) 
     CALL beemat(bee,deriv)
     sigma=MATMUL(dee,MATMUL(bee,eld)) 
     WRITE(11,'(I8,4X,3E12.4)')iel,gc
     WRITE(11,'(6E12.4)')sigma
   END DO gauss_pts_2
 END DO elements_3
STOP
END PROGRAM p56
