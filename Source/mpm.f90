module mpm

     SAVE
     CONTAINS

     SUBROUTINE locmp(iel,nipslo,m,n,mpart,npart,nx1,ny1,nx2,ny2,nyp,coord,m_coord,m_num,points,shape,slope,dir)

     ! This soubrutine locate the material points inside the mesh

      USE MAIN
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::coord(:,:),points(:,:)
      INTEGER,INTENT(IN)::m,n,mpart,npart,iel,nx1,ny1,nyp,nipslo
      LOGICAL,INTENT(IN):: shape,slope
      CHARACTER(LEN=1),INTENT(IN)::dir
      REAL(iwp),INTENT(OUT)::m_coord(:,:)
      INTEGER,INTENT(OUT)::m_num(:)
      REAL(iwp),ALLOCATABLE::s_coord(:,:),fun(:)
      INTEGER:: a(0:ny1)
      INTEGER::i,nip,ndim,ip,iq,nod,nx2,nx3,ny2,wt1,wt2,nxe,nye,nm1,nm2
      REAL(iwp)::aa(2),bb(2),cc(2),dd(2),pt5=0.5_iwp
      ! ===============================
      ! the added array n_z is for random part. number of elements in each row.
      !===============
      !   this part is changed specifically for the column problem, with a  ͹ shape.
         ! nx2=4;nx3=20;ny2=4
         ! the numbering direction for the elements and material points should be consistent.
      !===============
       IF(nipslo>0)nip=nipslo
       IF(nipslo==0)nip=UBOUND(points,1)
       !nip=UBOUND(points,1)
       ndim=UBOUND(points,2)
       nod=UBOUND(coord,1)
       !---------------------
       ALLOCATE(fun(nod),s_coord(nod,ndim))

      DO i=1,nip
       CALL shape_fun(fun,points,i)
       m_coord(i,:)=MATMUL(fun,coord)
      END DO
      DEALLOCATE(fun,s_coord)


      RETURN
      END SUBROUTINE locmp

     SUBROUTINE floc(coord,gm_coord,mpoints,i)

     ! local coordinate of the material points

      USE MAIN
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:)
      INTEGER,INTENT(IN)::i
      REAL(iwp),INTENT(IN OUT)::mpoints(:,:)
      REAL(iwp),ALLOCATABLE::guess_local(:,:),fun(:),guess_global(:,:), &
      der(:,:),jac(:,:)
      INTEGER::count,nod,ndim,nmps
      REAL(iwp)::tolerance1,tolerance2,beta,zero=0.0_iwp,pt9=0.9_iwp,one=1.0_iwp, &
                 tolerance
      LOGICAL:: find_local
        nod=UBOUND(coord,1)
        ndim=UBOUND(coord,2)
        nmps=UBOUND(gm_coord,2)
      ! --------------------
      ALLOCATE(fun(nod),guess_local(1,ndim),guess_global(1,ndim), &
      der(ndim,nod),jac(ndim,ndim))
      find_local=.FALSE.
      guess_local=zero
      CALL shape_fun(fun,guess_local,1)
      guess_global(1,:)=MATMUL(fun,coord)
      count=0
      DO
         count=count+1
         beta=one
         IF(count>=1) beta=pt9
         guess_global(1,:)=gm_coord(:,1)-guess_global(1,:)
         CALL shape_der(der,guess_local,1)
         jac=MATMUL(der,coord)
         CALL invert(jac)
     !guess_local(i,:)=guess_local(i,:)+beta*MATMUL(guess_global(i,:),jac)
         guess_local(1,:)=guess_local(1,:)+beta*MATMUL(TRANSPOSE(jac),guess_global(1,:))
         CALL shape_fun(fun,guess_local,1)
         guess_global(1,:)=MATMUL(fun,coord)
         tolerance1=ABS(MAXVAL(guess_global(1,:)-gm_coord(:,1)))
         tolerance2=ABS(MINVAL(guess_global(1,:)-gm_coord(:,1)))
         !find_local=(tolerance1<=1.0D-6)
         !find_local=(tolerance2<=1.0D-6)
         !IF(find_local) THEN
         IF(tolerance1<=1.0D-10.and.tolerance2<=1.0D-10)THEN
             mpoints(1,:)=guess_local(1,:)
             EXIT
         END IF
         IF(count>10) THEN
             write(*,*)" wrong numerical method"
             PAUSE
             EXIT
         END IF
      END DO
      DEALLOCATE(fun,guess_local,guess_global,der,jac)
      RETURN
      END SUBROUTINE floc

      SUBROUTINE paraview(input,realisation,argv,g_coord,g_num,nf,nels,nod,nn,nlen,  &
              diag,ddylds,d1x1,d2x1,gravlo,loads,normals,fcont,kv,mv,kdiag,vcm,f_fint)

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::g_coord(:,:),normals(:,:),kv(:),mv(:)
      REAL(iwp),INTENT(IN)::diag(:),ddylds(:),loads(:),gravlo(:),d1x1(:),d2x1(:),fcont(:),vcm(:),f_fint(:)
      INTEGER,INTENT(IN)::input,nels,nod,nn,nlen,g_num(:,:),nf(:,:),realisation,kdiag(:)
      CHARACTER(*),INTENT(IN)::argv
      INTEGER::i,iel,ss
      REAL(iwp):: zero=0.0_iwp,dis_vec(2,nn)
      character (len=6) ::  cnumber,cnumber1
      ss=input!+100
      write(cnumber,'(i6.6)') ss
      write(cnumber1,'(i6.6)') realisation
      ss=15
       !open(ss,FILE = argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
       OPEN(ss,FILE="Paraview1_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
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

!       case(8)
!         WRITE(ss,'(1A5,2I5)')"CELLS", nels, nels*(1+8)
!       DO iel=1, nels
!         WRITE(ss,'(9I5)')nod,  &
!          g_num(3,iel)-1, g_num(5,iel)-1, g_num(7,iel)-1, g_num(1,iel)-1, &
!                        g_num(4,iel)-1, g_num(6,iel)-1, g_num(8,iel)-1, g_num(2,iel)-1
!       END DO

!       WRITE(ss,'(a)')""
!       WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
!       DO iel = 1 , nels
!         WRITE(ss,*)"23"
!         END DO

       CASE DEFAULT
          WRITE(*,*)"wrong number of nodes input in paraview"
       end select

       WRITE(ss,'(a)')""
       !
       WRITE(ss,'(1A10,1I5)')"POINT_DATA", nn
     !!
     !  WRITE(ss,'(a)')"vectors displacement float "
     !  DO i=1,nn
     !  dis_vec(:,i)=dis(nf(:,i))
     !  WRITE(ss,'(3f9.4)') dis_vec(:,i), zero
     !  END DO
       !
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
         !WRITE(ss,'(3f9.4)')normals(nf(:,i)+1), zero
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
     END SUBROUTINE paraview

  SUBROUTINE paraview2(input,realisation,argv,nff,nelsy,g_coord,g_num,nf,nels,nod,nn,nlen,  &
              diag,ddylds,d1x1,d2x1,gravlo,loads,fcont,mv,kdiag,vcm,residual,damp,KVC,nels_2D,penpos,bod)

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::g_coord(:,:),mv(:)
      REAL(iwp),INTENT(IN)::diag(:),loads(:),d1x1(:),d2x1(:),fcont(:),vcm(:),residual(:),damp(:),KVC(:)
      REAL(iwp),INTENT(IN)::ddylds(:),gravlo(:)
      INTEGER,INTENT(IN)::input,nels,nod,nn,nlen,g_num(:,:),nf(:,:),realisation,kdiag(:),nff,nelsy,nels_2D,penpos(:),bod
      CHARACTER(*),INTENT(IN)::argv
      INTEGER::i,iel,ss,m,n
      REAL(iwp):: zero=0.0_iwp,dis_vec(2,nn),x
      character (len=8) ::  cnumber,cnumber1
      ss=input!+100
      write(cnumber,'(i8.6)') ss
      write(cnumber1,'(i8.6)') realisation
      ss=15
       !open(ss,FILE = argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
      IF(bod==1)THEN
          OPEN(ss,FILE="Output/Paraview2_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
      ELSE IF(bod==2)THEN
          OPEN(ss,FILE="Output/Paraview_2DM_1/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
      ELSE
          OPEN(ss,FILE="Output/Paraview_2DM_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
      END IF
      
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

       CASE(8)
         WRITE(ss,'(1A5,2I5)')"CELLS", nels, nels*(1+8)
       DO iel=1, nels
         WRITE(ss,'(9I5)')nod,  &
          g_num(3,iel)-1, g_num(5,iel)-1, g_num(7,iel)-1, g_num(1,iel)-1, &
                        g_num(4,iel)-1, g_num(6,iel)-1, g_num(8,iel)-1, g_num(2,iel)-1
       END DO

       WRITE(ss,'(a)')""
       WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
       DO iel = 1 , nels
         WRITE(ss,*)"23"
       END DO
       
      CASE(9)
         WRITE(ss,'(1A5,2I5)')"CELLS", (nels), (nels)*(1+9)
       DO iel=1, nels
         WRITE(ss,'(10I5)')nod,  &
          g_num(3,iel)-1, g_num(5,iel)-1, g_num(7,iel)-1, g_num(1,iel)-1, &
                        g_num(4,iel)-1, g_num(6,iel)-1, g_num(8,iel)-1, g_num(2,iel)-1, &
            g_num(9,iel)-1 
       END DO

       WRITE(ss,'(a)')""
       WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
       DO iel = 1 , nels
         WRITE(ss,*)"28"
         END DO

       CASE DEFAULT
          WRITE(*,*)"wrong number of nodes input in paraview"
      END SELECT

       WRITE(ss,'(a)')""
       !
       WRITE(ss,'(1A10,1I5)')"POINT_DATA", nn

        WRITE(ss,'(a)')"vectors Ftotal float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')loads(nf(1:2,i)+1), zero
        END DO

        WRITE(ss,'(a)')"vectors Fint_s float "
        DO  i=1, nn
         WRITE(ss,'(3f15.6)')ddylds(nf(1:2,i)+1), zero
        END DO
    

        WRITE(ss,'(a)')"vectors Fext_s float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')gravlo(nf(1:2,i)+1), zero
        END DO
        
      !  
      WRITE(ss,'(a)')"vectors fkin_s float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')vcm(nf(1:2,i)+1), zero
        END DO
        
              WRITE(ss,'(a)')"vectors f_damp float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')damp(nf(1:2,i)+1), zero
        END DO
        
        
         WRITE(ss,'(a)')"vectors R_s float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')residual(nf(1:2,i)+1), zero
        END DO
        
        
       !WRITE(ss,'(a)')"vectors normals float "
       ! DO i=1, nn
       !   WRITE(ss,'(3f15.6)')normals(1,i),normals(2,i), zero
       !  !WRITE(ss,'(3f9.4)')normals(nf(:,i)+1), zero
       ! END DO 
        
       WRITE(ss,'(a)')"vectors Fcont float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')fcont(nf(1:2,i)+1), zero
        END DO

       WRITE(ss,'(a)')"vectors Vel_s float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')d1x1(nf(1:2,i)+1), zero
        END DO
        

        WRITE(ss,'(a)')"vectors Acc_s float "
        DO i=1, nn
         WRITE(ss,'(3f15.6)')d2x1(nf(1:2,i)+1), zero
        END DO 
        

        WRITE(ss,'(a)')"vectors Mv_s float "
        DO i=1, nn
        IF(nf(1,i)==0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')zero,zero, zero
         IF(nf(1,i)>0.and.nf(2,i)==0)WRITE(ss,'(3f15.4)')mv(kdiag(nf(1,i))),zero, zero
         IF(nf(1,i)==0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')zero,mv(kdiag(nf(2,i))), zero
         IF(nf(1,i)>0.and.nf(2,i)>0)WRITE(ss,'(3f15.4)')mv(kdiag(nf(1,i))),mv(kdiag(nf(2,i))), zero
        END DO
        
        !WRITE(ss,'(a)')"vectors K_coupled float "
        !DO i=1, nn
        !   WRITE(ss,'(3f15.4)')KVC(penpos(i)),KVC(penpos(i)), zero
        !END DO
        
        

       close(ss)
       RETURN
                END SUBROUTINE paraview2
                
   SUBROUTINE point_viz_1D(input,realisation,argv,gm_coord,m_stress_efe,m_stress,a_ins,   &
							m_velocity,acc,nmps,nlen,void,young,bod)
   !---- SUBROUTINE used to save visualise outputs to Paraview format ---------------------

    IMPLICIT NONE  ! sr,pore,estress, removed from subroutine since not used now
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::gm_coord(:,:),m_stress(:,:),a_ins(:,:),m_velocity(:,:),   &
                          acc(:,:),m_stress_efe(:,:),void(:),young(:)
    !REAL(iwp),INTENT(IN)::sr(:),pore(:),estress(:)
    INTEGER,INTENT(IN)::input,nmps,nlen,bod
    integer,intent(in)::realisation
    CHARACTER(*),INTENT(IN)::argv
    INTEGER::i,j,ss
    REAL(iwp):: zero=0.0_iwp
    character (len=8) ::  cnumber,cnumber1
    ! fflag : flag of a point yielded or nor
    ! realisation: used for monte-carlo simulation
    ! adgama: incremental shear strain
    ! acstrain: accumulated shear strain
    ! sr: saturation ratio
    ! pore: pore pressure
    ! estress: elastic stress?!
     ss=input
     !ss=input/20
     !/100000
     write(cnumber,'(i8.6)') ss
     write(cnumber1,'(i8.6)') realisation
     !open(ss,FILE=argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp')
     ss=15
     IF(bod==2)OPEN(ss,FILE="Paraview_1D/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber
     IF(bod==3)OPEN(ss,FILE="Paraview_1D_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber
     write(ss,*) '#'
     write(ss,*) '# Simple AVS UCD File'
     write(ss,*) '#'
     write(ss,'(2i7.1,a)') nmps, nmps, ' 3 0 0' ! Number of nodes, number of cells, number of data items per node, number of data items per cell, model number (always 0)

    do j=1,nmps ! coordinates
    write(ss,'(i5.1,3g13.5)') j, gm_coord(:,j), zero
    ! MP number, x-coordinate y coordinate zero
    end do

    do j=1,nmps ! points
    write(ss,'(2i5.1,a,i5.1)') j, 0, '1 pt ', j
    end do

    write(ss,*) '7 2 4 4 2 2 1 1'              !1 1 1  for parameters not included (see below)
   ! 1st nr = number of variables you want to display, afterwards - how many components for each
    write(ss,*) 'displacement, meters'
    write(ss,*) 'Effective-Stress, kpa'
    write(ss,*) 'Total-Stress, kpa'
    write(ss,*) 'Acc, kpa'
    write(ss,*) 'velocity, m/s'
    write(ss,*) 'Void, n/a'
    write(ss,*) 'Young, m/s'
   !write(ss,*) 'pore, kpa'
   !write(ss,*) 'strain invariant, no-unit'
   !write(ss,*) 'sr, no-unit'

   ! 1 yield; 0 nor yield
    do j=1,nmps ! properties
    write(ss,'(i5.1,16f13.5)') j, a_ins(:,j),m_stress_efe(:,j),m_stress(:,j),acc(:,j), &
								m_velocity(:,j),void(j),young(j)
    end do

    CLOSE(ss)
    RETURN
  END SUBROUTINE point_viz_1D
                
    SUBROUTINE paraview3(input,realisation,argv,g_coord,g_num,nf,nels,nod,nnf,nn,nlen,  &
              x1)

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN)::g_coord(:,:)
      REAL(iwp),INTENT(IN)::x1(:)
      INTEGER,INTENT(IN)::input,nels,nod,nn,nnf,nlen,g_num(:,:),nf(:,:),realisation
      CHARACTER(*),INTENT(IN)::argv
      INTEGER::i,iel,ss
      REAL(iwp):: zero=0.0_iwp,dis_vec(2,nn)
      character (len=8) ::  cnumber,cnumber1
      ss=input!+100
      write(cnumber,'(i8.6)') ss
      write(cnumber1,'(i8.6)') realisation
      ss=15
       !open(ss,FILE = argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
       OPEN(ss,FILE="Paraview2_3/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.vtk')
       WRITE(ss,'(a)')'# vtk DataFile Version 3.0'
       WRITE(ss,'(a)')"vtk output"
       WRITE(ss,'(a)')"ASCII"
       WRITE(ss,'(a)')""
       WRITE(ss,'(a)')"DATASET UNSTRUCTURED_GRID"

       WRITE(ss,'(1A6,1I5,1A7)')"POINTS", nnf , "float"
       DO i=1, nnf
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

       case(8)
         WRITE(ss,'(1A5,2I5)')"CELLS", nels, nels*(1+8)
       DO iel=1, nels
         WRITE(ss,'(9I5)')nod,  &
          g_num(3,iel)-1, g_num(5,iel)-1, g_num(7,iel)-1, g_num(1,iel)-1, &
                        g_num(4,iel)-1, g_num(6,iel)-1, g_num(8,iel)-1, g_num(2,iel)-1
       END DO

       WRITE(ss,'(a)')""
       WRITE(ss,'(1A10,1I5)')"CELL_TYPES", nels
       DO iel = 1 , nels
         WRITE(ss,*)"23"
         END DO

       CASE DEFAULT
          WRITE(*,*)"wrong number of nodes input in paraview"
       end select

       WRITE(ss,'(a)')""
       !
       WRITE(ss,'(1A10,1I5)')"POINT_DATA", nnf
     !!
     !  WRITE(ss,'(a)')"vectors displacement float "
     !  DO i=1,nn
     !  dis_vec(:,i)=dis(nf(:,i))
     !  WRITE(ss,'(3f9.4)') dis_vec(:,i), zero
     !  END DO
       !
       
        WRITE(ss,'(a)')"vectors Pore_p float "
        DO i=1, nn
         IF(nf(3,i)>0)WRITE(ss,'(3f15.6)')x1(nf(3,i)),x1(nf(3,i)), zero
        END DO 
        
              

       close(ss)
       RETURN
     END SUBROUTINE paraview3

     SUBROUTINE point_viz(input,realisation,argv,gm_coord,m_stress,m_stress_ini,evpt,a_ins,Devstress,   &
							mstress,mpyield,cohesion,m_velocity,acc,nmps,nlen,eps_s,eps_f)
   !---- SUBROUTINE used to save visualise outputs to Paraview format ---------------------

    IMPLICIT NONE  ! sr,pore,estress, removed from subroutine since not used now
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::gm_coord(:,:),m_stress(:,:),a_ins(:,:),evpt(:),cohesion(:),m_velocity(:,:),   &
                          Devstress(:),acc(:,:),m_stress_ini(:,:),mstress(:),mpyield(:),eps_s(:),eps_f(:)
    !REAL(iwp),INTENT(IN)::sr(:),pore(:),estress(:)
    INTEGER,INTENT(IN)::input,nmps,nlen
    integer,intent(in)::realisation
    CHARACTER(*),INTENT(IN)::argv
    INTEGER::i,j,ss
    REAL(iwp):: zero=0.0_iwp
    character (len=6) ::  cnumber,cnumber1
    ! fflag : flag of a point yielded or nor
    ! realisation: used for monte-carlo simulation
    ! adgama: incremental shear strain
    ! acstrain: accumulated shear strain
    ! sr: saturation ratio
    ! pore: pore pressure
    ! estress: elastic stress?!
     ss=input
     !ss=input/20
     !/100000
     write(cnumber,'(i6.6)') ss
     write(cnumber1,'(i6.6)') realisation
     !open(ss,FILE=argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp')
     ss=15
     OPEN(ss,FILE="Paraview/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber
     write(ss,*) '#'
     write(ss,*) '# Simple AVS UCD File'
     write(ss,*) '#'
     write(ss,'(2i7.1,a)') nmps, nmps, ' 3 0 0' ! Number of nodes, number of cells, number of data items per node, number of data items per cell, model number (always 0)

    do j=1,nmps ! coordinates
    write(ss,'(i5.1,3g13.5)') j, gm_coord(:,j), zero
    ! MP number, x-coordinate y coordinate zero
    end do

    do j=1,nmps ! points
    write(ss,'(2i5.1,a,i5.1)') j, 0, '1 pt ', j
    end do

    write(ss,*) '10 2 1 4 4 2 2 1 1 1 1'              !1 1 1  for parameters not included (see below)
   ! 1st nr = number of variables you want to display, afterwards - how many components for each
    write(ss,*) 'displacement, meters'
    write(ss,*) 'plastic strains, meters'
    write(ss,*) 'Stress, kpa'
    write(ss,*) 'Strain_s, kpa'
    write(ss,*) 'Strain_f, kpa'
    write(ss,*) 'Acc, m/s2'
    write(ss,*) 'velocity, m/s'
    write(ss,*) 'cohesion, kpa'
   !write(ss,*) 'pore, kpa'
   !write(ss,*) 'strain invariant, no-unit'
   !write(ss,*) 'sr, no-unit'
   write(ss,*) 'Dev-stress, no-unit'
   write(ss,*) 'Mean-stress, no-unit'
   write(ss,*) 'Yield-value, no-unit'
   ! 1 yield; 0 nor yield
    do j=1,nmps ! properties
    write(ss,'(i5.1,17f13.3)') j, a_ins(:,j),evpt(j),m_stress(:,j),eps_s(j),eps_f(j),acc(:,j),    &
								m_velocity(:,j),cohesion(j),Devstress(j),mstress(j),mpyield(j)
    end do

    CLOSE(ss)
    RETURN
          END SUBROUTINE point_viz

 SUBROUTINE point_viz2(input,realisation,argv,gm_coord,m_stress,m_eps_acum,evpt,a_ins,Devstress,   &
							mstress,mpyield,cohesion,m_velocity,acc,nmps,nlen,bod)
   !---- SUBROUTINE used to save visualise outputs to Paraview format ---------------------

    IMPLICIT NONE  ! sr,pore,estress, removed from subroutine since not used now
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::gm_coord(:,:),m_stress(:,:),a_ins(:,:),evpt(:),cohesion(:),m_velocity(:,:),   &
                          Devstress(:),acc(:,:),m_eps_acum(:,:),mstress(:),mpyield(:)
    INTEGER,INTENT(IN)::input,nmps,nlen,bod
    integer,intent(in)::realisation
    CHARACTER(*),INTENT(IN)::argv
    INTEGER::i,j,ss
    REAL(iwp):: zero=0.0_iwp
    character (len=8) ::  cnumber,cnumber1
    ! fflag : flag of a point yielded or nor
    ! realisation: used for monte-carlo simulation
    ! adgama: incremental shear strain
    ! acstrain: accumulated shear strain
    ! sr: saturation ratio
    ! pore: pore pressure
    ! estress: elastic stress?!
     ss=input
     !ss=input/20
     !/100000
     write(cnumber,'(i8.6)') ss
     write(cnumber1,'(i8.6)') realisation
     !open(ss,FILE=argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp')
     ss=15
     IF(bod==1)THEN
      OPEN(ss,FILE="Output/Paraview2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber
     ELSE IF(bod==2)THEN
      OPEN(ss,FILE="Output/Paraview_2DP_1/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber   
     ELSE
      OPEN(ss,FILE="Output/Paraview_2DP_2/"//argv(1:nlen)//"_"//trim(adjustl(cnumber1))//"_"//trim(adjustl(cnumber))//'.inp') ! Creates a file with title argv_cnumber1_cnumber   
     END IF
     
     write(ss,*) '#'
     write(ss,*) '# Simple AVS UCD File'
     write(ss,*) '#'
     write(ss,'(2i7.1,a)') nmps, nmps, ' 3 0 0' ! Number of nodes, number of cells, number of data items per node, number of data items per cell, model number (always 0)

    do j=1,nmps ! coordinates
    write(ss,'(i5.1,3g13.5)') j, gm_coord(:,j), zero
    ! MP number, x-coordinate y coordinate zero
    end do

    do j=1,nmps ! points
    write(ss,'(2i5.1,a,i5.1)') j, 0, '1 pt ', j
    end do

    write(ss,*) '10 2 1 4 4 2 2 1 1 1 1'              !1 1 1  for parameters not included (see below)
   ! 1st nr = number of variables you want to display, afterwards - how many components for each
    write(ss,*) 'displacement, meters'
    write(ss,*) 'Void, meters'
    write(ss,*) 'Effective-Stress, kpa'
    write(ss,*) 'Total-strains, no-unit'
    write(ss,*) 'Acc, kpa'
    write(ss,*) 'velocity, m/s'
    write(ss,*) 'Young, kpa'
   !write(ss,*) 'pore, kpa'
   !write(ss,*) 'strain invariant, no-unit'
   !write(ss,*) 'sr, no-unit'
   write(ss,*) 'Dev-stress, no-unit'
   write(ss,*) 'Mean-stress, no-unit'
   write(ss,*) 'Pore-pressure, no-unit'

   ! 1 yield; 0 nor yield
    do j=1,nmps ! properties
    write(ss,'(i5.1,20f15.4)') j, a_ins(:,j),evpt(j),m_stress(:,j),m_eps_acum(:,j),acc(:,j), &
								m_velocity(:,j),cohesion(j),Devstress(j),mstress(j),mpyield(j)
    end do

    CLOSE(ss)
    RETURN
  END SUBROUTINE point_viz2

      SUBROUTINE couma(nels,nmps,a_ele,c_ele,k_ele,etype)

       ! count the number of material points inside an element.

       IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::a_ele(:)
      INTEGER,INTENT(IN)::nels,nmps
      INTEGER,INTENT(OUT)::c_ele(:),k_ele(0:),etype(:)
      INTEGER::i,iel,count
        k_ele=0
       DO iel=1,nels
         count=0
         DO i=1,nmps
             IF(iel==a_ele(i)) THEN
                 count=count+1
             END IF
         END DO
         c_ele(iel)=count
         k_ele(iel)=k_ele(iel-1)+c_ele(iel)
         IF(count==0) THEN
             etype(iel)=1
         ELSE
             etype(iel)=2
         END IF
       END DO
         RETURN
      END SUBROUTINE COUMA

      SUBROUTINE point_search(coord,gm_coord,i,converged)

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::i
      REAL(iwp),INTENT(IN)::coord(:,:),gm_coord(:,:)
      LOGICAL,INTENT(OUT)::converged
      REAL(iwp)::xa,ya,x1,y1,x2,y2,x3,y3,x4,y4
      REAL(iwp)::a1,a2,a3,a4
      REAL,PARAMETER::zero=0.0_iwp
      INTEGER::ndim,nod
      nod= UBOUND(coord,1)
      ndim=UBOUND(coord,2)
      CONVERGED=.FALSE.
      SELECT CASE(ndim)
      CASE(2)
      xa=gm_coord(1,i)
      ya=gm_coord(2,i)
      SELECT CASE(nod)
      CASE(4)
      x1=coord(1,1)
      y1=coord(1,2)
      x2=coord(4,1)
      y2=coord(4,2)
      x3=coord(3,1)
      y3=coord(3,2)
      x4=coord(2,1)
      y4=coord(2,2)
      a1=(x2-x1)*(ya-y1)-(xa-x1)*(y2-y1)
      a2=(x3-x2)*(ya-y2)-(xa-x2)*(y3-y2)
      a3=(x4-x3)*(ya-y3)-(xa-x3)*(y4-y3)
      a4=(x1-x4)*(ya-y4)-(xa-x4)*(y1-y4)
      IF(a1>zero.AND.a2>zero.AND.a3>zero.AND.a4>zero) THEN   ! my guess: feel free to change if you think a more appropriate algorithm works
     ! IF(a1>zero.AND.a2>zero.AND.a3>zero.AND.a4>zero.or.a1<zero.AND.a2<zero.AND.a3<zero.AND.a4<zero) THEN
      CONVERGED=.TRUE.
      END IF
      CASE(8)
      x1=coord(1,1)
      y1=coord(1,2)
      x2=coord(7,1)
      y2=coord(7,2)
      x3=coord(5,1)
      y3=coord(5,2)
      x4=coord(3,1)
      y4=coord(3,2)
      a1=(x2-x1)*(ya-y1)-(xa-x1)*(y2-y1)
      a2=(x3-x2)*(ya-y2)-(xa-x2)*(y3-y2)
      a3=(x4-x3)*(ya-y3)-(xa-x3)*(y4-y3)
      a4=(x1-x4)*(ya-y4)-(xa-x4)*(y1-y4)
      IF(a1>zero.AND.a2>zero.AND.a3>zero.AND.a4>zero) THEN
      CONVERGED=.TRUE.
      END IF
       CASE DEFAULT
        WRITE(*,*)"wrong number of nodes in point_search"
        END SELECT

        CASE DEFAULT
        WRITE(*,*)"wrong number of dimensions in point_search"
      END SELECT
      RETURN
      END SUBROUTINE point_search

      SUBROUTINE sort(vector)
       IMPLICIT NONE
       INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
       INTEGER::N,numbery,i,fil
       REAL(iwp)::zero=0.0_iwp,x
       INTEGER,INTENT(INOUT)::vector(:)
       INTEGER,ALLOCATABLE::vectoraux(:)
        fil=1
        N=zero
        numbery=ubound(vector,1)
        ALLOCATE(vectoraux(numbery))
        vectoraux=vector
        i=0

        Sorting:DO
         i=i+1
         IF(numbery>i)THEN
         IF(vector(i+1)<vector(i))THEN
           N=vector(i)
           vector(i)=vector(i+1)
           vector(i+1)=N
           i=0
         END IF
         END IF
         IF(i==numbery-1)EXIT
         IF(numbery==1)EXIT
        END DO Sorting

       RETURN
      END SUBROUTINE sort


          SUBROUTINE count2D(matrix,N)
           IMPLICIT NONE
           INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
           INTEGER::N,numberx,numbery,i,fil,col,total
           REAL(iwp)::zero=0.0_iwp
           INTEGER,INTENT(IN)::matrix(:,:)
            fil=1
            col=1
            N=zero
            numbery=ubound(matrix,1)
            numberx=ubound(matrix,2)
            total=numberx*numbery
            i=1
            DO WHILE (i<=total)
            IF(matrix(fil,col)>0.0)N=N+1
            IF(fil>=numbery)THEN
              fil=0
              col=col+1
            END IF
            fil=fil+1
            i=i+1
            END DO
           RETURN
          END SUBROUTINE count2D

      SUBROUTINE softening(ci,cp,cr,hd,sp,theta,phi,gamap1,gamap2,model_formulation)
       !--------------------------------------------------------------------------
       ! works out the plastic strains for YAP's cohesion-softening model
       !--------------------------------------------------------------------------
       implicit none
       INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
       integer, intent(in):: model_formulation
       real(iwp),intent(in):: ci,cp,cr,hd,sp,theta,phi
       real(iwp),intent(out):: gamap1,gamap2
      !local variables
       real(iwp):: sq3,csph,csth,snth,dprime
         lode_angle: select case (model_formulation)
            case (1)    ! YAP original formulation
               gamap1=(cp-ci)/hd                         !hd=hardening modulus
               gamap2=(cr-cp)/sp+gamap1                  !sp=softening parameter
            case default
               print*, "check subroutine 'SOFTENING' "
         end select lode_angle
       RETURN
      END SUBROUTINE softening

      !SUBROUTINE mocouf(phi,c,sigm,dsbar,theta,f)
      !!x
      !! This subroutine calculates the value of the yield function
      !! for a Mohr-Coulomb material (phi in degrees).
      !!
      ! IMPLICIT NONE
      ! INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      ! REAL(iwp),INTENT(IN)::phi,c,sigm,dsbar,theta
      ! REAL(iwp),INTENT(OUT)::f
      ! REAL(iwp)::phir,snph,csph,csth,snth,one=1.0_iwp,d3=3.0_iwp,d4=4.0_iwp,   &
      !   d180=180.0_iwp
      ! phir=phi*d4*ATAN(one)/d180
      ! snph=SIN(phir)
      ! csph=COS(phir)
      ! csth=COS(theta)
      ! snth=SIN(theta)
      ! f=snph*sigm+dsbar*(csth/SQRT(d3)-snth*snph/d3)-c*csph
      !RETURN
      !END SUBROUTINE mocouf


    SUBROUTINE bound_detect(iel,neighb,bound,values,c_ele,nf,g_num,aver)
     !
     !   This subroutine produces derivatives of shape functions withe respect
     !   to local coordinates for 16 nodes.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(OUT)::bound,values,aver
      INTEGER,INTENT(IN)::neighb(:,:),iel,c_ele(:),g_num(:,:),nf(:,:)
      INTEGER::i,CONT
      REAL,ALLOCATABLE::NEIG(:)
      ALLOCATE(NEIG(8))
      CONT=0
      NEIG=0
      bound=0

      !Next loop is to know if there are neighbour elements and if those elements have MPs inside
      DO i=1,8
        IF(neighb(iel,i)==0)THEN
          NEIG(i)=0
        ELSE
          IF(c_ele(neighb(iel,i))>=1)NEIG(i)=1
          IF(c_ele(neighb(iel,i))==0)NEIG(i)=0
        END IF
      END DO
     !DO i=1,8
     ! IF(neighb(iel,i)==0)THEN
     !      NEIG(i)=0
     ! ELSE
     !   IF(((nf(1,g_num(1,neighb(iel,i)))>0.or.nf(2,g_num(1,neighb(iel,i)))>0).and.  &
     !       (nf(1,g_num(2,neighb(iel,i)))>0.or.nf(2,g_num(2,neighb(iel,i)))>0).and.  &
     !       (nf(1,g_num(3,neighb(iel,i)))>0.or.nf(2,g_num(3,neighb(iel,i)))>0).and.  &
     !       (nf(1,g_num(4,neighb(iel,i)))>0.or.nf(2,g_num(4,neighb(iel,i)))>0)).or.  &
     !        c_ele(neighb(iel,i))>0)THEN
     !
     !           NEIG(i)=1
     !           ELSE
     !           NEIG(i)=0
     !    END IF
     !  END IF 
     !END DO

      IF(NEIG(1)==0.and.NEIG(2)==0.and.NEIG(3)==0.and.NEIG(4)>=1             &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)==0.and.NEIG(8)==0)THEN
        bound=1;aver=1 !upper left corner
      END IF
      IF(NEIG(1)==0.and.NEIG(2)==0.and.NEIG(3)==0.and.NEIG(4)>=1             &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=2;aver=2 !upper boundarie
      END IF
      IF(NEIG(1)==0.and.NEIG(2)==0.and.NEIG(3)==0.and.NEIG(4)==0             &
        .and.NEIG(5)==0.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=3;aver=1 !upper right
      END IF
      IF(NEIG(1)==0.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1             &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)==0.and.NEIG(8)==0)THEN
        bound=4;aver=2 !left boundarie
      END IF
      IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1             &
        .and.NEIG(5)>=1.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=5;aver=1 !center element
      END IF
      IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)==0.and.NEIG(4)==0             &
        .and.NEIG(5)==0.and.NEIG(6)>=1.and.NEIG(7)>=1.and.NEIG(8)>=1)THEN
        bound=6;aver=2 !right boundarie
      END IF
      IF(NEIG(1)==0.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1             &
        .and.NEIG(5)==0.and.NEIG(6)==0.and.NEIG(7)==0.and.NEIG(8)==0)THEN
        bound=7;aver=1 !lower left
      END IF
      IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)>=1.and.NEIG(4)>=1             &
        .and.NEIG(5)==0.and.NEIG(6)==0.and.NEIG(7)==0.and.NEIG(8)>=1)THEN
        bound=8;aver=2 !lower boundarie
      END IF
      IF(NEIG(1)>=1.and.NEIG(2)>=1.and.NEIG(3)==0.and.NEIG(4)==0             &
        .and.NEIG(5)==0.and.NEIG(6)==0.and.NEIG(7)==0.and.NEIG(8)>=1)THEN
        bound=9;aver=1 !lower right
      END IF
      IF(bound>=1.and.bound/=5)THEN
        values=9
      ELSE IF(bound==5)THEN
        values=16
      ELSE
        values=4
      END IF

     RETURN
     END SUBROUTINE bound_detect

     SUBROUTINE neigh_b(nnxe,nnye,nels,neighb,dir)

   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   INTEGER,INTENT(IN)::nnxe,nnye,nels
   INTEGER,INTENT(OUT)::neighb(:,:)
   CHARACTER(LEN=1),INTENT(IN)::dir
   INTEGER::row,col,i,neighbor

   row=1
   col=1
   IF(dir=='y')THEN
    DO i=1, nels
      IF(row==1.and.col==1)THEN
        neighb(i,1)=0
        neighb(i,2)=0
        neighb(i,3)=0
        neighb(i,7)=0
        neighb(i,8)=0
        neighb(i,4)=nnye*col+row
        neighb(i,5)=nnye*col+row+1
        neighb(i,6)=nnye*(col-1)+row+1
      ELSE IF(col==1.and.row>1.and.row<nnye)THEN
        neighb(i,1)=0
        neighb(i,7)=0
        neighb(i,8)=0
        neighb(i,2)=nnye*(col-1)+row-1
        neighb(i,3)=nnye*col+row-1
        neighb(i,4)=nnye*col+row
        neighb(i,5)=nnye*col+row+1
        neighb(i,6)=nnye*(col-1)+row+1
      ELSE IF(row==nnye.and.col==1)THEN
        neighb(i,1)=0
        neighb(i,5)=0
        neighb(i,6)=0
        neighb(i,7)=0
        neighb(i,8)=0
        neighb(i,2)=nnye*(col-1)+row-1
        neighb(i,3)=nnye*col+row-1
        neighb(i,4)=nnye*col+row
      ELSE IF(row==nnye.and.col>1.and.col<nnxe)THEN
        neighb(i,5)=0
        neighb(i,6)=0
        neighb(i,7)=0
        neighb(i,1)=nnye*(col-2)+row-1
        neighb(i,2)=nnye*(col-1)+row-1
        neighb(i,3)=nnye*col+row-1
        neighb(i,4)=nnye*col+row
        neighb(i,8)=nnye*(col-2)+row
      ELSE IF(row==nnye.and.col==nnxe)THEN
        neighb(i,3)=0
        neighb(i,4)=0
        neighb(i,5)=0
        neighb(i,6)=0
        neighb(i,7)=0
        neighb(i,1)=nnye*(col-2)+row-1
        neighb(i,2)=nnye*(col-1)+row-1
        neighb(i,8)=nnye*(col-2)+row
      ELSE IF(col==nnxe.and.row>1.and.row<nnye)THEN
        neighb(i,3)=0
        neighb(i,4)=0
        neighb(i,5)=0
        neighb(i,1)=nnye*(col-2)+row-1
        neighb(i,2)=nnye*(col-1)+row-1
        neighb(i,6)=nnye*(col-1)+row+1
        neighb(i,7)=nnye*(col-2)+row+1
        neighb(i,8)=nnye*(col-2)+row
      ELSE IF(row==1.and.col==nnxe)THEN
        neighb(i,1)=0
        neighb(i,2)=0
        neighb(i,3)=0
        neighb(i,4)=0
        neighb(i,5)=0
        neighb(i,6)=nnye*(col-1)+row+1
        neighb(i,7)=nnye*(col-2)+row+1
        neighb(i,8)=nnye*(col-2)+row
      ELSE IF(row==1.and.col>1.and.col<nnxe)THEN
        neighb(i,1)=0
        neighb(i,2)=0
        neighb(i,3)=0
        neighb(i,4)=nnye*col+row
        neighb(i,5)=nnye*col+row+1
        neighb(i,6)=nnye*(col-1)+row+1
        neighb(i,7)=nnye*(col-2)+row+1
        neighb(i,8)=nnye*(col-2)+row
      ELSE IF(row>1.and.row<nnye.and.col>1.and.col<nnxe)THEN
        neighb(i,1)=nnye*(col-2)+row-1
        neighb(i,2)=nnye*(col-1)+row-1
        neighb(i,3)=nnye*col+row-1
        neighb(i,4)=nnye*col+row
        neighb(i,5)=nnye*col+row+1
        neighb(i,6)=nnye*(col-1)+row+1
        neighb(i,7)=nnye*(col-2)+row+1
        neighb(i,8)=nnye*(col-2)+row
      END IF

      IF(row==nnye)THEN
        row=0
        col=col+1
      END IF
      row=row+1

    END DO
   ELSE
     DO i=1, nels
      IF(row==1.and.col==1)THEN
        neighb(i,1)=0
        neighb(i,2)=0
        neighb(i,3)=0
        neighb(i,7)=0
        neighb(i,8)=0
        neighb(i,4)=col+1
        neighb(i,5)=nnxe*row+col+1
        neighb(i,6)=nnxe*row+1
      ELSE IF(col==1.and.row>1.and.row<nnye)THEN
        neighb(i,1)=0
        neighb(i,7)=0
        neighb(i,8)=0
        neighb(i,2)=nnxe*(row-2)+1
        neighb(i,3)=nnxe*(row-2)+2
        neighb(i,4)=nnxe*(row-1)+2
        neighb(i,5)=nnxe*row+2
        neighb(i,6)=nnxe*row+1
      ELSE IF(row==nnye.and.col==1)THEN
        neighb(i,1)=0
        neighb(i,5)=0
        neighb(i,6)=0
        neighb(i,7)=0
        neighb(i,8)=0
        neighb(i,2)=nnxe*(row-2)+1
        neighb(i,3)=nnxe*(row-2)+2
        neighb(i,4)=nnxe*(row-1)+2
      ELSE IF(row==nnye.and.col>1.and.col<nnxe)THEN
        neighb(i,5)=0
        neighb(i,6)=0
        neighb(i,7)=0
        neighb(i,1)=nnxe*(row-2)+col-1
        neighb(i,2)=nnxe*(row-2)+col
        neighb(i,3)=nnxe*(row-2)+col+1
        neighb(i,4)=nnxe*(row-1)+col+1
        neighb(i,8)=nnxe*(row-1)+col-1
      ELSE IF(row==nnye.and.col==nnxe)THEN
        neighb(i,3)=0
        neighb(i,4)=0
        neighb(i,5)=0
        neighb(i,6)=0
        neighb(i,7)=0
        neighb(i,1)=nnxe*(row-2)+col-1
        neighb(i,2)=nnxe*(row-2)+col
        neighb(i,8)=nnxe*(row-1)+col-1
      ELSE IF(col==nnxe.and.row>1.and.row<nnye)THEN
        neighb(i,3)=0
        neighb(i,4)=0
        neighb(i,5)=0
        neighb(i,1)=nnxe*(row-2)+col-1
        neighb(i,2)=nnxe*(row-2)+col
        neighb(i,6)=nnxe*row+col
        neighb(i,7)=nnxe*row+col-1
        neighb(i,8)=nnxe*(row-1)+col-1
      ELSE IF(row==1.and.col==nnxe)THEN
        neighb(i,1)=0
        neighb(i,2)=0
        neighb(i,3)=0
        neighb(i,4)=0
        neighb(i,5)=0
        neighb(i,6)=nnxe*row+col
        neighb(i,7)=nnxe*row+col-1
        neighb(i,8)=nnxe*(row-1)+col-1
      ELSE IF(row==1.and.col>1.and.col<nnxe)THEN
        neighb(i,1)=0
        neighb(i,2)=0
        neighb(i,3)=0
        neighb(i,4)=col+1
        neighb(i,5)=nnxe*row+col+1
        neighb(i,6)=nnxe*row+col
        neighb(i,7)=nnxe*row+col-1
        neighb(i,8)=col-1
      ELSE IF(row>1.and.row<nnye.and.col>1.and.col<nnxe)THEN
        neighb(i,1)=nnxe*(row-2)+col-1
        neighb(i,2)=nnxe*(row-2)+col
        neighb(i,3)=nnxe*(row-2)+col+1
        neighb(i,4)=nnxe*(row-1)+col+1
        neighb(i,5)=nnxe*(row)+col+1
        neighb(i,6)=nnxe*(row)+col
        neighb(i,7)=nnxe*(row)+col-1
        neighb(i,8)=nnxe*(row-1)+col-1
      END IF

      IF(col==nnxe)THEN
        col=0
        row=row+1
      END IF
      col=col+1
    END DO
   END IF
 RETURN
 END SUBROUTINE neigh_b

   SUBROUTINE const(EQUATIONS,n,bound,m)

     ! this soubrutine is to compute the constants of a 9x9 or 16x16
     ! matrix, considering also the position of the element

      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(OUT)::EQUATIONS(:,:)
      INTEGER,INTENT(IN)::n,bound,m
      INTEGER::x,y,fil,col

      equations=0.0

      IF(bound==5)THEN
        x=-3;y=3;fil=1;col=1
        DO fil=1,16
         DO col=1,16
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
           IF(col==10)EQUATIONS(fil,col)=x**3
           IF(col==11)EQUATIONS(fil,col)=y**3
           IF(col==12)EQUATIONS(fil,col)=x**3*y
           IF(col==13)EQUATIONS(fil,col)=x*y**3
           IF(col==14)EQUATIONS(fil,col)=x**3*y**2
           IF(col==15)EQUATIONS(fil,col)=x**2*y**3
           IF(col==16)EQUATIONS(fil,col)=x**3*y**3
         END DO
         x=x+2
         IF(x>3)THEN
           x=-3
           y=y-2
         END IF
       END DO
      END IF

      IF((bound==1.or.bound==4.or.bound==2).and.m==1)THEN
        x=-1;y=1;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>3)THEN
           x=-1
           y=y-2
         END IF
       END DO
      END IF

    IF(bound==4.and.m==2)THEN
        x=-1;y=3;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2

         END DO
         x=x+2
         IF(x>3)THEN
           x=-1
           y=y-2
         END IF
       END DO
      END IF

     IF(bound==2.and.m==2)THEN
        x=-3;y=1;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>1)THEN
           x=-3
           y=y-2
         END IF
       END DO
      END IF

      IF((bound==3.or.bound==6).and.m==1)THEN
        x=-3;y=1;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>1)THEN
           x=-3
           y=y-2
         END IF
       END DO
      END IF

     IF(bound==6.and.m==2)THEN
        x=-3;y=3;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>1)THEN
           x=-3
           y=y-2
         END IF
       END DO
      END IF

      IF((bound==7.or.bound==8).and.m==1)THEN
        x=-1;y=3;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>3)THEN
           x=-1
           y=y-2
         END IF
       END DO
      END IF

     IF(bound==8.and.m==2)THEN
        x=-3;y=3;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>1)THEN
           x=-3
           y=y-2
         END IF
       END DO
      END IF

      IF(bound==9)THEN
        x=-3;y=3;fil=1;col=1
        DO fil=1,9
         DO col=1,9
           IF(col==1)EQUATIONS(fil,col)=1.0
           IF(col==2)EQUATIONS(fil,col)=x
           IF(col==3)EQUATIONS(fil,col)=y
           IF(col==4)EQUATIONS(fil,col)=y*x
           IF(col==5)EQUATIONS(fil,col)=x**2
           IF(col==6)EQUATIONS(fil,col)=y**2
           IF(col==7)EQUATIONS(fil,col)=y*x**2
           IF(col==8)EQUATIONS(fil,col)=x*y**2
           IF(col==9)EQUATIONS(fil,col)=x**2*y**2
         END DO
         x=x+2
         IF(x>1)THEN
           x=-3
           y=y-2
         END IF
       END DO
      END IF
      RETURN
     END SUBROUTINE const

     SUBROUTINE inv(matrix,values)
     !
     ! This subroutine inverts the values x values matrix
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      REAL(iwp),INTENT(IN OUT)::matrix(:,:)
      INTEGER,INTENT(IN)::values
      INTEGER,ALLOCATABLE::ipvt(:)
      REAL(iwp)::c,d,zero=0.0_iwp,one=1.0_iwp
      REAL(iwp),ALLOCATABLE::INVERSE(:,:),temp(:)
      INTEGER::l,k,N,M,j,imax(1)

          ALLOCATE(INVERSE(values,values),ipvt(values),temp(values))
          N=UBOUND(matrix,1)
          INVERSE=zero
          temp=zero
          ipvt=zero
          INVERSE = matrix
          ipvt = (/ (l, l = 1, N) /)
          DO k = 1,N
             imax = MAXLOC(ABS(INVERSE(k:n,k)))
             m = k-1+imax(1)
             IF (M /= k) THEN
                ipvt( (/M,k/) ) = ipvt( (/k,M/) )
                INVERSE((/M,k/),:) = INVERSE((/k,M/),:)
             END IF
             d = one/INVERSE(k,k)
             temp = INVERSE(:,k)
             DO j = 1, N
                c = INVERSE(k,j)*d
                INVERSE(:,j) = INVERSE(:,j)-temp*c
                INVERSE(k,j) = c
             END DO
             INVERSE(:,k) = temp*(-d)
             INVERSE(k,k) = d
          END DO
          matrix(:,ipvt) = INVERSE
     RETURN
     END SUBROUTINE inv
     
       subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
integer,INTENT(IN):: n
REAL(iwp),INTENT(IN OUT):: a(n,n), c(n,n)
REAL(iwp):: L(n,n), U(n,n), b(n), d(n), x(n)
REAL(iwp):: coeff
INTEGER:: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse


     SUBROUTINE shapefunct(points,i,CONSTANTS,SFUNCT)
     !
     !   This subroutine produces derivatives of shape functions withe respect
     !   to local coordinates for 16 and 9 nodes.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::i
      REAL(iwp)::eta,xi
      REAL(iwp),INTENT(IN)::points(:,:)
      REAL(iwp),INTENT(IN)::CONSTANTS(:,:)
      REAL(iwp),INTENT(OUT)::SFUNCT(:)
      REAL,PARAMETER::zero=0.0_iwp
      INTEGER::ndim,n
      xi=points(i,1)
      eta=points(i,2)
      ndim=UBOUND(CONSTANTS,1)

      IF(ndim==16)THEN
       DO n=1,ndim
        SFUNCT(n)=CONSTANTS(1,n)+CONSTANTS(2,n)*xi+eta*CONSTANTS(3,n)+xi*eta*CONSTANTS(4,n)+ &
                   CONSTANTS(5,n)*xi**2+CONSTANTS(6,n)*eta**2+CONSTANTS(7,n)*xi**2*eta+ &
                   CONSTANTS(8,n)*xi*eta**2+CONSTANTS(9,n)*xi**2*eta**2+CONSTANTS(10,n)*xi**3+ &
                   CONSTANTS(11,n)*eta**3+CONSTANTS(12,n)*eta*xi**3+CONSTANTS(13,n)*eta**3*xi+ &
                   CONSTANTS(14,n)*xi**3*eta**2+CONSTANTS(15,n)*xi**2*eta**3+CONSTANTS(16,n)*xi**3*eta**3
       END DO
      ELSE

       DO n=1,ndim
        SFUNCT(n)=CONSTANTS(1,n)+CONSTANTS(2,n)*xi+eta*CONSTANTS(3,n)+xi*eta*CONSTANTS(4,n)+ &
                   CONSTANTS(5,n)*xi**2+CONSTANTS(6,n)*eta**2+CONSTANTS(7,n)*xi**2*eta+ &
                   CONSTANTS(8,n)*xi*eta**2+CONSTANTS(9,n)*xi**2*eta**2
       END DO

      END IF
     RETURN
     END SUBROUTINE shapefunct

     SUBROUTINE derivatives(points,i,CONSTANTS,SIXTEEND)
     !
     !   This subroutine produces derivatives of shape functions withe respect
     !   to local coordinates for 16 and 9 nodes.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::i
      REAL(iwp)::eta,xi
      REAL(iwp),INTENT(IN)::points(:,:)
      REAL(iwp),INTENT(IN)::CONSTANTS(:,:)
      REAL(iwp),INTENT(OUT)::SIXTEEND(:,:)
      REAL,PARAMETER::zero=0.0_iwp
      INTEGER::m,ndim,n
      xi=points(i,1)
      eta=points(i,2)
      ndim=UBOUND(CONSTANTS,1)
      IF(ndim==16)THEN
       DO n=1,ndim
        DO m=1,2
          IF(m==1)SIXTEEND(m,n)=CONSTANTS(2,n)+eta*CONSTANTS(4,n)+2*xi*CONSTANTS(5,n)+ &
             2*CONSTANTS(7,n)*xi*eta+CONSTANTS(8,n)*eta**2+2*CONSTANTS(9,n)*xi*eta**2+ &
             3*CONSTANTS(10,n)*xi**2+3*CONSTANTS(12,n)*xi**2*eta+CONSTANTS(13,n)*eta**3+ &
             3*CONSTANTS(14,n)*xi**2*eta**2+2*CONSTANTS(15,n)*eta**3*xi+ &
             3*CONSTANTS(16,n)*xi**2*eta**3
          IF(m==2)SIXTEEND(m,n)=CONSTANTS(3,n)+xi*CONSTANTS(4,n)+2*eta*CONSTANTS(6,n)+ &
             CONSTANTS(7,n)*xi**2+2*CONSTANTS(8,n)*eta*xi+2*CONSTANTS(9,n)*xi**2*eta+ &
             3*CONSTANTS(11,n)*eta**2+CONSTANTS(12,n)*xi**3+3*CONSTANTS(13,n)*eta**2*xi+ &
             2*CONSTANTS(14,n)*xi**3*eta+3*CONSTANTS(15,n)*eta**2*xi**2+ &
             3*CONSTANTS(16,n)*xi**3*eta**2

        END DO
       END DO
      END IF
      IF(ndim==9)THEN
       DO n=1,ndim
        DO m=1,2
           IF(m==1)SIXTEEND(m,n)=CONSTANTS(2,n)+eta*CONSTANTS(4,n)+2*xi*CONSTANTS(5,n)+ &
              2*CONSTANTS(7,n)*xi*eta+CONSTANTS(8,n)*eta**2+2*CONSTANTS(9,n)*xi*eta**2
           IF(m==2)SIXTEEND(m,n)=CONSTANTS(3,n)+xi*CONSTANTS(4,n)+2*eta*CONSTANTS(6,n)+ &
              CONSTANTS(7,n)*xi**2+2*CONSTANTS(8,n)*eta*xi+2*CONSTANTS(9,n)*xi**2*eta
         END DO
       END DO
      END IF
     RETURN
     END SUBROUTINE derivatives

     SUBROUTINE beemat3(bee,deriv,values,SFUNCT)
     !
     ! This subroutine forms the bee matrix for 16 and 9 nodes.
     !
      IMPLICIT NONE
      INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
      INTEGER,INTENT(IN)::values
      REAL(iwp),INTENT(IN)::SFUNCT(:),deriv(:,:)
      REAL(iwp),INTENT(OUT)::bee(:,:)
      INTEGER::l,m,n,i,nod,k

      bee=0.0_iwp

     IF(values==16)THEN
       bee(1,1)=deriv(1,1)
       bee(3,1)=deriv(2,1)
       bee(2,2)=deriv(2,1)
       bee(3,2)=deriv(1,1)
       bee(1,3)=deriv(1,2)
       bee(3,3)=deriv(2,2)
       bee(2,4)=deriv(2,2)
       bee(3,4)=deriv(1,2)
       bee(1,5)=deriv(1,3)
       bee(3,5)=deriv(2,3)
       bee(2,6)=deriv(2,3)
       bee(3,6)=deriv(1,3)
       bee(1,7)=deriv(1,4)
       bee(3,7)=deriv(2,4)
       bee(2,8)=deriv(2,4)
       bee(3,8)=deriv(1,4)
       bee(1,9)=deriv(1,5)
       bee(3,9)=deriv(2,5)
       bee(2,10)=deriv(2,5)
       bee(3,10)=deriv(1,5)
       bee(1,11)=deriv(1,6)
       bee(3,11)=deriv(2,6)
       bee(2,12)=deriv(2,6)
       bee(3,12)=deriv(1,6)
       bee(1,13)=deriv(1,7)
       bee(3,13)=deriv(2,7)
       bee(2,14)=deriv(2,7)
       bee(3,14)=deriv(1,7)
       bee(1,15)=deriv(1,8)
       bee(3,15)=deriv(2,8)
       bee(2,16)=deriv(2,8)
       bee(3,16)=deriv(1,8)
       bee(1,17)=deriv(1,9)
       bee(3,17)=deriv(2,9)
       bee(2,18)=deriv(2,9)
       bee(3,18)=deriv(1,9)
       bee(1,19)=deriv(1,10)
       bee(3,19)=deriv(2,10)
       bee(2,20)=deriv(2,10)
       bee(3,20)=deriv(1,10)
       bee(1,21)=deriv(1,11)
       bee(3,21)=deriv(2,11)
       bee(2,22)=deriv(2,11)
       bee(3,22)=deriv(1,11)
       bee(1,23)=deriv(1,12)
       bee(3,23)=deriv(2,12)
       bee(2,24)=deriv(2,12)
       bee(3,24)=deriv(1,12)
       bee(1,25)=deriv(1,13)
       bee(3,25)=deriv(2,13)
       bee(2,26)=deriv(2,13)
       bee(3,26)=deriv(1,13)
       bee(1,27)=deriv(1,14)
       bee(3,27)=deriv(2,14)
       bee(2,28)=deriv(2,14)
       bee(3,28)=deriv(1,14)
       bee(1,29)=deriv(1,15)
       bee(3,29)=deriv(2,15)
       bee(2,30)=deriv(2,15)
       bee(3,30)=deriv(1,15)
       bee(1,31)=deriv(1,16)
       bee(3,31)=deriv(2,16)
       bee(2,32)=deriv(2,16)
       bee(3,32)=deriv(1,16)
       !bee(4,1:values*2-1:2)=SFUNCT(:)/gc(1)
       bee(4,1:values*2-1:2)=0.0

      ELSE
       bee(1,1)=deriv(1,1)
       bee(3,1)=deriv(2,1)
       bee(2,2)=deriv(2,1)
       bee(3,2)=deriv(1,1)
       bee(1,3)=deriv(1,2)
       bee(3,3)=deriv(2,2)
       bee(2,4)=deriv(2,2)
       bee(3,4)=deriv(1,2)
       bee(1,5)=deriv(1,3)
       bee(3,5)=deriv(2,3)
       bee(2,6)=deriv(2,3)
       bee(3,6)=deriv(1,3)
       bee(1,7)=deriv(1,4)
       bee(3,7)=deriv(2,4)
       bee(2,8)=deriv(2,4)
       bee(3,8)=deriv(1,4)
       bee(1,9)=deriv(1,5)
       bee(3,9)=deriv(2,5)
       bee(2,10)=deriv(2,5)
       bee(3,10)=deriv(1,5)
       bee(1,11)=deriv(1,6)
       bee(3,11)=deriv(2,6)
       bee(2,12)=deriv(2,6)
       bee(3,12)=deriv(1,6)
       bee(1,13)=deriv(1,7)
       bee(3,13)=deriv(2,7)
       bee(2,14)=deriv(2,7)
       bee(3,14)=deriv(1,7)
       bee(1,15)=deriv(1,8)
       bee(3,15)=deriv(2,8)
       bee(2,16)=deriv(2,8)
       bee(3,16)=deriv(1,8)
       bee(1,17)=deriv(1,9)
       bee(3,17)=deriv(2,9)
       bee(2,18)=deriv(2,9)
       bee(3,18)=deriv(1,9)
       !bee(4,1:values*2-1:2)=SFUNCT(:)/gc(1)
       bee(4,1:values*2-1:2)=0.0
       !pause
      END IF

     RETURN
     END SUBROUTINE beemat3


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

   SUBROUTINE alphavalue2(cohesion,phi,stress_int,stress_ext,a_slo,alpha_1)
   USE MAIN
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   REAL(iwp),INTENT(IN)::cohesion,phi
   REAL(iwp),INTENT(OUT)::alpha_1
   REAL(iwp),INTENT(IN)::stress_int(:),stress_ext(:),a_slo
   REAL(iwp)::sigm,dsbar,lode_theta,fa,fb,d3=3.0,ps1,ps2,ps3,Ktheta,avalMC,one=1.0_iwp
 
   CALL invar_MC(stress_int,sigm,dsbar,lode_theta)
  ps1=sigm+(2.0/3.0)*2.0_iwp/SQRT(3.0_iwp)*dsbar*sin(lode_theta+(2.0*3.1415926/3.0))  !2Pi/3 = 120 degrees
  ps2=sigm+(2.0/3.0)*2.0_iwp/SQRT(3.0_iwp)*dsbar*sin(lode_theta)
  ps3=sigm+(2.0/3.0)*2.0_iwp/SQRT(3.0_iwp)*dsbar*sin(lode_theta-(2.0*3.1415926/3.0))
  Ktheta=cos(lode_theta)-one/SQRT(3.0_iwp)*sin(phi*3.1415926/180.0_iwp)*sin(lode_theta)
  avalMC=a_slo*cohesion/tan(phi*3.1415926/180.0_iwp)
  fa=sigm+SQRT(dsbar**2*Ktheta**2+avalMC**2*sin(phi*3.1415926/180.0_iwp))-cohesion*cos(phi*3.1415926/180.0_iwp)!Hyperbolic Mohr-Coulomb
  

   CALL invar_MC(stress_ext,sigm,dsbar,lode_theta)
  ps1=sigm+(2.0/3.0)*2.0_iwp/SQRT(3.0_iwp)*dsbar*sin(lode_theta+(2.0*3.1415926/3.0))  !2Pi/3 = 120 degrees
  ps2=sigm+(2.0/3.0)*2.0_iwp/SQRT(3.0_iwp)*dsbar*sin(lode_theta)
  ps3=sigm+(2.0/3.0)*2.0_iwp/SQRT(3.0_iwp)*dsbar*sin(lode_theta-(2.0*3.1415926/3.0))
  Ktheta=cos(lode_theta)-one/SQRT(3.0_iwp)*sin(phi*3.1415926/180.0_iwp)*sin(lode_theta)
  avalMC=a_slo*cohesion/tan(phi*3.1415926/180.0_iwp)
  fb=sigm+SQRT(dsbar**2*Ktheta**2+avalMC**2*sin(phi*3.1415926/180.0_iwp))-cohesion*cos(phi*3.1415926/180.0_iwp)
  
      !CALL invar(stress_int,sigm,dsbar,lode_theta)
      !CALL mocouf(phi,cohesion,sigm,dsbar,lode_theta,fa)
      !fa=dsbar-SQRT(d3)*cohesion
      !CALL invar(stress_ext,sigm,dsbar,lode_theta)
      !CALL mocouf(phi,cohesion,sigm,dsbar,lode_theta,fb)
      !fb=dsbar-SQRT(d3)*cohesion
      IF(fa>0.0)THEN
        !IF(fa>0.001)PAUSE
        alpha_1=0.0
      ELSE
        alpha_1=((-1.0)*fa)/(fb-fa)
      END IF
      !PRINT*,'alpha_1:',alpha_1

   END SUBROUTINE alphavalue2
   
      SUBROUTINE alphavalue_Von(cohesion,stress_int,stress_ext,alpha_1)
   USE MAIN
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   REAL(iwp),INTENT(IN)::cohesion
   REAL(iwp),INTENT(OUT)::alpha_1
   REAL(iwp),INTENT(IN)::stress_int(:),stress_ext(:)
   REAL(iwp)::sigm,dsbar,lode_theta,fa,fb,d3=3.0
 
   CALL invar(stress_int,sigm,dsbar,lode_theta)
     fa=dsbar-SQRT(d3)*cohesion
     CALL invar(stress_ext,sigm,dsbar,lode_theta)
     fb=dsbar-SQRT(d3)*cohesion
      IF(fa>0.0)THEN
        !IF(fa>0.001)PAUSE
        alpha_1=0.0
      ELSE
        alpha_1=((-1.0)*fa)/(fb-fa)
      END IF
      !PRINT*,'alpha_1:',alpha_1

 END SUBROUTINE alphavalue_Von

 SUBROUTINE alphavalue_Mohr(phif,psi,a_slo,Tload,cohesion,stress_int,stress_ext,alpha_1)
   !--Subroutine to compute alpha_1 value. Alpha_1 is the factor that make the stresses
   !--lies on the yield surface f=0, using a forward Euler approach (pag. 243 Smith book)
   USE MAIN
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   REAL(iwp),INTENT(IN)::cohesion,phif,psi,a_slo,Tload
   REAL(iwp),INTENT(OUT)::alpha_1
   REAL(iwp),INTENT(IN)::stress_int(:),stress_ext(:)
   REAL(iwp)::sigm,dsbar,lode_theta,fa,fb,d3=3.0,KT

      CALL invar(stress_int,sigm,dsbar,lode_theta)
      CALL mocouf_Sloan(phif,cohesion,sigm,dsbar,lode_theta,fa,a_slo,psi,Tload,KT)
      CALL invar(stress_ext,sigm,dsbar,lode_theta)
      CALL mocouf_Sloan(phif,cohesion,sigm,dsbar,lode_theta,fb,a_slo,psi,Tload,KT)
      IF(fa>0.0)THEN   !--stressses already on the yield surface and so the stress increase is plastic
        alpha_1=0.0
      ELSE
        alpha_1=((-1.0)*fa)/(fb-fa)
      END IF   
  
 END SUBROUTINE alphavalue_Mohr

    subroutine formdlambda(flowf,flowg,dee,deps,a,dlambda)
   USE MAIN
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   REAL(iwp),INTENT(IN)::dee(:,:),deps(:),a,flowf(:),flowg(:)
   REAL(iwp),INTENT(OUT)::dlambda
   REAL(iwp)::top,bot
   INTEGER::nst
   REAL(iwp),ALLOCATABLE::caflow(:),dsigmae(:)
    nst=UBOUND(deps,1)
   ALLOCATE(caflow(nst),dsigmae(nst))
   caflow=MATMUL(dee,flowg)
   bot=dot_product(flowf,caflow)
   bot=bot+a
   dsigmae=matmul(dee,deps)
   top=dot_product(flowf,dsigmae)
   dlambda=top/bot
   deallocate(caflow,dsigmae)
   RETURN
    end subroutine formdlambda
    
        subroutine formdlambda2(flowf,flowg,dee,deps,a,dlambda)
   USE MAIN
   IMPLICIT NONE
   INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
   REAL(iwp),INTENT(IN)::dee(:,:),deps(:),a,flowf(:),flowg(:)
   REAL(iwp),INTENT(OUT)::dlambda
   REAL(iwp)::top,bot
   INTEGER::nst
   REAL(iwp),ALLOCATABLE::caflow(:),dsigmae(:)
    nst=UBOUND(deps,1)
   ALLOCATE(caflow(nst),dsigmae(nst))
   caflow=MATMUL(dee,flowg)
   bot=dot_product(flowf,caflow)
   bot=bot*a
   dsigmae=matmul(dee,deps)
   top=dot_product(flowf,dsigmae)
   dlambda=top/bot
   deallocate(caflow,dsigmae)
   RETURN
  end subroutine formdlambda2

!    SUBROUTINE ycorrec(softval,sigm,dsbar,lode_theta,stress,epsinv,cp,gamap2,cr,sp,dee,phif, &
!                     psi,psip,psires,a_slo,TLoad,coh)
!    USE MAIN
!     IMPLICIT NONE
!     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
!     REAL(iwp),INTENT(IN)::sp,cp,gamap2,cr,dee(:,:),a_slo,phif,TLoad, &
!     						psires,psip,coh,psi
!     INTEGER,INTENT(IN)::softval
!     REAL(iwp),INTENT(IN OUT)::stress(:),epsinv,sigm,dsbar,lode_theta
!     REAL(iwp)::bot,dnamda,dq1,dq2,dq3,a,f_trial,f0,ftol=1.0d-06,     &
!     			zero=0.0d0,d3=3.0d0,two=2.0d0,epsinv_0,epsinv_trial,h,cf,theta,pt5=0.5_iwp, &
!                fyield,dFdc,KT,dQ,fnew,psiaux
!     INTEGER::nst,iters,maxits=10
!     REAL(iwp),ALLOCATABLE::caflow(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),flowf(:),flowg(:), &
!     						stress0(:),stress_trial(:),vmfl(:)
!
!     
!     nst=UBOUND(stress,1)
!     ALLOCATE(caflow(nst),m1(nst,nst),m2(nst,nst),m3(nst,nst),flowf(nst),flowg(nst),flow(nst,nst), &
!     stress0(nst),stress_trial(nst),vmfl(nst))
!     stress0=stress
!     epsinv_0=epsinv
! 
!   psiaux=psi
!   cf=coh
!   iters=0
!   iters_yield: DO
!    iters=iters+1
!
!     CALL invar(stress0,sigm,dsbar,lode_theta)
!     CALL mocouf_Sloan(phif,cf,sigm,dsbar,lode_theta,fyield,a_slo,psiaux,Tload,KT)
!     CALL plasticmod(epsinv_0,gamap2,sp,h)
!     CALL mocouq_Sloan(cf,psiaux,dsbar,lode_theta,dq1,dq2,dq3,TLoad,a_slo)
!     CALL formm(stress0,m1,m2,m3)
!     flowf=MATMUL((m1*dq1+m2*dq2+m3*dq3),stress0)
!     flowg= flowf
!     dQ = pt5*((flowf(1) - flowf(2))**2 + (flowf(2)-flowf(4))**2 + (flowf(4)-flowf(1))**2) + 0.75_iwp*flowf(3)**2
!     CALL form_Ai_Mohr(softval,dQ,h,a_slo,dsbar,sigm,cf,phif,psiaux,KT,lode_theta,TLoad,a,dFdc)
!     caflow=MATMUL(dee,flowg)
!     bot=dot_product(flowf,caflow)
!     dnamda=fyield/(a+bot)
!
!     stress_trial=stress0-dnamda*caflow
!     epsinv_trial=epsinv_0+dnamda*(2.0/3.0)*sqrt(dQ)
!     CALL invar(stress_trial,sigm,dsbar,lode_theta)
!     CALL mocouf_Sloan(phif,cf,sigm,dsbar,lode_theta,f_trial,a_slo,psiaux,Tload,KT)
!     
!     
!!$$$$$$      IF(softval==1)cf=epsinv_0*h+cp
!!$$$$$$      IF(softval==2)psiaux=epsinv_0*h+psip
!!$$$$$$      IF(cf<cr)cf=cr;IF(psiaux<psires)psiaux=psires
!!$$$$$$      IF(cf>cp)cf=cp;IF(psiaux>psip)psiaux=psip
!     
! 
!     
!    
!!$$$$$$      CALL plasticmod(epsinv_trial,gamap2,sp,h)
!!$$$$$$      IF(softval==1)cf=epsinv_trial*h+cp
!!$$$$$$      IF(softval==2)psiaux=epsinv_trial*h+psip
!!$$$$$$      IF(cf<cr)cf=cr;IF(psiaux<psires)psiaux=psires
!!$$$$$$      IF(cf>cp)cf=cp;IF(psiaux>psip)psiaux=psip
!     
!     
!      IF(abs(f_trial)>abs(fyield))THEN
!       bot=dot_product(flowf,flowf)
!       dnamda=fyield/bot
!       stress_trial=stress0-dnamda*flowf
!       epsinv_trial=epsinv_0
!!$$$$$$        CALL invar(stress_trial,sigm,dsbar,lode_theta)
!!$$$$$$        CALL plasticmod(epsinv_trial,gamap2,sp,h)
!!$$$$$$         IF(softval==1)cf=epsinv_trial*h+cp
!!$$$$$$      IF(softval==2)psiaux=epsinv_trial*h+psip
!!$$$$$$      IF(cf<cr)cf=cr;IF(psiaux<psires)psiaux=psires
!!$$$$$$      IF(cf>cp)cf=cp;IF(psiaux>psip)psiaux=psip
!!$$$$$$        CALL mocouf_Sloan(phif,cf,sigm,dsbar,lode_theta,f_trial,a_slo,psiaux,Tload,KT)  
!      END IF
!      stress0=stress_trial
!      epsinv_0=epsinv_trial
!      
!      IF(abs(f_trial)<=ftol)EXIT
!     
!      IF(iters==maxits)THEN 
!        write(*,*) "correction stress not convergent"
!      END IF  
!      IF(iters==maxits)EXIT
!     END DO iters_yield
!
!      stress=stress_trial;epsinv=epsinv_trial
!      ! k0
!!$$$$$$      CALL invar(stress,sigm,dsbar,lode_theta)
!!$$$$$$      CALL plasticmod(epsinv_trial,gamap2,sp,h)
!!$$$$$$      IF(softval==1)coh=epsinv_trial*h+cp
!!$$$$$$      IF(softval==2)psi=epsinv_trial*h+psip
!!$$$$$$      IF(coh<cr)coh=cr;IF(psi<psires)psi=psires
!!$$$$$$      IF(coh>cp)coh=cp;IF(psi>psip)psi=psip
! 
!     DEALLOCATE(caflow,m1,m2,m3,flowf,flowg,flow,stress0,stress_trial)
!     RETURN
!   END SUBROUTINE ycorrec

   SUBROUTINE yieldcorrection(stress,k,ci,gamap1,cp,gamap2,cr,hp,sp,dee,ft)
    USE MAIN
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::hp,sp,ci,gamap1,cp,gamap2,cr,dee(:,:)
     REAL(iwp),INTENT(IN OUT)::stress(:),k,ft
     REAL(iwp)::bot,dnamda,sigm,dsbar,lode_theta,dq1,dq2,dq3,a,f_trial,f0,ftol=1.0d-06, &
     zero=0.0d0,d3=3.0d0,two=2.0d0,k0,k_trial,h,cf,theta
     INTEGER::nst,iters,maxits=10
     REAL(iwp),ALLOCATABLE::caflow(:),m1(:,:),m2(:,:),m3(:,:),flow(:,:),flowf(:),flowg(:), &
     stress0(:),stress_trial(:),vmfl(:)
     nst=UBOUND(stress,1)
     ALLOCATE(caflow(nst),m1(nst,nst),m2(nst,nst),m3(nst,nst),flowf(nst),flowg(nst),flow(nst,nst), &
     stress0(nst),stress_trial(nst),vmfl(nst))
     stress0=stress
     k0=k
     iters=0
     !CALL formh(k0,gamap1,gamap2,hp,sp,h)
     CALL plasticmod(k0,gamap2,sp,h)
     a=sqrt(d3)*h
     iters_2: DO
      iters=iters+1
      !CALL deemat(dee,e,v)
      CALL invar(stress0,sigm,dsbar,lode_theta)
      !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k0)
      cf=k0*sp+cp
      IF(cf<cr)cf=cr
      f0=dsbar-sqrt(d3)*cf
      CALL vmflow(stress0,dsbar,vmfl)
      flowf=vmfl
      flowg= flowf
      !CALL formh(k0,gamap1,gamap2,hp,sp,h)
      CALL plasticmod(k0,gamap2,sp,h)
      a=sqrt(d3)*h
      caflow=MATMUL(dee,flowg)
      bot=dot_product(flowf,caflow)
      dnamda=f0/(a+bot)
      stress_trial=stress0-dnamda*caflow
      k_trial=k0+dnamda
      CALL invar(stress_trial,sigm,dsbar,lode_theta)
      !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k_trial)
      cf=k_trial*sp+cp
      IF(cf<cr)cf=cr
      f_trial=dsbar-sqrt(d3)*cf   ! cf=k_trial ;;;;;;  my assumption not sure...
      IF(abs(f_trial)>abs(f0))THEN
       bot=dot_product(flowf,flowf)
       dnamda=f0/bot
       stress_trial=stress0-dnamda*flowf
       k_trial=k0
       CALL invar(stress_trial,sigm,dsbar,lode_theta)
       !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k_trial)
       cf=k_trial*sp+cp
       IF(cf<cr)cf=cr
       f_trial=dsbar-sqrt(d3)*cf
      END IF
      IF(abs(f_trial)<=ftol) EXIT
      stress0=stress_trial
      k0=k_trial
      IF(iters==maxits) write(*,*) "correction stress not convergent"
     END DO iters_2

      stress=stress_trial;k=k_trial ! k0
      CALL invar(stress,sigm,dsbar,lode_theta)
      !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k)
      cf=k*sp+cp
      IF(cf<cr)cf=cr
      ft=dsbar-sqrt(d3)*cf
      DEALLOCATE(caflow,m1,m2,m3,flowf,flowg,flow,stress0,stress_trial)
      RETURN
   END SUBROUTINE yieldcorrection
   
   SUBROUTINE yieldcorrection_Sloan(stress,c,phi,dee,a_slo,Theta_t,k,ft)
    USE MAIN
     IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
     REAL(iwp),INTENT(IN)::c,dee(:,:),phi,a_slo,Theta_t
     REAL(iwp),INTENT(IN OUT)::stress(:),ft,k
     REAL(iwp)::bot,dnamda,sigm,dsbar,lode_theta,dq1,dq2,dq3,a,f_trial,f0,ftol=1.0d-06, &
     zero=0.0_iwp,d3=3.0_iwp,two=2.0_iwp,k0,k_trial,h,cf,theta,Ktheta,DTDT,avalMC
     INTEGER::nst,iters,maxits=10
     REAL(iwp),ALLOCATABLE::caflow(:),m1_s(:),m2_s(:),m3_s(:),flow(:,:),flowf(:),flowg(:), &
     stress0(:),stress_trial(:),vmfl(:)
     nst=UBOUND(stress,1)
     ALLOCATE(caflow(nst),m1_s(nst),m2_s(nst),m3_s(nst),flowf(nst),flowg(nst),flow(nst,nst), &
     stress0(nst),stress_trial(nst),vmfl(nst))
     stress0=stress
     k0=k
     iters=0
     !CALL formh(k0,gamap1,gamap2,hp,sp,h)
     !CALL plasticmod(k0,gamap2,sp,h)
     !a=sqrt(d3)*h
     a=zero
     iters_2: DO
      iters=iters+1
      CALL invar_MC(stress0,sigm,dsbar,lode_theta)
      CALL Sloan_KT(lode_theta,Theta_t,phi,Ktheta,DTDT)!Ktheta=cos(lode_theta)-one/SQRT(3.0_iwp)*sin(mbod(2)%phi_mp(i)*3.1415926/180.0_iwp)*sin(lode_theta)
      avalMC=a_slo*c/tan(phi*3.1415926_iwp/180.0_iwp)
      f0=sigm+SQRT(dsbar**2*Ktheta**2+avalMC**2*sin(phi*3.1415926_iwp/180.0_iwp))-c*cos(phi*3.1415926_iwp/180.0_iwp)!Hyperbolic Mohr-Coulomb
      !CALL deemat(dee,e,v)
      !CALL invar(stress0,sigm,dsbar,lode_theta)
      !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k0)
      !cf=k0*sp+cp
      !IF(cf<cr)cf=cr
      !f0=dsbar-sqrt(d3)*cf
      !CALL vmflow(stress0,dsbar,vmfl)
      !flowf=vmfl
      CALL formm_Sloan(stress0,dsbar,m1_s,m2_s,m3_s)
      CALL mocouq_Sloan(c,phi,dsbar,lode_theta,dq1,dq2,dq3,Ktheta,a_slo,DTDT)
      flowf=(m1_s*dq1+m2_s*dq2+m3_s*dq3)
      flowg= flowf
      !CALL formh(k0,gamap1,gamap2,hp,sp,h)
      !CALL plasticmod(k0,gamap2,sp,h)
      !a=sqrt(d3)*h
      a=zero
      caflow=MATMUL(dee,flowg)
      bot=dot_product(flowf,caflow)
      dnamda=f0/(a+bot)
      stress_trial=stress0-dnamda*caflow
      k_trial=k0+dnamda
      CALL invar_MC(stress_trial,sigm,dsbar,lode_theta)
      CALL Sloan_KT(lode_theta,Theta_t,phi,Ktheta,DTDT)!Ktheta=cos(lode_theta)-one/SQRT(3.0_iwp)*sin(mbod(2)%phi_mp(i)*3.1415926/180.0_iwp)*sin(lode_theta)
      avalMC=a_slo*c/tan(phi*3.1415926_iwp/180.0_iwp)
      f_trial=sigm+SQRT(dsbar**2*Ktheta**2+avalMC**2*sin(phi*3.1415926_iwp/180.0_iwp))-c*cos(phi*3.1415926_iwp/180.0_iwp)
      
      !CALL invar(stress_trial,sigm,dsbar,lode_theta)
      !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k_trial)
      !cf=k_trial*sp+cp
      !IF(cf<cr)cf=cr
      !f_trial=dsbar-sqrt(d3)*cf   ! cf=k_trial ;;;;;;  my assumption not sure...
      IF(abs(f_trial)>abs(f0))THEN
       bot=dot_product(flowf,flowf)
       dnamda=f0/bot
       stress_trial=stress0-dnamda*flowf
       k_trial=k0
       CALL invar_MC(stress_trial,sigm,dsbar,lode_theta)
       CALL Sloan_KT(lode_theta,Theta_t,phi,Ktheta,DTDT)!Ktheta=cos(lode_theta)-one/SQRT(3.0_iwp)*sin(mbod(2)%phi_mp(i)*3.1415926/180.0_iwp)*sin(lode_theta)
      avalMC=a_slo*c/tan(phi*3.1415926_iwp/180.0_iwp)
      f_trial=sigm+SQRT(dsbar**2*Ktheta**2+avalMC**2*sin(phi*3.1415926_iwp/180.0_iwp))-c*cos(phi*3.1415926_iwp/180.0_iwp)
       !CALL invar(stress_trial,sigm,dsbar,lode_theta)
       !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k_trial)
       !cf=k_trial*sp+cp
       !IF(cf<cr)cf=cr
       !f_trial=dsbar-sqrt(d3)*cf
      END IF
      IF(abs(f_trial)<=ftol) EXIT
      stress0=stress_trial
      k0=k_trial
      IF(iters==maxits) write(*,*) "correction stress not convergent"
     END DO iters_2

      stress=stress_trial;k=k_trial ! k0
      CALL invar_MC(stress_trial,sigm,dsbar,lode_theta)
       CALL Sloan_KT(lode_theta,Theta_t,phi,Ktheta,DTDT)!Ktheta=cos(lode_theta)-one/SQRT(3.0_iwp)*sin(mbod(2)%phi_mp(i)*3.1415926/180.0_iwp)*sin(lode_theta)
      avalMC=a_slo*c/tan(phi*3.1415926_iwp/180.0_iwp)
      ft=sigm+SQRT(dsbar**2*Ktheta**2+avalMC**2*sin(phi*3.1415926_iwp/180.0_iwp))-c*cos(phi*3.1415926_iwp/180.0_iwp)
      !CALL invar(stress,sigm,dsbar,lode_theta)
      !CALL vm_hsr(ci,gamap1,cp,gamap2,cr,cf,sigm,dsbar,theta,k)
      !cf=k*sp+cp
      !IF(cf<cr)cf=cr
      !ft=dsbar-sqrt(d3)*cf
      DEALLOCATE(caflow,m1_s,m2_s,m3_s,flowf,flowg,flow,stress0,stress_trial)
      RETURN
    END SUBROUTINE yieldcorrection_Sloan

    FUNCTION euclid_distance(A,B) RESULT(C)
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::A(:),B(:)
        REAL(iwp) ::C
        C = sqrt(sum((A-B)**2))
    END FUNCTION

!** FUNCTION double_dot takes the double dot product of two tensors
    FUNCTION double_dot(A,B) RESULT (C)
        IMPLICIT none
        integer,parameter:: iwp=selected_real_kind(15)
        real(iwp),intent(in):: A(:,:),B(:,:)
        integer::i,j
        real(iwp):: C
        C=0.0_iwp
        DO i=1,size(A,1)
            DO j=1,size(A,2)
                C=C+A(i,j)*B(i,j)
            END DO
        END DO
    RETURN
    END FUNCTION

!** FUNCTION bilinreg to perform a bilinear regression and fit a plane to data
    FUNCTION bilinreg(x,y,z) RESULT(res)
        ! The function performs an ordinary least-squares linear regression in 2 dimensions for input data x,y,z
        ! It is used to calculate the gradient of density in the support of a contact node in the resolution of contact
        ! The normal direction (normal to the body surface) is taken as the opposite direction to the gradient of density
        ! How it works:
        ! X = [x y 1]; beta = [a b c]'
        ! X.beta=z => X'.X.beta=X'.z => beta = inv(X'.X).X'.z
        USE main
        IMPLICIT none
        INTEGER,PARAMETER:: iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::x(:),y(:),z(:)
        INTEGER::i
        REAL(iwp)::avg_x2,avg_xy,avg_x,avg_y2,avg_y,avg_xz,avg_yz,avg_z
        REAL(iwp)::mat(3,3),vec(3),res(3)
        !Calculate terms for X'X matrix and X'z vector
        avg_x2=sum(x**2)/size(x)
        avg_xy=sum(x*y)/size(x)
        avg_x=sum(x)/size(x)
        avg_y2=sum(y**2)/size(y)
        avg_y=sum(y)/size(y)
        avg_xz=sum(x*z)/size(x)
        avg_yz=sum(y*z)/size(y)
        avg_z=sum(z)/size(z)
        !Assemble X'X matrix
        mat=transpose(reshape((/avg_x2, avg_xy, avg_x, avg_xy, avg_y2, avg_y, avg_x, avg_y, 1.0_iwp/),(/size(mat,1),size(mat,2)/)))
        !Assemble X'z vector
        vec=(/avg_xz, avg_yz, avg_z/)
        CALL invert(mat)
        res=matmul(mat,vec)
    RETURN
    END FUNCTION

!** FUNCTION b_spline used to calculate the quadratic b-spline for interpolating densities to support cell centres of a boundary node
    FUNCTION b_spline(point,center,cellsize) RESULT(S)
        !Funtion to interpolate the mass to the centres of cells in the support of boundary nodes-----------
        !Inputs:
        ! point - coordinate of the material point
        ! center - coordinate of the cell center
        ! cellsize - size of the cell
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::point,center,cellsize
        REAL(iwp)::distance
        REAL(iwp)::S
        distance=center-point
            IF ((-3/(2*cellsize).le.distance) .and. (distance .le.-1/(2*cellsize))) THEN
                S = 1*(distance**2)/(2*(cellsize**2)) + 3*distance/(2*cellsize) + 9/8
            ELSEIF ((-1/(2*cellsize).lt.distance) .and. (distance .lt.1/(2*cellsize))) THEN
                S = -1*(distance**2)/(cellsize**2) + 3/4
            ELSEIF ((1/(2*cellsize).le.distance) .and. (distance .le.3/(2*cellsize))) THEN
                S = 1*(distance**2)/(2*(cellsize**2)) - 3*distance/(2*cellsize) + 9/8
            ELSE
                S = 0
            END IF
    RETURN
    END FUNCTION b_spline

!** SUBROUTINE ball_geometry used to define a ball with local subdivision
    SUBROUTINE ball_geometry(el_in_r)
        IMPLICIT NONE
        INTEGER:: el_in_r

    END SUBROUTINE ball_geometry

!** SUBROUTINE normals used to calculate the body normals at the contact nodes
    SUBROUTINE contact_node_normals(normal,contact_node_index,node_support,gm_coord,g_coord,g_num, &
        m_mass,a_ele,cellsize,theta)
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::gm_coord(:,:),g_coord(:,:),m_mass(:),cellsize,theta
        INTEGER,INTENT(IN)::contact_node_index,node_support(:,:),g_num(:,:),a_ele(:)
        INTEGER::i,j,iel
        INTEGER::support(size(node_support,2)),active_support(size(node_support,2)),num(size(g_num,1)),nmps
        REAL(iwp)::centers(size(node_support,2),2),densities(size(node_support,2)),coord(size(g_coord,1),2), &
        masses(size(node_support,2)),volume(size(node_support,2)),normal(:),coeffs(3)
        centers=0.0_iwp; volume=0.0_iwp; normal=0.0_iwp
        active_support=0
        densities=0.0_iwp
        nmps=size(m_mass)
        centers=0.0_iwp
        support = node_support(contact_node_index,:) !Find the elements in the nodal support of the current node
        !Interpolate density onto support cell centres
        DO i = 1,size(support,1)
            iel = support(i) ! Get element number
            IF (iel/=0) THEN
                centers(i,1)=sum(g_coord(1,g_num(:,iel)))/size(support,1) !Get x-coord of current element center
                centers(i,2)=sum(g_coord(2,g_num(:,iel)))/size(support,1) !Get y-coord of current element center
                DO j = 1,nmps
                    !IF (abs(centers(i,1)-gm_coord(1,j)).le.1.5*cellsize*theta .AND. &
                    !    abs(centers(i,2)-gm_coord(2,j)).le.1.5*cellsize) then
                    IF (a_ele(j)==iel) THEN
                        densities(i)=densities(i)+(b_spline(gm_coord(1,j),centers(i,1),cellsize*theta)* &
                        b_spline(gm_coord(2,j),centers(i,2),cellsize)*m_mass(j))/(theta*cellsize**2)
                        !masses(i)=masses(i)+(b_spline(gm_coord(1,j),centers(i,1),cellsize*theta)* &
                        !b_spline(gm_coord(2,j),centers(i,2),cellsize)*m_mass(j))
                    END IF
                END DO
            END IF
        END DO
        !densities=masses/(theta*(cellsize**2))
        !Perform bilinear regression to fit plane to density data
        coeffs=bilinreg(centers(:,1),centers(:,2),densities)
        !Normalise normal vector
        normal(1)=-coeffs(1)/sqrt(sum(coeffs(1:2)**2))
        normal(2)=-coeffs(2)/sqrt(sum(coeffs(1:2)**2))
    END SUBROUTINE contact_node_normals


   SUBROUTINE tang_dir(m,node_support,mpoints,m_velocity,tangent)

    USE main
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::mpoints(:,:),m_velocity(:,:)
    INTEGER,INTENT(IN)::node_support(:,:),m
    REAL(iwp),INTENT(IN OUT)::tangent(:,:)



    
   END SUBROUTINE tang_dir

!** SUBROUTINE mass_gradient_normal used to find the body normal at contact nodes from mass gradient
    SUBROUTINE mass_gradient_normal(normal,contact_node,node_support,g_num,mpoints1,mpoints2, &
    								a_ele,elemcont,d_ele,c_ele,nf)

    USE main
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::mpoints1(:,:),mpoints2(:,:)
    INTEGER,INTENT(IN)::g_num(:,:),contact_node,node_support(:,:),a_ele(:),elemcont(:),d_ele(:),c_ele(:),nf(:,:)
    REAL(iwp),INTENT(IN OUT)::normal(:,:)
    REAL(iwp)::magnorm,zero=0.0_iwp
    REAL(iwp)::der(SIZE(normal,1),SIZE(g_num,1)),fun(4),shapeacum(4) !Holds the shape function derivatives
    REAL(iwp),ALLOCATABLE::loccoord(:,:)
    INTEGER::i,ind,num(SIZE(g_num,1)),iel,auxel,contel,m,ielint,nmps,k,mp,ndim,direction,elcount,minpoints
    INTEGER::saveel(4)
    LOGICAL::freebound
    ALLOCATE(loccoord(1,2))
    loccoord(1,:)=0.0_iwp
    der=0.0_iwp
    m=contact_node
    auxel=0
    nmps=SIZE(mpoints2,1)
    IF(SIZE(mpoints1,1)>SIZE(mpoints1,1))nmps=SIZE(mpoints1,1)
    shapeacum=0.0_iwp
    ndim=SIZE(mpoints1,2)
    saveel=0
 

    !-Detect the contact condition(center contact, boundary, inside the body, outside, etc)
    freebound=.false.
    elcount=0
    k=1
    DO i=1,2**ndim
     iel=node_support(m,i) !Element of material point
      IF(iel>0)THEN
     IF((nf(1,g_num(1,iel))>=1.or.nf(2,g_num(1,iel))>=1).and.(nf(1,g_num(2,iel))>=1.or.nf(2,g_num(2,iel))>=1).and. &
       (nf(1,g_num(3,iel))>=1.or.nf(2,g_num(3,iel))>=1).and.(nf(1,g_num(4,iel))>=1.or.nf(2,g_num(4,iel))>=1))THEN !-Find if the element is active
    
          IF(d_ele(iel)==1)THEN
            elcount=elcount+1
            saveel(k)=iel
            k=k+1
          END IF 
           
       IF(d_ele(iel)==0)freebound=.true.
         
      END IF 
     END IF  
    END DO
    !--If freebound=false then the node is inside the body
    
  IF(freebound)THEN
    DO i=1,2**ndim
     iel=node_support(m,i) !Element of material point 
     IF(iel>=1.and.c_ele(iel)>=1)THEN
       num=g_num(:,iel)
       auxel=iel
       DO ind=1,SIZE(num,1)
         IF(num(ind)==m) THEN
            EXIT
         END IF
       END DO
         CALL shape_der(der,loccoord,1)
         IF(nf(1,m)==0)der(1,:)=0.0 !-This option erase the normal vector if the node is fix in this direction
         IF(nf(2,m)==0)der(2,:)=0.0
         normal(:,m)=normal(:,m)+der(:,ind)!*m_mass(i) !Add contribution of current mp to the node mass
      END IF
    END DO
     
  !-If magnorm>0.0 then the bodies are in different elemenmts, otherwise the bodies are in the same element
    magnorm=SQRT(SUM(normal(:,m)**2))
    IF(magnorm>0.0)THEN
      normal(:,m)=normal(:,m)/SQRT(SUM(normal(:,m)**2)) !Normalise normal vector
    END IF
     
  ELSE 
         
   !-Bodies in the same element and the node is surounded by one bodie
   IF(elcount==4)THEN !if elcount==4 the node is surounded by material points in the four elements
     IF(c_ele(node_support(m,1))+c_ele(node_support(m,2))>                           &
                        c_ele(node_support(m,3))+c_ele(node_support(m,4)))THEN

             normal(2,m)=-1.0!There are more material points in the two upper elements
     ELSE IF(c_ele(node_support(m,1))+c_ele(node_support(m,2))<                           &
                        c_ele(node_support(m,3))+c_ele(node_support(m,4)))THEN
               
  			normal(2,m)=1.0!There are more material points in the two lower elements
                        
     END IF

     IF(c_ele(node_support(m,1))+c_ele(node_support(m,3))>                           &
                        c_ele(node_support(m,2))+c_ele(node_support(m,4)))THEN

             normal(1,m)=1.0!There are more material points in the two left elements
     ELSE IF(c_ele(node_support(m,1))+c_ele(node_support(m,3))<                           &
                        c_ele(node_support(m,2))+c_ele(node_support(m,4)))THEN
               
  			normal(1,m)=-1.0!There are more material points in the two right elements
                        
     END IF

   ELSE  
    direction=2
    IF(saveel(1)+1==saveel(2))direction=1
     IF(direction==1)THEN
        IF(c_ele(saveel(1))>c_ele(saveel(2)))THEN
         normal(1,m)=1.0!There are more material points in the left elements
        ELSE! IF(c_ele(node_support(m,2))>c_ele(node_support(m,1)))THEN 
         normal(1,m)=-1.0!There are more material points in the right elements
        END IF  
        IF(nf(1,m)==0)normal(1,m)=0.0
      ELSE 
 		IF(c_ele(saveel(1))>c_ele(saveel(2)))THEN
         normal(2,m)=-1.0!There are more material points in the left elements
        ELSE! IF(c_ele(node_support(m,2))>c_ele(node_support(m,1)))THEN 
         normal(2,m)=1.0!There are more material points in the right elements
        END IF 
        IF(nf(2,m)==0)normal(2,m)=0.0
     END IF  
   END IF  
 END IF
     
  END SUBROUTINE mass_gradient_normal

  SUBROUTINE normal_dir(normal,contact_node,node_support,g_num,mpoints1,mpoints2, &
    								a_ele1,a_ele2,nf,c_ele1,c_ele2)

    USE main
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),INTENT(IN)::mpoints1(:,:),mpoints2(:,:)
    INTEGER,INTENT(IN)::g_num(:,:),contact_node,node_support(:,:),nf(:,:)
    INTEGER,INTENT(IN)::a_ele1(:),a_ele2(:),c_ele1(:),c_ele2(:)
    REAL(iwp),INTENT(IN OUT)::normal(:,:)
    REAL(iwp)::magnorm,zero=0.0_iwp,div,normal_acum1(size(nf,1)),normal_acum2(size(nf,1)),normal_acum(size(nf,1))
    REAL(iwp)::derb1(SIZE(normal,1),SIZE(g_num,1)),fun(4),shapeacum(4) !Holds the shape function derivatives
    REAL(iwp),ALLOCATABLE::loccoord(:,:)
    INTEGER::i,ind,num(SIZE(g_num,1)),auxel,contel,m,ielint,nmps,k,mp,ndim,direction,minpoints
    INTEGER::iel1,iel2,countp,nmp1,nmp2
    INTEGER::saveel(4)
    LOGICAL::freebound

    derb1=0.0_iwp
    m=contact_node
    normal_acum=zero
    normal_acum1=zero
    normal_acum2=zero
    countp=0
    nmp1=0
    nmp2=0

    DO i=1,4
    nmp1=nmp1+c_ele1(node_support(m,i)) !-Total amount of material points of body 1 in the suport domain
    nmp2=nmp2+c_ele2(node_support(m,i)) !-Total amount of material points of body 2 in the suport domain
    END DO
    
    Loop_1:DO i=1,SIZE(a_ele1,1)
     iel1=a_ele1(i)  !Element of material point 
     IF(ANY(node_support(m,:)==iel1)) THEN
        
        num=g_num(:,iel1)
        DO ind=1,SIZE(num,1)
          IF(num(ind)==m) THEN
             EXIT
          END IF
        END DO
          CALL shape_der2(derb1,mpoints1,i)    !Calculate shape function derivatives for each material point; der stores the shape function derivative values
          normal_acum1=normal_acum1+derb1(:,ind) !Add contribution of current mp to the node mass
          IF(nf(1,m)==0)normal_acum(1)=0.0 !-This option erase the normal vector if the node is fix in this direction
          IF(nf(2,m)==0)normal_acum(2)=0.0
          countp=countp+1
          IF(countp==nmp1)EXIT
      END IF
    END DO Loop_1

    magnorm=SQRT(SUM(normal_acum1(:)**2))
    IF(magnorm>0.0)THEN
    normal_acum1=normal_acum1/magnorm
    ELSE
      normal_acum1=zero
    END IF  
    countp=0
   Loop_2:DO i=1,SIZE(a_ele2,1)
     iel2=a_ele2(i)  !Element of material point 
     IF(ANY(node_support(m,:)==iel2)) THEN
      
        num=g_num(:,iel2)
        DO ind=1,SIZE(num,1)
          IF(num(ind)==m) THEN
             EXIT
          END IF
        END DO
          CALL shape_der2(derb1,mpoints2,i)    !Calculate shape function derivatives for each material point; der stores the shape function derivative values
          normal_acum2=normal_acum2-derb1(:,ind) !Add contribution of current mp to the node mass
          IF(nf(1,m)==0)normal_acum(1)=0.0 !-This option erase the normal vector if the node is fix in this direction
          IF(nf(2,m)==0)normal_acum(2)=0.0
          countp=countp+1
          IF(countp==nmp2)EXIT
      END IF
    END DO Loop_2

   magnorm=SQRT(SUM(normal_acum2(:)**2))
   IF(magnorm>0.0)THEN
   normal_acum2=normal_acum2/magnorm
   ELSE
     normal_acum2=zero
   END IF  

    normal_acum=normal_acum1+normal_acum2
     
    !-If magnorm>0.0 then the bodies are in different elemenmts, otherwise the bodies are in the same element
    magnorm=SQRT(SUM(normal_acum(:)**2))
    IF(magnorm>0.0)THEN
      normal(:,m)=normal_acum(:)/magnorm !Normalise normal vector
    END IF
     
  
     
  END SUBROUTINE normal_dir

  

!** SUBROUTINE ini_stress used to initialise stress in the bodies (hydrostatic/K0-procedure)
    SUBROUTINE ini_stress(gm_coord,m_stress,rho,nu,g_matrix,datum,K0)
        !This subroutine initialises vertical and horizontal geostatic stresses based on an optional input of K0, for a body
        ! w.r.t. a datum (top of the body).
        !ASSUMES UNIFORM MATERIAL PROPERTIES FOR THE BODY
        !GRAVITY ACTING IN VERTICAL DIRECTION ONLY
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        INTEGER::i
        REAL(iwp),INTENT(IN)::datum,g_matrix(:),gm_coord(:,:),nu(:),rho(:)
        REAL(iwp)::m_stress(:,:)
        REAL(iwp),OPTIONAL::K0
        m_stress=0.0_iwp
        IF (PRESENT(K0)) THEN
            DO i=1,size(m_stress,2)
                m_stress(2,i)=abs(gm_coord(2,i)-datum)*rho(i)*g_matrix(2)
                m_stress(1,i)=K0*m_stress(2,i)
            END DO
        ELSE
            DO i=1,size(m_stress,2)
                m_stress(2,i)=abs(gm_coord(2,i)-datum)*rho(i)*g_matrix(2)
                m_stress(1,i)=((1-nu(i))/nu(i))*m_stress(2,i)
            END DO
        END IF
    RETURN
    END SUBROUTINE ini_stress

!** SUBORUTINE deemat_2 used to assemble the plane-stress stress-strain matrix
    SUBROUTINE deemat_2(dee,e,v)
        ! This subroutine returns the elastic dee matrix for ih=3 (plane strain),
        ! ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6
        ! (three dimensions).
        IMPLICIT NONE
        INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        REAL(iwp),INTENT(IN)::e,v
        REAL(iwp),INTENT(OUT)::dee(:,:)
        REAL(iwp)::v1,v2,c,vv,zero=0.0_iwp,pt5=0.5_iwp,one=1.0_iwp,two=2.0_iwp
        INTEGER::i,ih
        dee=zero
        ih=UBOUND(dee,1)
        v1=one-v
        c=e/(1-v**2)
        SELECT CASE(ih)
            CASE(3)
                dee(1,1)=1*c
                dee(2,2)=1*c
                dee(1,2)=v*c
                dee(2,1)=v*c
                dee(3,3)=v1
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
    END SUBROUTINE deemat_2

!** SUBROUTINE read_mps used to read mp locations from a file
    SUBROUTINE read_mps(filename,mp_coord,mp_num)
        IMPLICIT NONE
        INTEGER, PARAMETER::iwp=selected_real_kind(15)
        CHARACTER(*),intent(in)::filename
        INTEGER,INTENT(IN) :: mp_num
        INTEGER::i=1
        REAL(iwp),INTENT(OUT)::mp_coord(:,:)
        OPEN(850,FILE=filename//'.dat')
        DO i=1,mp_num
            READ(850,*) mp_coord(:,i)
        END DO
        CLOSE(850)
        RETURN
    END SUBROUTINE


end module mpm
