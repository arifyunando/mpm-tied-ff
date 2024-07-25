module utilities
  implicit none
  INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
  !** Data structures
  TYPE::mpm_body
    CHARACTER(10):: name
    INTEGER:: A,Bi,emps,newnodes,nyp,newel,nmps,nn,slopeopt,                    &
            slopeel,slope1,tep
              
    REAL(iwp)::w1,s1,h1
    REAL(iwp)::Young,Poiss,frictfact
    INTEGER:: nex,ney,nels,locx,locy,ale,neq,np_types,nprops,dist_x,dist_y,     &
              nx1,nx2,ny1,ny2
    INTEGER,ALLOCATABLE:: b(:),etype(:),g(:),g_g(:,:),g_num(:,:),nf(:,:),       &
    num(:),MPPE(:),neighb(:,:),nodecont(:),newels(:),valuesg(:),boundnod(:),    &
    kdiag(:),kconst(:),nodecont_f(:),tied_nn(:,:),base_nn(:)
    !Variables to track material points
    INTEGER,ALLOCATABLE:: a_ele(:),c_ele(:),d_ele(:),b_ele(:),dp_ele(:),        &
    flag(:),GIMP_nodes(:,:),GIMP_node_mp(:,:),k_ele(:),elemmpoints(:,:)
    !Material properties
    REAL(iwp),ALLOCATABLE:: dee(:,:),deeinv(:,:),epsinv_acum(:),g_matrix(:),    &
    mpcp(:),prop(:,:)
    !Variables for each material point
    REAL(iwp),ALLOCATABLE:: accp(:,:),ins_acum(:,:),Devstress(:),eps_acum(:,:), &
      gm_coord(:,:),g_coord(:,:),ini_density(:),ini_volume(:),ins(:,:),         &
      mweights(:),m_volume(:),m_coord(:,:),m_stress(:,:), m_velocity(:,:),      &
      m_acc(:,:),mpoints(:,:),m_mass(:),m_stress_ini(:,:),accb(:,:),            &  
      m_stress_change(:,:),vccp(:,:),mpyield(:),mean_stress(:),lp_mp(:,:),      &
      m_stress_prev(:,:),ground_acc(:),mp_dens(:) 
    !Single body field
    REAL(iwp),ALLOCATABLE:: a_field(:,:),ddylds(:),diag(:),d1x1(:),d2x1(:),     &
      eps_m(:,:),eps_1(:,:),eps_2(:,:),fnorm(:),fcont(:),fdamp(:),Freact(:),    &
      f_ff(:),gravlo(:),kv(:),kp(:),kp_2(:),kinup_d2x1(:),kinup_d1x1(:),        &
      loads(:),mv(:),cv(:),m_mvp(:),m_mva(:),m_field(:),m_phi(:),normal(:,:),x1(:),   &
      temp_d1x1(:),temp_d2x1(:),tangent(:,:),v_field(:,:),vcm(:),vel_change(:), &
      x1_orig(:),x1_change(:),x1_ini(:),f_earth(:),kinup_Ground_d2x1(:),cdamp(:),mvkv(:)
    !Contact-specific variables
    REAL(iwp):: phi,gimptol
    LOGICAL::activate=.false.

    !Additional Variables
    INTEGER*8           :: pt(64)
    INTEGER,ALLOCATABLE :: iparm(:)
    INTEGER:: mpart,npart,skylength,nn_2D,ndof
    INTEGER,ALLOCATABLE :: ia(:),ja(:),ia_aux(:),ja_aux(:)
    REAL(iwp),ALLOCATABLE:: g_coord_aux(:,:),m_pore(:,:),ground_loads(:),m_stress_efe(:,:),m_dens(:)
    REAL(iwp),ALLOCATABLE:: c_fact(:),P_ext(:),Cc(:,:),mvis(:),km(:,:),         &
      x1_acum(:), P_matrix(:),kv_CSR(:),kv_CSR_aux(:),c_damp(:),c_force(:),     &
      residual(:), loads_ini(:),loads_end(:),c_matrix(:,:),d1p1(:),mf_matrix(:,:),&
      mf_force(:), loads_base(:),KGC(:,:),MMS(:,:),CCS(:,:),MOD_MTX(:,:)
    INTEGER,ALLOCATABLE::penpos(:),eq_pen(:),eq_pen_v(:),penpos_v(:)
    REAL(iwp),ALLOCATABLE::ddsddt(:),drplde(:),drpldt(:),stran(:),              &
      time(:),predef(:),dpred(:),props(:),drot(:,:),dfgrd0(:,:),dfgrd1(:,:),    &
      stressumat(:),epsumat(:),deeumat(:,:),statev(:,:)
    REAL(iwp),ALLOCATABLE::damping(:),fk(:)
  END TYPE
  
contains

  SUBROUTINE allocate_body(mbod)
    IMPLICIT NONE
    CLASS(mpm_body)::mbod
    ! Allocate kinematics state variables
    ALLOCATE(                                 &
      mbod%x1(0:mbod%neq),          &
      mbod%x1_orig(0:mbod%neq),     &
      mbod%x1_ini(0:mbod%neq),      &
      mbod%x1_acum(0:mbod%neq),     &
      mbod%x1_change(0:mbod%neq),   &
      mbod%d1x1(0:mbod%neq),        &
      mbod%d2x1(0:mbod%neq),        &
      mbod%temp_d2x1(0:mbod%neq),   &
      mbod%temp_d1x1(0:mbod%neq),   &
      mbod%kinup_d1x1(0:mbod%neq),  &
      mbod%kinup_d2x1(0:mbod%neq),  &                                    
      mbod%kinup_Ground_d2x1(0:mbod%neq))

    ! Allocate nodal forces and momentum
    ALLOCATE(                                 &  
      mbod%m_mvp(0:mbod%neq),       &
      mbod%m_mva(0:mbod%neq),       &
      mbod%gravlo(0:mbod%neq),      &
      mbod%ddylds(0:mbod%neq),      &
      mbod%vcm(0:mbod%neq),         &
      mbod%f_ff(0:mbod%neq),        &
      mbod%f_earth(0:mbod%neq),     &
      mbod%ground_loads(0:mbod%neq),&
      mbod%cdamp(0:mbod%neq),       &
      mbod%c_force(0:mbod%neq),     &
      mbod%mf_force(0:mbod%neq),    &
      mbod%residual(0:mbod%neq),    &
      mbod%loads_ini(0:mbod%neq),   &
      mbod%loads_end(0:mbod%neq),   &
      mbod%loads_base(0:mbod%neq),  &
      mbod%loads(0:mbod%neq))

    ! allocate unused kinematics and forces variables
    ALLOCATE(                       &
      mbod%d1p1(0:mbod%neq),        &
      mbod%P_ext(0:mbod%neq),       &
      mbod%c_damp(0:mbod%neq))

    ALLOCATE(                       &
      mbod%fcont(0:mbod%neq),       &
      mbod%fnorm(0:mbod%neq),       &
      mbod%fdamp(0:mbod%neq),       &
      mbod%Freact(0:mbod%neq),      &
      mbod%vel_change(0:mbod%neq),  &
      mbod%diag(0:mbod%neq),        &
      mbod%kdiag(mbod%neq),         &
      mbod%mvkv(0:mbod%neq))

    ! Penalty positions
    ALLOCATE(                       &
      mbod%penpos(mbod%neq),        &
      mbod%penpos_v(mbod%neq),      &
      mbod%eq_pen(mbod%neq),        &
      mbod%eq_pen_v(mbod%neq))
  END SUBROUTINE


  SUBROUTINE deallocate_body(mbod)
    IMPLICIT NONE
    CLASS(mpm_body)::mbod
    ! Allocate kinematics state variables
    DEALLOCATE(                                 &
      mbod%x1,          &
      mbod%x1_orig,     &
      mbod%x1_ini,      &
      mbod%x1_acum,     &
      mbod%x1_change,   &
      mbod%d1x1,        &
      mbod%d2x1,        &
      mbod%temp_d2x1,   &
      mbod%temp_d1x1,   &
      mbod%kinup_d1x1,  &
      mbod%kinup_d2x1,  &                                    
      mbod%kinup_Ground_d2x1)

    ! Allocate nodal forces and momentum
    DEALLOCATE(                                 &  
      mbod%m_mvp,       &
      mbod%m_mva,       &
      mbod%gravlo,      &
      mbod%ddylds,      &
      mbod%vcm,         &
      mbod%f_ff,        &
      mbod%f_earth,     &
      mbod%ground_loads,&
      mbod%cdamp,       &
      mbod%c_force,     &
      mbod%mf_force,    &
      mbod%residual,    &
      mbod%loads_ini,   &
      mbod%loads_end,   &
      mbod%loads_base,  &
      mbod%loads)

    ! allocate unused kinematics and forces variables
    DEALLOCATE(         &
      mbod%d1p1,        &
      mbod%P_ext,       &
      mbod%c_damp)

    DEALLOCATE(         &
      mbod%fcont,       &
      mbod%fnorm,       &
      mbod%fdamp,       &
      mbod%Freact,      &
      mbod%vel_change,  &
      mbod%diag,        &
      mbod%kdiag,       &
      mbod%mvkv)

    ! Penalty positions
    DEALLOCATE(         &
      mbod%penpos,      &
      mbod%penpos_v,    &
      mbod%eq_pen,      &
      mbod%eq_pen_v)
  END SUBROUTINE
  
  SUBROUTINE IO_POINT_VIZ(input,coord,a_ins,evpt,m_stress,m_stress_inc,acc,     &
    velocity,cohesion,devstress,meanstress,mpyield,directory,argv)
    !
    ! SUBROUTINE used to save visualise outputs to Paraview format
    ! https://dav.lbl.gov/archive/NERSC/Software/express/help6.1/help/reference/dvmac/UCD_Form.htm
    !
    IMPLICIT NONE
    INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
    REAL(iwp), PARAMETER :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp
    INTEGER, INTENT(IN) :: input
    REAL(iwp), INTENT(IN) :: coord(:,:)
    REAL(iwp), INTENT(IN) :: a_ins(:,:), evpt(:,:), m_stress(:,:),              &
      m_stress_inc(:,:), acc(:,:), velocity(:,:), cohesion(:), devstress(:),    &
      meanstress(:), mpyield(:)
    CHARACTER(*), OPTIONAL, INTENT(IN):: argv, directory
    
    CHARACTER(128) :: cnumber, s_directory="Results/", s_argv="Particles"
    INTEGER :: i, nmps, unit
    !--- File setup
    if (present(directory)) s_directory = directory
    if (present(argv)) s_argv = argv
    ! CALL EXECUTE_COMMAND_LINE("mkdir "//trim(s_directory))
    unit=1000
    WRITE(cnumber, '(i8.8)') input
    OPEN(unit,FILE=trim(s_directory)//trim(s_argv)//"_"//trim(cnumber)//'.inp')

    !--- Variables
    nmps = ubound(coord, 2)

    !--- File Comments
    write(unit,*) '#'
    write(unit,*) '# Simple AVS UCD File'
    write(unit,*) '#'

    !--- General Structure
    ! Number of nodes, number of cells, number of data items per node, 
    ! number of data items per cell, model number (always 0)
    write(unit,'(2i7.1,a)') nmps, nmps, ' 3 0 0' 

    !--- Node Coordinates
    do i=1,nmps ! coordinates
      write(unit,'(i5.1,3g13.5)') i, coord(:,i), zero
    end do
    
    !--- Cell Descriptions
    do i=1,nmps ! points
        write(unit,'(2i5.1,a,i5.1)') i, 0, '1 pt ', i
    end do

    !--- Node Based Data Descriptions
    write(unit,*) '10 2 4 4 4 2 2 1 1 1 1'
    write(unit,*) 'Displacement, meters'    ! ins_acum
    write(unit,*) 'Strains, meters'         ! evpt
    write(unit,*) 'Stress, kPa'             ! m_stress
    write(unit,*) 'Stress Increase, kPa'    ! m_stress_inc
    write(unit,*) 'Acceleration, m2/s'      ! acc
    write(unit,*) 'Velocity, m/s'           ! velocity
    write(unit,*) 'Cohesion, kPa'           ! cohesion
    write(unit,*) 'Dev-stress, no-unit'     ! devstress
    write(unit,*) 'Mean-stress, no-unit'    ! meanstress
    write(unit,*) 'Yield-value, no-unit'
    
    !--- Node Based Data
    DO i=1,nmps
      WRITE(unit,'(i5.1,22f13.3)') i, a_ins(:,i), evpt(:,i), m_stress(:,i),         &
        m_stress_inc(:,i), acc(:,i), velocity(:,i), cohesion(i), devstress(i),  &
        meanstress(i), mpyield(i)
    END DO

    CLOSE(unit)
  END SUBROUTINE IO_POINT_VIZ


  SUBROUTINE IO_PARAVIEW(input,node_type,coord,num,nf,kdiag,diag,loads,         &
    ddylds,gravlo,vcm,cdamp,f_earth,f_ff,x1,d1x1,d2x1,kv,mv,directory,argv)
    !
    ! Subroutine used to save visualization output of computational mesh
    ! https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html
    !
    IMPLICIT NONE
    INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
    REAL(iwp), PARAMETER :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp
    INTEGER, INTENT(IN)    :: input, node_type
    INTEGER, INTENT(IN)    :: num(:,:)
    REAL(iwp), INTENT(IN)  :: coord(:,:)
    CHARACTER(*), OPTIONAL :: directory, argv
    INTEGER, OPTIONAL, INTENT(IN)   :: nf(:,:), kdiag(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: loads(:), ddylds(:), gravlo(:), vcm(:), cdamp(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: diag(:), f_earth(:), f_ff(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: kv(:), mv(:)
    REAL(iwp), OPTIONAL, INTENT(IN) :: x1(:), d1x1(:), d2x1(:)

    CHARACTER(128) :: cnumber, s_directory="Results/", s_argv="Paraview"
    INTEGER :: i, iel, nels, nn, unit

    !--- File setup
    if (present(directory)) s_directory = directory
    if (present(argv)) s_argv = argv
    ! CALL EXECUTE_COMMAND_LINE("mkdir "//trim(s_directory))
    unit=1000
    WRITE(cnumber, '(i8.8)') input
    OPEN(unit,FILE=trim(s_directory)//trim(s_argv)//"_"//trim(cnumber)//'.vtk')

    !--- Variables
    nn = ubound(coord, 2)
    nels = ubound(num, 2)
    
    !--- File Description Parts (cf. VTK File Format)
    WRITE(unit,'(a)')'# vtk DataFile Version 3.0' ! Header
    WRITE(unit,'(a)')"vtk output"                 ! Title
    WRITE(unit,'(a)')"ASCII"                      ! Data type
    WRITE(unit,'(a)')""                           !
    WRITE(unit,'(a)')"DATASET UNSTRUCTURED_GRID"  ! Geometry/Topology

    !--- Datasets: Geometry
    ! Node Coordinates
    WRITE(unit,'(1A6,1I5,1A7)')"POINTS", nn, "float"
    DO i=1, nn
        WRITE(unit,'(3f9.4)') coord(:,i), zero
    END DO
    ! Cell Constructs
    WRITE(unit,'(a)')""
    SELECT CASE(node_type)
      CASE(4)
        WRITE(unit,'(1A6,2I6)')"CELLS", nels, nels*(1+4)
        DO iel=1, nels
            WRITE(unit,'(5I5)') node_type,  &
              num(1,iel)-1, num(4,iel)-1, num(3,iel)-1, num(2,iel)-1
        END DO

        WRITE(unit,'(a)')""
        WRITE(unit,'(1A10,1I5)')"CELL_TYPES", nels
        DO iel = 1 , nels
            WRITE(unit,*)"9"
        END DO

      CASE DEFAULT
        WRITE(*,*)"wrong number of nodes input in paraview"
    END SELECT

    IF (.not. present(nf)) RETURN
    !--- Datasets: Node Values
    WRITE(unit,'(a)')""
    WRITE(unit,'(1A10,1I5)')"POINT_DATA", nn
    NODAL_MASS: IF (present(diag)) THEN
      WRITE(unit,'(a)')"vectors mass float "
      DO i=1, nn
          IF(nf(1,i)==0 .and. nf(2,i)==0) WRITE(unit,'(3f15.4)')zero, zero, zero
          IF(nf(1,i)>0  .and. nf(2,i)==0) WRITE(unit,'(3f15.4)')diag(nf(1,i)+1), zero, zero
          IF(nf(1,i)==0 .and. nf(2,i)>0 ) WRITE(unit,'(3f15.4)')zero, diag(nf(2,i)+1), zero
          IF(nf(1,i)>0  .and. nf(2,i)>0 ) WRITE(unit,'(3f15.4)')diag(nf(1,i)+1), diag(nf(2,i)+1), zero
      END DO
    END IF NODAL_MASS

    F_TOTAL: IF (present(loads)) THEN
      WRITE(unit,'(a)')"vectors F_Total float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')loads(nf(:,i)+1), zero
      END DO
    END IF F_TOTAL

    F_INTERNAL: IF (present(ddylds)) THEN
      WRITE(unit,'(a)')"vectors F_Internal float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')ddylds(nf(:,i)+1), zero
      END DO
    END IF F_INTERNAL

    F_EXTERNAL: IF (present(gravlo)) THEN
      WRITE(unit,'(a)')"vectors F_External float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')gravlo(nf(:,i)+1), zero
      END DO
    END IF F_EXTERNAL

    F_KINETIC: IF (present(vcm)) THEN
      WRITE(unit,'(a)')"vectors F_Kinetic float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')vcm(nf(:,i)+1), zero
      END DO
    END IF F_KINETIC
    
    F_DAMPING: IF (present(cdamp)) THEN
      WRITE(unit,'(a)')"vectors F_Damping float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')cdamp(nf(:,i)+1), zero
      END DO
    END IF F_DAMPING

    F_EARTHQUAKE: IF (present(f_earth)) THEN
      WRITE(unit,'(a)')"vectors f_earth float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')f_earth(nf(:,i)+1), zero
      END DO
    END IF F_EARTHQUAKE

    F_FREEFIELD: IF (present(f_ff)) THEN
      WRITE(unit,'(a)')"vectors f_freefield float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')f_ff(nf(:,i)+1), zero
      END DO
    END IF F_FREEFIELD
    
    NODAL_DISPLACEMENT: IF (present(x1)) THEN
      WRITE(unit,'(a)')"vectors Displacement float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')x1(nf(:,i)+1), zero
      END DO
    END IF NODAL_DISPLACEMENT

    NODAL_VELOCITY: IF (present(d1x1)) THEN
      WRITE(unit,'(a)')"vectors Velocity float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')d1x1(nf(:,i)+1), zero
      END DO
    END IF NODAL_VELOCITY
    
    NODAL_ACC: if (present(d2x1)) THEN
      WRITE(unit,'(a)')"vectors Acceleration float "
      DO i=1, nn
          WRITE(unit,'(3f15.6)')d2x1(nf(:,i)+1), zero
      END DO 
    END IF NODAL_ACC

    STIFFNESS_MATRIX: if (present(kv) .and. present(kdiag)) THEN
      WRITE(unit,'(a)')"vectors KM_Stiffness float "
      DO i=1, nn
          IF(nf(1,i)==0 .and. nf(2,i)==0) WRITE(unit,'(3f15.4)') zero, zero, zero
          IF(nf(1,i)>0  .and. nf(2,i)==0) WRITE(unit,'(3f15.4)') kv(kdiag(nf(1,i))), zero, zero
          IF(nf(1,i)==0 .and. nf(2,i)>0 ) WRITE(unit,'(3f15.4)') zero, kv(kdiag(nf(2,i))), zero
          IF(nf(1,i)>0  .and. nf(2,i)>0 ) WRITE(unit,'(3f15.4)') kv(kdiag(nf(1,i))), kv(kdiag(nf(2,i))), zero
      END DO
    END IF STIFFNESS_MATRIX

    MASS_MATRIX: if (present(mv) .and. present(kdiag)) THEN
      WRITE(unit,'(a)')"vectors Mv_Mass float "
      DO i=1, nn
          IF(nf(1,i)==0 .and. nf(2,i)==0) WRITE(unit,'(3f15.4)') zero, zero, zero
          IF(nf(1,i)>0  .and. nf(2,i)==0) WRITE(unit,'(3f15.4)') mv(kdiag(nf(1,i))), zero, zero
          IF(nf(1,i)==0 .and. nf(2,i)>0 ) WRITE(unit,'(3f15.4)') zero, mv(kdiag(nf(2,i))), zero
          IF(nf(1,i)>0  .and. nf(2,i)>0 ) WRITE(unit,'(3f15.4)') mv(kdiag(nf(1,i))),mv(kdiag(nf(2,i))), zero
      END DO
    END IF MASS_MATRIX

    ! Close File
    CLOSE(unit)
  END SUBROUTINE IO_PARAVIEW
  
end module utilities