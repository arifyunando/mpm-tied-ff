MODULE MPM_CORE
  USE CLASS_PARTICLE
  USE CLASS_MATERIAL
  USE CLASS_MESH
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp
  
  TYPE::mpm_body
    CLASS(particles),ALLOCATABLE::particles
    CLASS(material_model),ALLOCATABLE::material
    CLASS(mesh),ALLOCATABLE::mesh

    INTEGER,ALLOCATABLE::a_ele(:) ! Element index each MP resides 
    INTEGER,ALLOCATABLE::c_ele(:) ! Number of MPs inside a given element
    LOGICAL,ALLOCATABLE::d_ele(:) ! Active element (1:active; 0:deactive)
    INTEGER,ALLOCATABLE::k_ele(:) ! Accumulation of c_ele
    INTEGER::n_active_ele

    CONTAINS

    PROCEDURE :: SET_MESH => p_SET_MESH
    PROCEDURE :: LOCATE_PARTICLES => p_LOCATE_PARTICLES
    PROCEDURE :: GET_PARTICLE_SUPPORT_SIZE => p_GET_PARTICLE_SUPPORT_SIZE
    PROCEDURE :: GET_SUPPORT_ELEMENTS => p_GET_SUPPORT_ELEMENTS
    PROCEDURE :: GET_SUPPORT_NODES => p_GET_SUPPORT_NODES
    PROCEDURE :: SET_MATERIAL => p_SET_MATERIAL
    ! PROCEDURE :: CONSTRUCT_DMMPM_STIFFNESS => MPMCORE_CONSTRUCT_DMMPM_STIFFNESS
  END TYPE

  CONTAINS

!****************************** PUBLIC FUNCTIONS ******************************!
!                                                                              !

  SUBROUTINE p_SET_MESH(this, grid)
    !
    ! Assign a material body onto a background grid
    !
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    CLASS(mesh),INTENT(IN)::grid
    this%mesh=grid
    CALL m_INITIATE_MPM_BODY(this)
  END SUBROUTINE p_SET_MESH


  SUBROUTINE p_LOCATE_PARTICLES(this)
    !
    ! Look for all the elements that are affected by particles and mark it
    ! in a_ele, c_ele, and d_ele.
    !
    USE FUNCTIONS
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    INTEGER,ALLOCATABLE::v_neighbour_id(:)
    REAL(iwp),ALLOCATABLE::mp_coord(:)
    INTEGER::i,j,iel,iel_neigh,n_neighbour_id
    SELECT CASE (this%mesh%mesh_ind)
    CASE(4) ! Quadrilateral Cartesian 4 Nodes
      this%n_active_ele=0
      this%c_ele=0
      this%d_ele=.false.
      
      MPS: DO i=1,this%particles%nmps
        ! check particle has been located before
        iel=this%a_ele(i)
        IF (iel > 0) THEN
          ! check particle in previous known cell
          IF ( POINT_IN_CARTESIAN_CELL(this,iel,mp_coord) ) THEN
            this%a_ele(i) = iel
            this%c_ele(iel)=this%c_ele(iel)+1
            this%d_ele(iel)=.true.
            CYCLE MPS
          END IF

          ! check particle in neighbouring cell
          v_neighbour_id = this%mesh%v_neighbour_id(:,iel)
          n_neighbour_id = this%mesh%n_neighbour_id(iel)
          DO j=1,n_neighbour_id
            iel_neigh=v_neighbour_id(j)
            IF ( POINT_IN_CARTESIAN_CELL(this,iel_neigh,mp_coord) ) THEN
              this%a_ele(i) = iel
              this%c_ele(iel)=this%c_ele(iel)+1
              this%d_ele(iel)=.true.
              CYCLE MPS
            END IF
          END DO
        END IF

        ! check particle in all cell
        DO iel=1,this%mesh%nels
          mp_coord=this%particles%gm_coord(:,i)
          IF ( POINT_IN_CARTESIAN_CELL(this,iel,mp_coord) ) THEN
            this%a_ele(i)=iel
            this%c_ele(iel)=this%c_ele(iel)+1
            this%d_ele(iel)=.true.
          END IF
        END DO
      END DO MPS

    CASE DEFAULT
      WRITE(*, *) 'indicated mesh type not found. Aborting program'
      PAUSE
      CALL EXIT(-1)
    END SELECT
  END SUBROUTINE p_LOCATE_PARTICLES


  SUBROUTINE p_GET_PARTICLE_SUPPORT_SIZE(this)
    !
    ! Calculate lp for each material points
    ! Note: Must be called after FLAG_ELEMENTS()
    !       to make sure c_ele exists
    !
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    INTEGER::i
    REAL(iwp)::nmps_in_cell
    DO i=1,this%particles%nmps
      nmps_in_cell=SQRT(REAL(this%c_ele(i)))
      this%particles%lp(:,i)=this%mesh%cellsize(:,this%a_ele(i))/nmps_in_cell
    END DO
  END SUBROUTINE p_GET_PARTICLE_SUPPORT_SIZE


  SUBROUTINE p_SET_MATERIAL(this,input_json)
    !
    ! Load Material to body according to particles.generator.material_id
    !
    USE JSON_MODULE
    IMPLICIT NONE
    CLASS(mpm_body), INTENT(INOUT) :: this
    TYPE(json_file), INTENT(INOUT) :: input_json
    LOGICAL :: found=.true.
    INTEGER :: index
    CHARACTER(1) :: i_char
    ! obtain material id
    WRITE(i_char, '(I1)') this%particles%id
    CALL input_json%GET('particles('//i_char//').generator.material_id', index)
    ! check material existance. load if material available
    WRITE(i_char, '(I1)') index
    CALL input_json%GET('materials('//i_char//').id', index, found)
    IF (found) THEN
      ALLOCATE(this%material)
      CALL this%material%LOAD(index, input_json)
    END IF
  END SUBROUTINE p_SET_MATERIAL


  ! SUBROUTINE MPMCORE_CONSTRUCT_DMMPM_STIFFNESS(this,nip,nst)
  !   !
  !   ! Construct Stiffness Matrix According to DM-MPM Formulation
  !   ! cf. Gonzales Acosta, 2021. page 54 Eq.3.9
  !   !
  !   USE FUNCTIONS
  !   IMPLICIT NONE
  !   CLASS(mpm_body),INTENT(INOUT)::this
  !   INTEGER,INTENT(IN)::nip,nst
  !   ! Local Variable
  !   INTEGER::i,i_mp,i_element,iel
  !   INTEGER,ALLOCATABLE::num(:),coord(:,:),g(:),node_num(:,:)
  !   REAL(iwp)::w_average,factor_dee,det
  !   REAL(iwp),ALLOCATABLE::func_gauss(:),deriv_gauss(:,:)
  !   REAL(iwp),ALLOCATABLE::func_gimp(:),deriv_gimp(:,:),deriv(:,:)
  !   REAL(iwp),ALLOCATABLE::scaled_dee(:,:),jac(:,:),bee(:,:),km(:,:)
  !   REAL(iwp),ALLOCATABLE::weights(:),points(:,:)
  !   CHARACTER(13)::element='quadrilateral'
  !   ! Local Memory Allocations
  !   ALLOCATE(                                     &
  !     num(this%mesh%nod),                         &
  !     node_num(this%mesh%nod,1),                  &
  !     coord(this%mesh%nod, this%mesh%ndim),       &
  !     g(this%mesh%nodof),                         &
  !     func_gauss(this%mesh%nod),                  &
  !     deriv_gauss(this%mesh%ndim,this%mesh%nod),  &
  !     func_gimp(this%mesh%nod),                   &
  !     deriv_gimp(this%mesh%ndim,this%mesh%nod),   &
  !     deriv(this%mesh%ndim,this%mesh%nod),        &
  !     scaled_dee(nst,nst),                        &
  !     jac(this%mesh%ndim,this%mesh%ndim),         &
  !     bee(nst,this%mesh%ndof),                    &
  !     km(this%mesh%ndof,this%mesh%ndof),          &
  !     weights(nip),                               &
  !     points(nip,this%mesh%ndim)                  &
  !   )
  !   ! Loop over all MPs
  !   DO i_mp=1,this%nmps
  !     this%mv=zero
  !     DO i_element=1,4 ! 4 is th maximum number of supporting element for a particle
  !       iel = this%member_elements(i_mp,i_element)
  !       GET_STIFFNESS: IF (iel > 0) THEN
  !         num = this%mesh%g_num(:,iel)
  !         coord = TRANSPOSE(this%mesh%g_coord(:,num))
  !         g = this%g_g(:,iel)

  !         km = zero   
  !         w_average = this%mweights(i_mp)
          
  !         ! Double Mapping
  !         node_num(:,1) = num
  !         CALL GIMP_GET_SHAPE_FUNCTION_AND_DERIVATIVES(                        &
  !           func=func_gimp,                                                    &
  !           deriv=deriv_gimp,                                                  &
  !           s=i,                                                               &
  !           iel=iel,                                                           &
  !           lm_coord=this%lm_coord,                                            &
  !           lp=this%lp,                                                        &
  !           cellsize=this%mesh%cellsize,                                       &
  !           g_coord=this%mesh%g_coord,                                         &
  !           gm_coord=this%gm_coord,                                            &
  !           support_nodes=node_num,                                            &
  !           a_ele=this%a_ele                                                   &
  !         )
  !         CALL SAMPLE_GAUSS_IP(element, points, weights)
  !         GAUSS_INTEGRATION: DO i=1,nip
  !           CALL SHAPE_FUN(func_gauss, points, i)
  !           CALL SHAPE_DER(deriv_gauss, points, i)
  !           factor_dee = SUM(func_gimp*func_gauss*w_average)
  !           scaled_dee = factor_dee*this%dee
  !           jac = MATMUL(deriv_gauss, coord)
  !           det = DETERMINANT(jac)
  !           CALL INVERT(jac)
  !           deriv = MATMUL(jac, deriv_gauss)
  !           CALL BEEMAT(bee, deriv)

  !           IF(this%c_ele(iel) < 1) THEN
  !             ! adjust double mapping weight
  !             w_average=1.5_iwp 
  !           END IF

  !           km = km + MATMUL(MATMUL(TRANSPOSE(bee),scaled_dee),bee)*det*weights(i)
  !         END DO GAUSS_INTEGRATION
  !         CALL FSPARV(this%kv, km, g, this%kdiag)
  !       END IF GET_STIFFNESS
  !     END DO
  !   END DO
  ! END SUBROUTINE MPMCORE_CONSTRUCT_DMMPM_STIFFNESS


  ! SUBROUTINE MPMCORE_CONSTRUCT_GIMP_MASS_MATRIX(this)
  !   !
  !   !
  !   !
  !   IMPLICIT NONE
  !   CLASS(mpm_body), INTENT(INOUT) :: this

  !   ! For each Material Points:
  !   !   Calculate GIMP Shape Function
  !   !   Determine Streering Vector for building the Mass Matrix 
  !   !   NB: The steering vector is contructed based on the supporting domain of
  !   !       of the particles
  !   !   For each supporting nodes:
  !   !     Construct Local Mass Matrices
  !   !     Sum over all supporting nodes
  !   !   Add local mass matrix to global mass matrix (fsparv)
  
  ! END SUBROUTINE MPMCORE_CONSTRUCT_GIMP_MASS_MATRIX


  SUBROUTINE p_GET_SUPPORT_ELEMENTS(this)
    !
    ! Activate element by assigning active element variable to body
    !
    USE FUNCTIONS
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    CALL GIMP_GET_SUPPORT_ELEMENTS(                   &
      member_elements=this%particles%member_elements, &
      g_coord=this%mesh%g_coord,                      &
      gm_coord=this%particles%gm_coord,               &
      c_ele=this%c_ele,                               &
      nf=this%mesh%nf,                                &
      g_num=this%mesh%g_num,                          &
      nels=this%mesh%nels,                            &
      lp=this%particles%lp                            &
    )
  END SUBROUTINE p_GET_SUPPORT_ELEMENTS


  SUBROUTINE p_GET_SUPPORT_NODES(this,gimptol)
    !     
    ! Subroutine to save al the nodes (number of the nodes) inside the 
    ! support domain for each material point with GIMP 
    ! node numbering should be in y direction
    !
    USE FUNCTIONS
    IMPLICIT NONE
    INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
    ! dummy variables
    CLASS(mpm_body),INTENT(INOUT)::this
    REAL(iwp),OPTIONAL,INTENT(IN)::gimptol
    ! local variables
    REAL(iwp)::minx,maxx,miny,maxy,span_x,span_y,tol=1.0e-6
    REAL(iwp),ALLOCATABLE::g_coord(:,:),gm_coord(:,:)
    INTEGER,ALLOCATABLE::nodes(:),neigh(:)
    INTEGER,ALLOCATABLE::gimp_nodes_tmp(:,:),gimp_nodes(:,:),gimp_nodes_count(:)
    INTEGER::s,i,count,iel,nn,nmps,neigh_count
    
    IF(present(gimptol)) tol=gimptol
    nn=this%mesh%nn
    nmps=this%particles%nmps
    g_coord=this%mesh%g_coord
    gm_coord=this%particles%gm_coord

    ALLOCATE(gimp_nodes_count(nmps),gimp_nodes_tmp(30,nmps))
    DO s=1,nmps
      iel=this%a_ele(s)
      span_x=this%particles%lp(1,s)+this%mesh%cellsize(1,iel)
      span_y=this%particles%lp(2,s)+this%mesh%cellsize(2,iel)
      
      minx=gm_coord(1,s)-span_x+tol ! Left boundarie of the support domain
      maxx=gm_coord(2,s)+span_y-tol ! Upper boundarie of the support domain
      miny=gm_coord(1,s)+span_x+tol ! Right boundarie of the support domain
      maxy=gm_coord(2,s)-span_y-tol ! Lower boundarie of the support domain
      
      neigh=this%mesh%v_neighbour_id(:,iel)
      neigh_count=this%mesh%n_neighbour_id(iel)
      CALL GET_NEIGHBOUR_NODES(nodes,neigh,neigh_count,this%mesh%g_num)

      count=0
      DO i=1,size(nodes)
        IF(g_coord(1,nodes(i))>=minx.and.g_coord(1,nodes(i))<=maxx)THEN 
        IF(g_coord(2,nodes(i))<=maxy.and.g_coord(2,nodes(i))>=miny)THEN  
          count=count+1 
          gimp_nodes_tmp(count,s)=nodes(i)
        END IF; END IF
      END DO
      gimp_nodes_count(s)=count
    END DO
    
    count=MAXVAL(gimp_nodes_count)
    ALLOCATE(gimp_nodes(count,nmps))
    gimp_nodes(:,:)=gimp_nodes_tmp(:count,:)
    CALL MOVE_ALLOC(gimp_nodes, this%particles%v_support_nodes)
    CALL MOVE_ALLOC(gimp_nodes_count, this%particles%n_support_nodes)
  END SUBROUTINE p_GET_SUPPORT_NODES

!                                                                              !
!****************************** PRIVATE FUNCTIONS *****************************!
!                                                                              !
  
  SUBROUTINE m_INITIATE_MPM_BODY(this)
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    !--- Check initiation
    IF (.not. (this%mesh%is_initiated .or. this%particles%is_initiated)) THEN
      WRITE(*,*) "Mesh initiatiation is ", this%mesh%is_initiated 
      WRITE(*,*) "Particles initiatiation is ", this%particles%is_initiated
      WRITE(*,*) "MPM Body initiation failed. Make sure all members are properly initiated"
      PAUSE
      CALL EXIT(-1)
    END IF

    !--- Allocate Tracking Variables
    ALLOCATE(                            &
      this%a_ele(this%particles%nmps),   &
      this%c_ele(this%mesh%nels),        &
      this%d_ele(this%mesh%nels),        &
      this%k_ele(0:this%mesh%nels)       &
    )
    this%a_ele = 0 ! element id for each material point
    this%c_ele = 0 ! total of MP inside each element
    this%k_ele = 0 ! total of MP in the domain
    this%d_ele = 0 ! activated element array (1 active/0 deactive)
  END SUBROUTINE m_INITIATE_MPM_BODY


  LOGICAL FUNCTION POINT_IN_CARTESIAN_CELL(this,iel,point_coord)
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    INTEGER,INTENT(IN)::iel
    REAL(iwp),ALLOCATABLE,INTENT(IN)::point_coord(:)
    ! Local variable
    INTEGER,ALLOCATABLE::num(:)
    REAL(iwp),ALLOCATABLE::coord(:,:)
    REAL(iwp)::el_coordx_max,el_coordx_min,el_coordy_max,el_coordy_min
    ! determine boundary
    num=this%mesh%g_num(:,iel)
    coord=TRANSPOSE(this%mesh%g_coord(:,num))
    el_coordx_max=MAXVAL(coord(:,1))
    el_coordx_min=MINVAL(coord(:,1))
    el_coordy_max=MAXVAL(coord(:,2))
    el_coordy_min=MINVAL(coord(:,2))
    ! check point location
    POINT_IN_CARTESIAN_CELL = .false.
    IF(point_coord(1)>el_coordx_min)THEN; IF(point_coord(1)<el_coordx_max)THEN
    IF(point_coord(2)>el_coordy_min)THEN; IF(point_coord(2)<el_coordy_max)THEN
      POINT_IN_CARTESIAN_CELL = .true.
    END IF; END IF; END IF; END IF
  END FUNCTION

END MODULE MPM_CORE

