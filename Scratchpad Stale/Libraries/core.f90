MODULE MPM_CORE
  USE CLASS_PARTICLE
  USE CLASS_MATERIAL
  USE CLASS_MESH
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp
  
  TYPE, EXTENDS(particles)::mpm_body
    !
    ! MPM Particle Base Class
    !
    CLASS(material_model), ALLOCATABLE :: material
    CLASS(mpm_grid), ALLOCATABLE :: mesh

    CONTAINS

    PROCEDURE :: SET_MESH => MPMCORE_SET_MESH
    PROCEDURE :: SET_MATERIAL => MPMCORE_SET_MATERIAL
    PROCEDURE :: CONSTRUCT_DMMPM_STIFFNESS => MPMCORE_CONSTRUCT_DMMPM_STIFFNESS
    PROCEDURE :: FLAG_ELEMENTS => MPMCORE_FLAG_ELEMENTS
    PROCEDURE :: GET_SUPPORT_ELEMENTS => MPMCORE_GET_SUPPORT_ELEMENTS
    PROCEDURE :: GET_PARTICLE_SUPPORT_SIZE => MPMCORE_GET_PARTICLE_SUPPORT_SIZE
  END TYPE

  CONTAINS

  SUBROUTINE MPMCORE_SET_MESH(this, grid)
    !
    ! Assign a material body onto a background grid
    !
    IMPLICIT NONE
    CLASS(mpm_body), INTENT(INOUT) :: this
    CLASS(mpm_grid), INTENT(INOUT) :: grid

    this%mesh = grid ! Connect grid to body
    ALLOCATE(                                                                   &
      this%c_ele(1:grid%nels),                                                  &
      this%d_ele(1:grid%nels),                                                  &
      this%k_ele(0:grid%nels)                                                   &
    )

    this%a_ele = 0 ! element id for each material point
    this%c_ele = 0 ! total of MP inside each element
    this%k_ele = 0 ! total of MP in the domain
    this%d_ele = 0 ! activated element array (1 active/0 deactive)
  END SUBROUTINE MPMCORE_SET_MESH

  SUBROUTINE MPMCORE_GET_PARTICLE_SUPPORT_SIZE(this)
    !
    ! Calculate lp for each material points
    ! Note: Must be called after FLAG_ELEMENTS()
    !       to make sure c_ele exists
    !
    IMPLICIT NONE
    CLASS(mpm_body), INTENT(INOUT) :: this
    INTEGER::i
    DO i=1,this%nmps
      this%lp(:,i) = this%mesh%cellsize(:,this%a_ele(i)) / INT(SQRT(REAL(this%c_ele(i))))
    END DO
  END SUBROUTINE MPMCORE_GET_PARTICLE_SUPPORT_SIZE

  SUBROUTINE MPMCORE_SET_MATERIAL(this, input_json)
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
    WRITE(i_char, '(I1)') this%id
    CALL input_json%GET('particles('//i_char//').generator.material_id', index)
    ! check material existance. load if material available
    WRITE(i_char, '(I1)') index
    CALL input_json%GET('materials('//i_char//').id', index, found)
    IF (found) THEN
      ALLOCATE(this%material)
      CALL this%material%LOAD(index, input_json)
    END IF
  END SUBROUTINE MPMCORE_SET_MATERIAL

  ! TODOs
  SUBROUTINE MPMCORE_CONSTRUCT_DMMPM_STIFFNESS(this,nip,nst)
    !
    ! Construct Stiffness Matrix According to DM-MPM Formulation
    ! cf. Gonzales Acosta, 2021. page 54 Eq.3.9
    !
    USE FUNCTIONS
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    INTEGER,INTENT(IN)::nip,nst
    ! Local Variable
    INTEGER::i,i_mp,i_element,iel
    INTEGER,ALLOCATABLE::num(:),coord(:,:),g(:),node_num(:,:)
    REAL(iwp)::w_average,factor_dee,det
    REAL(iwp),ALLOCATABLE::func_gauss(:),deriv_gauss(:,:)
    REAL(iwp),ALLOCATABLE::func_gimp(:),deriv_gimp(:,:),deriv(:,:)
    REAL(iwp),ALLOCATABLE::scaled_dee(:,:),jac(:,:),bee(:,:),km(:,:)
    REAL(iwp),ALLOCATABLE::weights(:),points(:,:)
    CHARACTER(13)::element='quadrilateral'
    ! Local Memory Allocations
    ALLOCATE(                                     &
      num(this%mesh%nod),                         &
      node_num(this%mesh%nod,1),                  &
      coord(this%mesh%nod, this%mesh%ndim),       &
      g(this%mesh%nodof),                         &
      func_gauss(this%mesh%nod),                  &
      deriv_gauss(this%mesh%ndim,this%mesh%nod),  &
      func_gimp(this%mesh%nod),                   &
      deriv_gimp(this%mesh%ndim,this%mesh%nod),   &
      deriv(this%mesh%ndim,this%mesh%nod),        &
      scaled_dee(nst,nst),                        &
      jac(this%mesh%ndim,this%mesh%ndim),         &
      bee(nst,this%mesh%ndof),                    &
      km(this%mesh%ndof,this%mesh%ndof),          &
      weights(nip),                               &
      points(nip,this%mesh%ndim)                  &
    )
    ! Loop over all MPs
    DO i_mp=1,this%nmps
      this%mv=zero
      DO i_element=1,4 ! 4 is th maximum number of supporting element for a particle
        iel = this%member_elements(i_mp,i_element)
        GET_STIFFNESS: IF (iel > 0) THEN
          num = this%mesh%num(:,iel)
          coord = TRANSPOSE(this%mesh%g_coord(:,num))
          g = this%g_g(:,iel)

          km = zero   
          w_average = this%mweights(i_mp)
          
          ! Double Mapping
          node_num(:,1) = num
          CALL GIMP_GET_SHAPE_FUNCTION_AND_DERIVATIVES(                        &
            func=func_gimp,                                                    &
            deriv=deriv_gimp,                                                  &
            s=i,                                                               &
            iel=iel,                                                           &
            lm_coord=this%lm_coord,                                            &
            lp=this%lp,                                                        &
            cellsize=this%mesh%cellsize,                                       &
            g_coord=this%mesh%g_coord,                                         &
            gm_coord=this%gm_coord,                                            &
            support_nodes=node_num,                                            &
            a_ele=this%a_ele                                                   &
          )
          CALL SAMPLE_GAUSS_IP(element, points, weights)
          GAUSS_INTEGRATION: DO i=1,nip
            CALL SHAPE_FUN(func_gauss, points, i)
            CALL SHAPE_DER(deriv_gauss, points, i)
            factor_dee = SUM(func_gimp*func_gauss*w_average)
            scaled_dee = factor_dee*this%dee
            jac = MATMUL(deriv_gauss, coord)
            det = DETERMINANT(jac)
            CALL INVERT(jac)
            deriv = MATMUL(jac, deriv_gauss)
            CALL BEEMAT(bee, deriv)

            IF(this%c_ele(iel) < 1) THEN
              ! adjust double mapping weight
              w_average=1.5_iwp 
            END IF

            km = km + MATMUL(MATMUL(TRANSPOSE(bee),scaled_dee),bee)*det*weights(i)
          END DO GAUSS_INTEGRATION
          CALL FSPARV(this%kv, km, g, this%kdiag)
        END IF GET_STIFFNESS
      END DO
    END DO
  END SUBROUTINE MPMCORE_CONSTRUCT_DMMPM_STIFFNESS


  SUBROUTINE MPMCORE_CONSTRUCT_GIMP_MASS_MATRIX(this)
    !
    !
    !
    IMPLICIT NONE
    CLASS(mpm_body), INTENT(INOUT) :: this

    ! For each Material Points:
    !   Calculate GIMP Shape Function
    !   Determine Streering Vector for building the Mass Matrix 
    !   NB: The steering vector is contructed based on the supporting domain of
    !       of the particles
    !   For each supporting nodes:
    !     Construct Local Mass Matrices
    !     Sum over all supporting nodes
    !   Add local mass matrix to global mass matrix (fsparv)
  
  END SUBROUTINE MPMCORE_CONSTRUCT_GIMP_MASS_MATRIX


  SUBROUTINE MPMCORE_FLAG_ELEMENTS(this)
      !
      ! Look for all the elements that are affected by particles and mark it
      ! in a_ele, c_ele, and d_ele. Also, determine the local coordinates
      ! in lp_coord
      !
      USE FUNCTIONS
      IMPLICIT NONE
      CLASS(mpm_body),INTENT(INOUT)::this
      INTEGER,ALLOCATABLE::num(:)
      REAL(iwp)::mp_coordx,mp_coordy
      REAL(iwp)::el_coordx_max,el_coordx_min,el_coordy_max,el_coordy_min
      REAL(iwp),ALLOCATABLE::sp_coord(:,:),lp_coord(:,:)
      REAL(iwp),ALLOCATABLE::coord(:,:)
      INTEGER::i,iel,count_mp
      ALLOCATE(sp_coord(this%mesh%ndim,1),lp_coord(1,this%mesh%ndim))
      SELECT CASE (this%mesh%mesh_ind)
      CASE(4) ! Quadrilateral Cartesian 4 Element
        ALLOCATE(num(this%mesh%nn))
        this%n_active_ele=1
        this%d_ele=0
        this%c_ele=0
        count_mp=1
        DO iel=1,this%mesh%nels
          ! get element coordinates
          num=this%mesh%num(:,iel)
          coord=TRANSPOSE(this%mesh%g_coord(:,num))
          el_coordx_max=MAXVAL(coord(:,1))
          el_coordx_min=MINVAL(coord(:,1))
          el_coordy_max=MAXVAL(coord(:,2))
          el_coordy_min=MINVAL(coord(:,2))
          
          ! check for any MPs inside it
          DO i=1,this%nmps
            mp_coordx=this%gm_coord(1,i)
            mp_coordy=this%gm_coord(2,i)
            ! Check whether paticle 'i' is in element 'iel'
            IF(mp_coordx>el_coordx_min)THEN; IF(mp_coordx<el_coordx_max)THEN
            IF(mp_coordy>el_coordy_min)THEN; IF(mp_coordy<el_coordy_max)THEN
              this%a_ele(i)=iel
              this%c_ele(iel)=this%c_ele(iel)+1
              this%d_ele(iel)=1
              
              ! calculate local coordinate of the mp
              sp_coord(:,1)=this%gm_coord(:,i)
              CALL LOCATE_MP_LOCAL(lp_coord,coord,sp_coord)
              this%lm_coord(:,i)=lp_coord(1,:)
              
              ! count mp for early stop
              count_mp = count_mp + 1
            END IF; END IF; END IF; END IF
          END DO
          IF (count_mp>this%nmps) EXIT ! early stop loop if all mps are found
        END DO
      CASE DEFAULT
        WRITE(*, *) 'indicated mesh type not found. Aborting program'
        PAUSE
        CALL EXIT(1)
      END SELECT
  END SUBROUTINE MPMCORE_FLAG_ELEMENTS


  SUBROUTINE MPMCORE_GET_SUPPORT_ELEMENTS(this)
    !
    ! Activate element by assigning active element variable to body
    !
    USE FUNCTIONS
    IMPLICIT NONE
    CLASS(mpm_body),INTENT(INOUT)::this
    CALL GIMP_GET_SUPPORT_ELEMENTS(         &
      member_elements=this%member_elements, &
      g_coord=this%mesh%g_coord,            &
      gm_coord=this%gm_coord,               &
      c_ele=this%c_ele,                     &
      nf=this%mesh%nf,                      &
      g_num=this%mesh%num,                  &
      nels=this%mesh%nels,                  &
      lp=this%lp                            &
    )
  END SUBROUTINE MPMCORE_GET_SUPPORT_ELEMENTS

END MODULE MPM_CORE