MODULE MPM_CORE
  USE CLASS_PARTICLE
  USE CLASS_MATERIAL
  USE CLASS_MESH
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)
  REAL(iwp), PARAMETER, PRIVATE :: zero=0.0_iwp, one=1.0_iwp, two=2.0_iwp
  
  TYPE, EXTENDS(mpm_particles)::mpm_body
    !
    ! MPM Particle Base Class
    !
    CLASS(material_model), ALLOCATABLE :: material
    CLASS(mpm_grid), ALLOCATABLE :: mesh

    CONTAINS

    PROCEDURE :: SET_MESH => MPMCORE_SET_MESH
    PROCEDURE :: LOAD_MATERIAL => MPMCORE_LOAD_MATERIAL
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


  SUBROUTINE MPMCORE_LOAD_MATERIAL(this, input_json)
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
  END SUBROUTINE MPMCORE_LOAD_MATERIAL

  ! TODO
  SUBROUTINE MPMCORE_ACTIVATE_NODE()
  !   !
  !   !
  !   !
  !   IMPLICIT NONE
  !   Flags: DO bod=1,size(mbod)
  !   IF(mbod(bod)%nmps>1)THEN  
  !       mbod(bod)%ale=1
  !       mbod(bod)%d_ele=0
  !       mbod(bod)%c_ele=0
  !       DO i=1, mbod(bod)%nmps
  !           ! search in every elements
  !           inner: DO iel=1,nels
  !             ! calculate element column number based on current MP location
  !             colx=mbod(bod)%gm_coord(1,i)/cellsize+1
  !             ! calculate element row number based on current MP location
  !             rowy=ABS(mbod(bod)%gm_coord(2,i))/cellsize+1
  !             ! calculate element number of the current MP
  !             ielloc=(rowy-1.0)*nx1+colx
  !             ! assign element number to a_le
  !             mbod(bod)%a_ele(i)=ielloc
  !             ! get node coordinates of the element
  !             num=g_num(:,ielloc)
  !             coord=TRANSPOSE(g_coord(:,num))
              
  !             ! get current global MP coordinate
  !             sp_coord(:,1)=mbod(bod)%gm_coord(:,i)
  !             ! get previous local MP coordinate
  !             lp_coord(1,:)=mbod(bod)%mpoints(i,:)
  !             ! calculate local coordinate of the mp
  !             CALL floc(coord,sp_coord,lp_coord,i)
  !             ! reassign local coordinate back to mpoints
  !             mbod(bod)%mpoints(i,:)=lp_coord(1,:)

  !             ! mark activated node in d_ele
  !             mbod(bod)%d_ele(ielloc)=1

  !             IF(i>1)THEN
  !               DO j=1,i-1
  !                 IF(mbod(bod)%a_ele(j)==mbod(bod)%a_ele(i)) THEN
  !                   EXIT inner
  !                 END IF
  !               END DO
  !               ! tally the numbers of activated elements
  !               mbod(bod)%ale=mbod(bod)%ale+1
  !             END IF
              
  !             EXIT inner
  !           END DO inner

  !       END DO
  !       ! Count MP in each cells
  !       CALL couma(nels,mbod(bod)%nmps,mbod(bod)%a_ele,mbod(bod)%c_ele,mbod(bod)%k_ele,etype)

  !       ! ------------- list the corresponding activated elements --------------

  !       k=1
  !       mbod(bod)%b=0
  !       DO i=1,mbod(bod)%nmps
  !         ! Get element number of current MP
  !         iel=mbod(bod)%a_ele(i)
  !         ! Loop through the MP of inhabited element
  !         DO j=1,mbod(bod)%c_ele(iel)
  !             ! 
  !             IF(mbod(bod)%b(mbod(bod)%k_ele(iel-1)+j)==0) THEN
                
  !               mbod(bod)%b(mbod(bod)%k_ele(iel-1)+j)=i
  !               EXIT
  !             END IF
  !         END DO
  !       END DO
  !     END IF
  !   END DO Flags
  END SUBROUTINE MPMCORE_ACTIVATE_NODE
END MODULE MPM_CORE