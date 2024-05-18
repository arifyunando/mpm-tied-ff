MODULE CLASS_MESH
  USE JSON_MODULE
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)

  TYPE::mpm_grid
    INTEGER :: id
    ! Constants
    INTEGER :: nnodes, nels

    ! Tracking variables
    REAL(iwp), ALLOCATABLE :: g_coord(:,:)  ! Node Coordinates
    INTEGER, ALLOCATABLE :: num(:,:)        ! Node Index that build an element
    INTEGER, ALLOCATABLE :: nf(:,:)         ! Track degree of freedoms

    CONTAINS

    PROCEDURE :: LOAD_MESH => CLASS_MESH_LOAD_MESH
    PROCEDURE :: FORM_GLOBAL_NF => CLASS_MESH_FORM_GLOBAL_NF
  END TYPE

  CONTAINS

  SUBROUTINE CLASS_MESH_LOAD_MESH(this, directory, input_json, nodof, nst)
    USE JSON_MODULE
    USE IO
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN)       :: directory
    TYPE(json_file), INTENT(INOUT) :: input_json
    CLASS(mpm_grid), INTENT(INOUT) :: this
    INTEGER, OPTIONAL, INTENT(INOUT) :: nodof, nst
    CHARACTER(:), ALLOCATABLE :: filename, entity_filename
    TYPE(json_file) :: entity_json
    INTEGER :: def_nodof=2, def_nst=4

    ! Get Mesh Filename
    CALL input_json%GET('mesh.mesh', filename)

    ! Get mesh nodal locations
    CALL IO_LOAD_MESH(trim(directory), trim(filename), this%g_coord, this%num)
    if (present(nodof)) def_nodof = 2
    if (present(nst)) def_nst = 4

    ! Determine Mesh Variable
    this%nnodes = ubound(this%g_coord, 2)
    this%nels = ubound(this%num, 2)

    ALLOCATE(this%nf(def_nodof, this%nnodes))

  END SUBROUTINE CLASS_MESH_LOAD_MESH


  SUBROUTINE CLASS_MESH_FORM_GLOBAL_NF(this, directory, input_json)
    USE JSON_MODULE
    USE IO
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN)       :: directory
    TYPE(json_file), INTENT(INOUT) :: input_json
    CLASS(mpm_grid), INTENT(INOUT) :: this
    TYPE(json_file) :: entity_json
    LOGICAL :: found=.true.
    INTEGER :: i, nset_id, dir
    INTEGER, ALLOCATABLE :: nodes(:)
    CHARACTER(:), ALLOCATABLE :: entity_filename
    CHARACTER(1) :: ic

    CALL input_json%GET('mesh.entity_sets', entity_filename)
    CALL IO_LOAD_JSON(trim(directory), trim(entity_filename), entity_json)

    this%nf = 1
    i = 1
    DO WHILE(found)
      write(ic, '(I1)') i
      CALL input_json%GET(                                                      &
      'mesh.boundary_conditions.displacement_constraints('//ic//').nset_id',    &
      nset_id, found                                                            &
      )
      CALL input_json%GET(                                                      &
      'mesh.boundary_conditions.displacement_constraints('//ic//').dir',        &
      dir, found                                                                &
      )
      IF (found) THEN
        write(ic, '(I1)') nset_id
        CALL entity_json%GET('node_sets('//ic//').set', nodes)
        this%nf(dir,nodes) = 0
      END IF
      i = i+1
    END DO
  END SUBROUTINE CLASS_MESH_FORM_GLOBAL_NF

END MODULE CLASS_MESH