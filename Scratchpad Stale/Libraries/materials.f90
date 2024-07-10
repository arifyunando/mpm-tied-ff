MODULE CLASS_MATERIAL
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE::iwp=SELECTED_REAL_KIND(15)

  TYPE::material_model
    INTEGER :: id
    CHARACTER(:), ALLOCATABLE :: type

    ! Material Properties
    REAL(iwp) :: density, cohesion, phi, dilation
    
    ! Elasticity Parameters
    REAL(iwp) :: modulus, poisson_ratio

    CONTAINS

    PROCEDURE :: LOAD => MATERIAL_LOAD
  END TYPE

  CONTAINS

  SUBROUTINE MATERIAL_LOAD(this, index, input_json)
    USE JSON_MODULE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: index
    TYPE(json_file), INTENT(INOUT) :: input_json
    CLASS(material_model), INTENT(OUT) :: this
    CHARACTER(1) :: i_char
    
    WRITE(i_char, '(I1)') index
    CALL input_json%GET('materials('//i_char//').id', this%id)
    CALL input_json%GET('materials('//i_char//').type', this%type)
    
    !--- Load Common Material Parameter
    ! Material Properties
    CALL input_json%GET('materials('//i_char//').density', this%density)
    CALL input_json%GET('materials('//i_char//').cohesion', this%cohesion)
    CALL input_json%GET('materials('//i_char//').friction', this%phi)
    CALL input_json%GET('materials('//i_char//').dilation', this%dilation)

    ! Elastic Parameters
    CALL input_json%GET('materials('//i_char//').poisson_ratio', this%poisson_ratio)
    CALL input_json%GET('materials('//i_char//').youngs_modulus', this%modulus)
  END SUBROUTINE

END MODULE CLASS_MATERIAL