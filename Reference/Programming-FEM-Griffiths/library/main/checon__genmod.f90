        !COMPILER-GENERATED INTERFACE MODULE: Sun May 19 04:01:57 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECON__genmod
          INTERFACE 
            SUBROUTINE CHECON(LOADS,OLDLDS,TOL,CONVERGED)
              REAL(KIND=8), INTENT(IN) :: LOADS(0:)
              REAL(KIND=8), INTENT(INOUT) :: OLDLDS(0:)
              REAL(KIND=8), INTENT(IN) :: TOL
              LOGICAL(KIND=4), INTENT(OUT) :: CONVERGED
            END SUBROUTINE CHECON
          END INTERFACE 
        END MODULE CHECON__genmod
