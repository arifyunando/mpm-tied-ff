        !COMPILER-GENERATED INTERFACE MODULE: Sun May 19 04:00:12 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VMFLOW__genmod
          INTERFACE 
            SUBROUTINE VMFLOW(STRESS,DSBAR,VMFL)
              REAL(KIND=8), INTENT(IN) :: STRESS(:)
              REAL(KIND=8), INTENT(IN) :: DSBAR
              REAL(KIND=8), INTENT(OUT) :: VMFL(:)
            END SUBROUTINE VMFLOW
          END INTERFACE 
        END MODULE VMFLOW__genmod
