        !COMPILER-GENERATED INTERFACE MODULE: Sun May 19 04:01:50 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FMRMAT__genmod
          INTERFACE 
            SUBROUTINE FMRMAT(VMFL,DSBAR,DLAM,DEE,RMAT)
              REAL(KIND=8), INTENT(IN) :: VMFL(:)
              REAL(KIND=8), INTENT(IN) :: DSBAR
              REAL(KIND=8), INTENT(IN) :: DLAM
              REAL(KIND=8), INTENT(IN) :: DEE(:,:)
              REAL(KIND=8), INTENT(OUT) :: RMAT(:,:)
            END SUBROUTINE FMRMAT
          END INTERFACE 
        END MODULE FMRMAT__genmod
