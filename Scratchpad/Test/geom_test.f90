PROGRAM TEST_GEOMETRY_GENERATOR
  IMPLICIT NONE
  CALL TEST_GEOMETRY_GENERATOR()

  CONTAINS

  SUBROUTINE TEST_OFFSET_GENERATION
    USE GEOMETRY_GENERATOR
    USE IO
    IMPLICIT NONE
    INTEGER, PARAMETER :: iwp=SELECTED_REAL_KIND(15)
    REAL(iwp),ALLOCATABLE::coord(:,:)
    INTEGER,ALLOCATABLE::num(:,:)
    INTEGER::nn,nels
  
    
    CALL RECTANGULAR_2D(coord,num,nn,nels,20,15,60.0_iwp,45.0_iwp,4)
    PRINT*, nn, nels
    DO nn=1,nels
      print*,num(:,nn)
    END DO
    CALL IO_PARAVIEW(0,4,coord,num)
    CALL RECTANGULAR_2D(coord,num,nn,nels,20,15,60.0_iwp,45.0_iwp,4,offsetx=25)
    CALL IO_PARAVIEW(0,4,coord,num,argv="offset")
  END SUBROUTINE TEST_OFFSET_GENERATION
END PROGRAM