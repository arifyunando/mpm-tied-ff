PROGRAM Hydrostatic_Column
  USE MPMCore
  IMPLICIT NONE

  !*** PREPARATION PHASE ***!
  !--- Initiate Program Settings
  !--- Initiate Domain
  !--- Set Material Properties
  !--- Build Computational Grid
  !--- Build Geometry
  !--- Infuse Material Points
  !--- Setup Initial Conditions and Boundary Conditions
  !--- Calculate Initial Properties of Material Points
  !--- Initiate Time Stepping Variables
  
  !*** TIME STEPPING ***!
  !--- Body Solution
  !--- Stiffness
  !--- Combined Stiffness
  !--- Plastic Strain Iteration
  !--- Reverse Mapping (Node -> MP)
  !--- Save Material Point Data (VTK)
  !--- Reallocate Memory
  !--- Determine and Activate New Element for Next Time Step
    
END PROGRAM Hydrostatic_Column