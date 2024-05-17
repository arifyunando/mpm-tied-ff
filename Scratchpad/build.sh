#! /bin/bash

# initiate folder
mkdir -p build
cd build

# build static library
ifx -g2 -O0 -c ../Libraries/json_module.f90
ifx -g2 -O0 -c ../Libraries/logging.f90

ifx -g2 -O0 -c ../Libraries/mpm_io.f90                 

ifx -g2 -O0 -c ../Libraries/mpm_core.f90                 
ifx -g2 -O0 -c ../Libraries/fem_core.f90                 


# build application
ifx -g2 -O0 -static-intel ../Tied-FF-MPM.f90 ./*.o -o Tied-FF-MPM.run
mv *.mod ../modules/
rm *.o