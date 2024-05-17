#! /bin/bash

ifx -g2 -c ../library/main/*.f90
ifx -g2 -c ../library/geom/*.f90
ar -r femlib.a ./*.o
rm *.o 