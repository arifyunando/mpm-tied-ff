# MPM_Site_response
Development of a MPM code to perform site response analysis.

The main goal of this repository is to enhance a basic implicit MPM code to simulate earthquakes. 

NOTES:

The codes included in the branches are those used to analyse the examples of the thesis 
"Investigation of MPM inaccuracies, contact simulation and robust implementation for geotechnical problems".

The codes follow similar structures (see theses pages 21, 22, and 78) and most of the subroutines are detailed 
in "Programming the Finite Element Method" by Smith, D.V. Griffiths, L. Margetts.

The codes have been successfully compiled and tested using visual studio Community 2019 (recommended) 
(https://my.visualstudio.com/Downloads?q=visual%20studio%202019&wt.mc_id=o~msft~vscom~older-downloads) 
plus intel oneAPI base toolkit and intel oneAPI HPC Toolkit (https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html),
and it has been observed that other compiling procedures and interfaces can also be used (e.g., CodeBlock). Nevertheless, some codes use
MKL libraries (pardiso solver), making the use of Visual Studio Community 2019, oneAPI base toolkit and intel oneAPI HPC Toolkit mandatory. Also, it
is necessary to activate MKL libraries in Visual Studio preferences to allow the compilation process.


The codes have been cleaned to avoid confusion. Nevertheless, variables, vectors or matrices can still be present in the codes and not be used.

When running the codes, a couple of folders are necessary (i.e. Paraview2 or Paraview) to save the data created. 
Check the subroutines "paraview2" and "point_viz2" used to create the results files and to know which folders are needed to save the data.

Jose Leon Gonzalez Acosta lgonzalez.a87@gmail.com 20-06-2023
