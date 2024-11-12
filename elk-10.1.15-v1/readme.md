This is a modified verison of Elk-10.1.15. 

Modifications:
 1. A global magnetic field now can be applied in ULR calculations
 2. ULR magnetic field, include random magnetic field, is now initialised in R-mesh rahter than Q-mesh.
 3. Add new subroutines to read and write external vector potential in R-mesh


New Subroutines:
 1. readpotvecr.f90
 2. writepotvecr.f90
 3. bfuinit.f90

New global variables:
 1. trdbfcr

Changed src files
 1. readpotvecr.f90
 2. writepotvecr.f90
 3. bfuinit.f90
 4. init0.f90
 5. gndstulr.f90
 6. potuinit.f90
 7. readinput.f90
 8. modulr.f90
 9. ./src/Makefile
