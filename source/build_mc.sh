#!/bin/bash
#
# build_pecube.sh - compiles pecube executable
#

# Standard include statements for building detrital monte carlo code
KP_INCLUDE=.

#
# Portland Group Fortran compiler
#pgf90 -O3 -c nrtype.f90 nrutil.f90
#pgf90 -O3 -c *.f90
#pgf90 -O3 -o dmc *.o

# Portland Group Fortran compiler w/ debug/warnings
#pgf90 -C -w -c nrtype.f90 nrutil.f90
#pgf90 -C -w -c *.f90
#pgf90 -C -w -o dmc *.o

# Intel Fortran compiler
#ifort -O3 -c nrtype.f90 nrutil.f90
#ifort -O3 -c *.f90
#ifort -O3 -o dmc *.o

# Intel Fortran compiler
#ifort -c nrtype.f90 nrutil.f90
#ifort -c *.f90
#ifort -o dmc *.o

# Intel Fortran compiler optimized for Core 2 Duos
#ifort -fast -c *.f90 *.f
#ifort -fast -o pecube *.o

# Intel Fortran compiler w/ debug/warning flags
#ifort -c -C -w -debug all nrtype.f90 nrutil.f90
#ifort -c -C -w -debug all *.f90
#ifort -o dmc -C -w -debug all *.o

# Latest GNU Fortran compiler (slower than ifort or pgf90, but free)
#gfortran -O3 -c nrtype.f90 nrutil.f90
#gfortran -O3 -c *.f90
#gfortran -O3 -o dmc *.o

# Latest GNU Fortran compiler (slower than ifort or pgf90, but free)
gfortran -O2 -c nrtype.f90 nrutil.f90
gfortran -O2 -c *.f90
gfortran -O2 -o dmc *.o

# Latest GNU Fortran compiler w/ debug/warnings
#gfortran -C -w -c nrtype.f90 nrutil.f90
#gfortran -C -w -c *.f90
#gfortran -C -w -o dmc *.o

# Another open source Fortran compiler
#g95 -c *.f*
#g95 -o pecube *.o

# Remove object files after compilation
rm -f *.o
