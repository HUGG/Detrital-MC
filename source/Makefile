### Start of the makefile ---------------------------------------------------###
NAME=det_mc

### Define compilers --------------------------------------------------------###
# GNU
#F90    = /usr/local/bin/gfortran
#F77    = /usr/local/bin/gfortran
#CC     = /usr/local/bin/gcc

# Intel
F90    = /usr/bin/ifort
F77    = /usr/bin/ifort
CC     = /usr/bin/icc

# Define compiler flags -----------------------------------------------------###
#F90FLAGS  = -O0
#F77FLAGS  = -O0
#CFLAGS  = -O0
#F90FLAGS  = -O0 -g -w -traceback
#F77FLAGS  = -O0 -g -w -traceback
#CFLAGS  = -O0 -g -w -traceback
#F90FLAGS  = -O0 -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace
#F77FLAGS  = -O0 -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace
#CFLAGS  = -O0 -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace
#F90FLAGS  = -C -w
#F77FLAGS  = -C -w
#CFLAGS  = -C -w
#F90FLAGS  = -C -w -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
#F77FLAGS  = -C -w -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
#CFLAGS  = -C -w -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
F90FLAGS  = -O2
F77FLAGS  = -O2
CFLAGS  = -O2
#F90FLAGS  = -O3
#F77FLAGS  = -O3
#CFLAGS  = -O3
#INT8FLAG = -fdefault-integer-8

#AR=ar
BITS=_v1.1
#BITS=_v1.1_int8

# Compile rules
#
.SUFFIXES:
.SUFFIXES: .o .c .f .f90
.f90.o:
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<
.f.o:
	$(F77) $(F77FLAGS) $(INCLUDES) -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<

BIN=$(NAME)$(BITS)

OBJECTS = \
nrtype.o \
nrutil.o \
indexx.o \
init_random_seed.o \
kptwo.o \
sort2.o \
pdf_functions.o \
detrital_mc.o
#detrital_mc_int8.o

all: $(BIN)
$(BIN): $(OBJECTS)
	$(F90) $(FFLAGS) $(OPTIONS) $(OBJECTS) $(LIBS) -o $(BIN)
	#$(F90) $(FFLAGS) $(OPTIONS) $(OBJECTS) $(INCLUDE) $(LIBS) -o $(BIN)

.PHONY: install
install:
	mv -f $(BIN) ../bin

.PHONY: objclean
objclean:
	rm -f *.o *.mod

.PHONY: distclean
distclean: objclean
	rm -f $(BIN)

# End of the makefile
