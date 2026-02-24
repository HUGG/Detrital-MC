### Start of the makefile ---------------------------------------------------###
NAME=det_mc

### Use Fortran stdlib
# This is for building on macOS with the Fortran stdlib installed using homebrew
# Adjust accordingly!
stdlib_install_dir := /opt/homebrew
libdir := $(stdlib_install_dir)/lib
moduledir := $(stdlib_install_dir)/include/fortran_stdlib/GNU-15.2.0

### Define compilers --------------------------------------------------------###
# GNU
F90    = gfortran

# Intel
#F90    = /usr/bin/ifort

# Define compiler flags -----------------------------------------------------###
#FFLAGS  = -O2
FFLAGS  = -O2 -ffree-line-length-160
#FFLAGS  = -O0 -g -w -traceback
#FFLAGS  = -O0 -g -w -fbacktrace -ffree-line-length-160 -Wall -Wextra
#FFLAGS  = -O0 -g -Wall -Wextra -pedantic -fimplicit-none -fbounds-check -fbacktrace
#FFLAGS  = -C -w
#FFLAGS  = -C -w -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic
#FFLAGS  = -O2
#FFLAGS  = -O3

#AR=ar
#BITS=_v3.3
BITS =

# Fortran stdlib stuff
# NOTE: You may need to adjust the LDLIBS depending on how the Fortran stdlib
#       was installed
FFLAGS  = -I$(moduledir)
LDFLAGS = -L$(libdir)
LDLIBS  = -lfortran_stdlib_sorting \
-lfortran_stdlib_bitsets \
-lfortran_stdlib_strings \
-lfortran_stdlib_core \
-lfortran_stdlib \
-lm

# Compile rules
#
.SUFFIXES:
.SUFFIXES: .o .c .f .f90
.f90.o:
	$(F90) $(FFLAGS) -c $<

BIN=$(NAME)$(BITS)

OBJECTS = \
nrtype.o \
nrutil.o \
module_definitions.o \
module_pdf_functions.o \
indexx.o \
init_random_seed.o \
kptwo.o \
read_input_file.o \
sort2.o \
detrital_mc.o

all: $(BIN)
$(BIN): $(OBJECTS)
	$(F90) $(FFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $(BIN)

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
