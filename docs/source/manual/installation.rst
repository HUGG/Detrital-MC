Installation
============

Dependencies
------------

To compile the ``det_mc`` executable, you need the following software on your system:

`input file <https://github.com/HUGG/Detrital-MC/blob/master/input/det_mc_input.txt>`__
`latest source code release <https://github.com/HUGG/Detrital-MC/releases/>`__


- A Fortran compiler (`gfortran <https://gcc.gnu.org/fortran/>`__ is recommended)
- The `Fortran Standard Library <https://fortran-lang.github.io/stdlib/>`__
- `GNU Make <https://www.gnu.org/software/make/>`__ or `CMake <https://cmake.org/>`__

### Downloading the Detrital MC source code

To build Detrital MC, you should first download the `latest source code release <https://github.com/HUGG/Detrital-MC/releases/>`__ as a ``.zip`` or ``.tar.gz`` file and extract its contents. For those who are familiar, you can also clone the Detrital MC git repository and build using that.

### Building using GNU Make

1. Edit the ``source/Makefile`` file to ensure the installation location of the Fortran Standard Library is correct for your system.
2. In a terminal, navigate to the ``source`` subdirectory of Detrital MC.
3. Run 

   .. code-block:: bash
   
      make
      make install
    
   to compile Detrital MC and install it in the ``bin`` subdirectory.

### Building using CMake

1. Edit the ``source/CMakeLists.txt`` file to ensure the installation location of the Fortran Standard Library is correct for your system.
2. In a terminal, navigate to the base directory of Detrital MC.
3. Run 

   .. code-block:: bash

      mkdir build
      cmake -B build -S source
      cmake --build build
      cmake --install build --prefix .
    
   to compile Detrital MC and install it in the ``bin`` subdirectory.