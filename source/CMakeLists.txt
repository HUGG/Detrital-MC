cmake_minimum_required(VERSION 3.15)

# Define project name
project(det_mc VERSION 3.3
DESCRIPTION "Program for plotting and comparing detrital age observations and predictions as PDFs, CDFs or ECDFs. Monte Carlo sampling of predicted ages can be performed as well."
LANGUAGES Fortran)
enable_language(Fortran)

# Find installed fortran_stdlib
# NOTE: This assumes installation using homebrew for macOS. Adjust accordingly!
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/lib")
find_package(fortran_stdlib REQUIRED)

# add the executable
add_executable(det_mc indexx.f90 kptwo.f90 module_definitions.f90 nrtype.f90 read_input_file.f90 detrital_mc.f90 init_random_seed.f90 module_pdf_functions.f90 nrutil.f90 sort2.f90)

# Link against stdlib
target_link_libraries(det_mc PRIVATE fortran_stdlib::fortran_stdlib_sorting)

# Define programs to install
install(TARGETS det_mc DESTINATION bin)