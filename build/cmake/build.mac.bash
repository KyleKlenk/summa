#!/bin/bash

# build SUMMA on a Mac using Bash, from cmake directory run this as ./build.mac.bash
# Environment variables that must be set for cmake: FC, LINK_DIRS, INCLUDES_DIRS, and LIBRARY_LINKS

# Mac Example using MacPorts:
export FC=gfortran                                            # Fortran compiler family
export LINK_DIRS=/opt/local/lib                               # Link directories for cmake
export INCLUDES_DIRS='/opt/local/include;/opt/local/lib'      # directories for INCLUDES cmake variable (cmake uses semicolons as separators)
export LIBRARY_LINKS='-llapack;-lgfortran;-lnetcdff;-lnetcdf' # list of library links (cmake uses semicolons as separators)

cmake -B ../cmake_build -S . -DCMAKE_BUILD_TYPE=Sundials
cmake --build ../cmake_build --target all