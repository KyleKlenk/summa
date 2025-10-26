#!/bin/bash

# build SUMMA on a Mac using Bash, from cmake directory run this as ./build.mac.bash
# Environment variables may be set within this script (see examples below) or in the terminal environment before executing this script
# Actual settings may vary

# Mac Example using MacPorts:
export FC=/opt/homebrew/bin/gfortran                             # Fortran compiler family
#export FLAGS_OPT="-flto=1"                                   # -flto=1 is slow to compile, but might want to use
export LIBRARY_LINKS='-llapack'                               # list of library links
export SUNDIALS_DIR=$HOME/local/sundials

cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DSPECIFY_LAPACK_LINKS=ON -DCMAKE_BUILD_TYPE=Release
#cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON -DSPECIFY_LAPACK_LINKS=ON -DCMAKE_BUILD_TYPE=Debug
#cmake -B ../cmake_build -S ../. -DCMAKE_BUILD_TYPE=Debug \
#  -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -Wall -Wextra -Wuninitialized -Wmaybe-uninitialized -fcheck=all -finit-real=snan -finit-integer=2147483647"
cmake --build ../cmake_build --target all -j
