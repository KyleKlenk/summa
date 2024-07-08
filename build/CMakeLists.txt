cmake_minimum_required(VERSION 3.21 FATAL_ERROR)
project(summa VERSION 4.0.0 LANGUAGES Fortran C CXX)

if(DEFINED ENV{FC})
    set (CMAKE_Fortran_COMPILER $ENV{FC})  #for NextGen need, works for others
else()
    set (CMAKE_Fortran_COMPILER gfortran)
endif()

#=========================================================================================
# SUMMA can be compiled with debug flags by specifiying -DCMAKE_BUILD_TYPE=Debug
#
# There are multiple options for building SUMMA, including using the
# Sundials suite of solvers and the NextGen framework.
# The options are set by specifying -DOPTION=ON when running cmake.
# For example to compile SUMMA with Sundials, use
# cmake -B ../cmake_build -S ../. -DUSE_SUNDIALS=ON
# cmake --build ../cmake_build --target all -j
#
# To compile SUMMA with Actors, see the Summa-Actors repo:
# https://github.com/uofs-simlab/Summa-Actors
#=========================================================================================

# Add options for build type
set(CMAKE_CONFIGURATION_TYPES Release Debug)

# Options: Enable each with cmake -DOPTION=ON
option(USE_SUNDIALS "Use IDA solver from SUNDIALS suite" OFF)
option(USE_NEXTGEN "Use NextGen Framework" OFF)
set(EXT_TARGETS) # list of external targets to link to

# Set Default Executable Name
set(EXEC_NAME summa.exe)
# Set top-level directory
if (USE_NEXTGEN)
    set(F_MASTER "${CMAKE_CURRENT_SOURCE_DIR}/summa")
else()
    set(F_MASTER "${CMAKE_CURRENT_SOURCE_DIR}/..")
endif()

set(PARENT_DIR "${F_MASTER}")
set(EXEC_DIR "${F_MASTER}/bin")

# Handle optional external libraries
if (USE_SUNDIALS)
    message("ENABLING SUNDIALS")
    find_package(SUNDIALS REQUIRED)
    list(APPEND EXT_TARGETS SUNDIALS::fida_mod_shared
                            SUNDIALS::fkinsol_mod_shared)
    add_compile_definitions(SUNDIALS_ACTIVE)
    set(EXEC_NAME summa_sundials.exe)
endif()


if (USE_NEXTGEN)
    message("ENABLING NEXTGEN")
    add_compile_definitions(NGEN_ACTIVE BMI_ACTIVE NGEN_FORCING_ACTIVE
                            NGEN_OUTPUT_ACTIVE)
    add_subdirectory(../iso_c_fortran_bmi ${CMAKE_BINARY_DIR}/iso_c_bmi)
    list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/summa/build/cmake/")
    list(INSERT CMAKE_PREFIX_PATH 0 "${CMAKE_CURRENT_LIST_DIR}/summa/build/cmake/")
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran)
    set(SUMMA_LIB_NAME_CMAKE summabmi)
    set(SUMMA_LIB_DESC_CMAKE "Summa-Sundials BMI Module Shared Library")
else()
    LIST(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake/")
    LIST(INSERT CMAKE_PREFIX_PATH 0 "${CMAKE_SOURCE_DIR}/cmake/")
endif()

get_filename_component(F_MASTER "${F_MASTER}" REALPATH)
get_filename_component(PARENT_DIR "${PARENT_DIR}" REALPATH)
get_filename_component(EXEC_DIR "${EXEC_DIR}" REALPATH)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXEC_DIR})

#=========================================================================================
# SET THE PATHS TO THE LIBRARIES AND INCLUDE FILES
#=========================================================================================

# NetCDF is found with a custom FindNetCDF.cmake file
find_package(NetCDF REQUIRED)
list(APPEND EXT_TARGETS NetCDF::NetCDF)

# OpenBLAS
set(BLA_VENDOR OpenBLAS) # MKL
find_package(OpenBLAS REQUIRED)
list(APPEND EXT_TARGETS OpenBLAS::OpenBLAS)

# Set compiler flags
set(FLAGS_OPT $ENV{FLAGS_OPT}) # get optional user-specified flags from environment variables
if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("\nSetting SUMMA Debug Options")
    add_compile_definitions(DEBUG)
    set(FLAGS_NOAH -g -O0 -fbacktrace -fbounds-check -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors ${FLAGS_OPT})
    set(FLAGS_ALL -g -O0 -fbacktrace -fbounds-check -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors -cpp ${FLAGS_OPT})
    set(FLAGS_CXX -g -O0 -fbounds-check -Wfatal-errors -std=c++17 ${FLAGS_OPT})
else()
    message("\nSetting SUMMA Release Options")
    set(FLAGS_NOAH -O3 -ffree-form -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors ${FLAGS_OPT})
    set(FLAGS_ALL  -O3 -ffree-line-length-none -fmax-errors=0 -fPIC -Wfatal-errors -cpp ${FLAGS_OPT})
    set(FLAGS_CXX  -O3 -Wfatal-errors -std=c++17 ${FLAGS_OPT})
endif()

#=========================================================================================
# COMPILE PART 1: Define directory paths
#=========================================================================================

# Define directories that contains source code
set(DRIVER_DIR ${F_MASTER}/build/source/driver)
set(DSHARE_DIR ${F_MASTER}/build/source/dshare)
set(ENGINE_DIR ${F_MASTER}/build/source/engine)
set(HOOKUP_DIR ${F_MASTER}/build/source/hookup)
set(NETCDF_DIR ${F_MASTER}/build/source/netcdf)
set(NOAHMP_DIR ${F_MASTER}/build/source/noah-mp)

#=========================================================================================
# COMPILE PART 2: Assemble all of the SUMMA sub-routines
#=========================================================================================

# SUMMA Source Files are defined in the CMakeLists.txt file in the subdirectory
add_subdirectory(${F_MASTER}/build/source/)


#=========================================================================================
# COMPILE PART 3: Collect the subroutines into build groups depending on build type
#=========================================================================================

set(COMM_ALL ${NRPROC} ${HOOKUP} ${DATAMS} ${UTILMS})
set(SUMMA_ALL ${NETCDF} ${PRELIM} ${MODRUN} ${SOLVER} ${DRIVER})
set(MAIN_SUMMA ${DRIVER_DIR}/summa_driver.f90)

if (USE_NEXTGEN)
    set(SUMMA_ALL ${SUMMA_ALL} ${DRIVER_NEXTGEN})
endif()

if (USE_SUNDIALS)
    set(COMM_ALL ${COMM_ALL} ${DATAMS_SUNDIALS} ${UTILMS_SUNDIALS})
    set(SUMMA_ALL ${SUMMA_ALL} ${MODRUN_SUNDIALS} ${SOLVER_SUNDIALS})
endif()


# Define version number, not working correctly
set(VERSIONFILE     ${DRIVER_DIR}/summaversion.inc)
execute_process(COMMAND "    ${GIT_EXECUTABLE} tag | tail -n 1" OUTPUT_VARIABLE VERSION)
execute_process(COMMAND "date" OUTPUT_VARIABLE BULTTIM)
execute_process(COMMAND "    ${GIT_EXECUTABLE} describe --long --all --always | sed -e's/heads\///'" OUTPUT_VARIABLE GITBRCH)
execute_process(COMMAND "    ${GIT_EXECUTABLE} rev-parse HEAD" OUTPUT_VARIABLE GITHASH)


#=========================================================================================
# COMPILE PART 4: Do the compilation
#=========================================================================================

# update version information, not working correctly
file(WRITE  ${VERSIONFILE} "character(len=64), parameter     :: summaVersion = '${VERSION}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: buildTime = ''\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: gitBranch = '${GITBRCH}'\n")
file(APPEND ${VERSIONFILE} "character(len=64), parameter     :: gitHash = '${GITHASH}'")

# Build SUMMA_NOAHMP Object
add_library(SUMMA_NOAHMP OBJECT ${NOAHMP} ${NRUTIL})
target_compile_options(SUMMA_NOAHMP PRIVATE ${FLAGS_NOAH})

# Build SUMMA_COMM Object
add_library(SUMMA_COMM OBJECT ${COMM_ALL})
target_compile_options(SUMMA_COMM PRIVATE ${FLAGS_ALL})
target_link_libraries(SUMMA_COMM PUBLIC SUMMA_NOAHMP ${EXT_TARGETS}) # added flags to the link step

# For NextGen, build SUMMA Shared Library and add the outside BMI libraries
if(USE_NEXTGEN)
    if(WIN32)
        add_library(summabmi ${SUMMA_ALL})
    else()
        add_library(summabmi SHARED ${SUMMA_ALL})
    endif()
    target_compile_options(summabmi PRIVATE ${FLAGS_ALL})
	target_link_libraries(summabmi PUBLIC ${EXT_TARGETS} iso_c_bmi SUMMA_NOAHMP SUMMA_COMM)

    set_target_properties(summabmi PROPERTIES VERSION ${PROJECT_VERSION})
    include(GNUInstallDirs)

    install(TARGETS summabmi
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    configure_file(summabmi.pc.in summabmi.pc @ONLY)
    install(FILES ${CMAKE_BINARY_DIR}/summabmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

else()
	add_library(summa SHARED ${SUMMA_ALL})
	target_compile_options(summa PRIVATE ${FLAGS_ALL})
	target_link_libraries(summa PUBLIC ${EXT_TARGETS} SUMMA_NOAHMP SUMMA_COMM)
    add_executable(${EXEC_NAME} ${MAIN_SUMMA})
    set_property(TARGET ${EXEC_NAME} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${EXEC_NAME} summa ${EXT_TARGETS}) # added flags to the link step
endif()