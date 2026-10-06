# CMake build of PP-RMTA against a CMake-built and *installed* Quantum ESPRESSO
#
#  rm -rf build
#  FC=mpiifx CC=mpiicx cmake -S . -B build -DQE_ROOT=<QE install prefix>
#  cmake --build build
#
# Counterpart of `include $(QE_ROOT)/make.inc`: the build configuration of QE
# (preprocessor definitions, HDF5/libxc/ScaLAPACK/... usage) is read from QE's
# own installed <prefix>/lib*/cmake/qe/qeTargets.cmake, so nothing here needs
# editing when QE is rebuilt with different options. If QE's dependencies live
# outside the default search paths and QE did not record them, pass them via
# -DCMAKE_PREFIX_PATH=...
#
# QE's exported package cannot simply be find_package()'d (relative include
# paths, undefined MPI::/Mbd::/Wannier90:: targets), hence the parsing.
#
cmake_minimum_required(VERSION 3.20)
project(pp-rmta LANGUAGES Fortran C) # C: needed by FindBLAS/FindHDF5 probes

set(QE_ROOT "$ENV{QE_ROOT}" CACHE PATH
    "Prefix of the CMake-installed Quantum ESPRESSO")
file(GLOB QE_TARGETS_FILE ${QE_ROOT}/lib*/cmake/qe/qeTargets.cmake)
if(NOT QE_TARGETS_FILE)
  message(FATAL_ERROR "No QE install found in QE_ROOT='${QE_ROOT}' "
          "(expected lib*/cmake/qe/qeTargets.cmake; "
          "run 'cmake --install' for QE)")
endif()

# read QE's configuration ======================================================
file(READ ${QE_TARGETS_FILE} _qe)
# brackets would break list handling
string(REGEX REPLACE "[][]" "" _qe "${_qe}")
# keep each quoted value in one piece
string(REPLACE ";" "|" _qe "${_qe}")

# collect all values of the given INTERFACE_* property into a list
function(qe_property prop out)
  string(REGEX MATCHALL "${prop} \"[^\"]*\"" _m "${_qe}")
  set(_vals "")
  foreach(_x IN LISTS _m)
    string(REGEX REPLACE "^${prop} \"(.*)\"$" "\\1" _x "${_x}")
    string(REPLACE "|" ";" _x "${_x}")
    list(APPEND _vals ${_x})
  endforeach()
  list(FILTER _vals EXCLUDE REGEX "^\\\\?\\$<")   # drop generator expressions
  list(REMOVE_DUPLICATES _vals)
  set(${out} "${_vals}" PARENT_SCOPE)
endfunction()

qe_property(INTERFACE_COMPILE_DEFINITIONS QE_DEFS) # __MPI; __HDF5;...
qe_property(INTERFACE_INCLUDE_DIRECTORIES QE_INCS)
qe_property(INTERFACE_LINK_LIBRARIES      QE_LIBS)
# keep external (absolute) include dirs and libraries only
list(FILTER QE_INCS INCLUDE REGEX "^/")
list(FILTER QE_LIBS INCLUDE REGEX "(^/.*\\.(so|a|dylib)|^-l)")
list(FILTER QE_LIBS EXCLUDE REGEX "^${QE_ROOT}")
message(STATUS "QE definitions: ${QE_DEFS}")

function(qe_uses def out)
  if("${def}" IN_LIST QE_DEFS)
    set(${out} TRUE PARENT_SCOPE)
  else()
    set(${out} FALSE PARENT_SCOPE)
  endif()
endfunction()

# QE sources: modules and archives =============================================
# one module directory per QE component
file(GLOB QE_MODDIRS LIST_DIRECTORIES true ${QE_ROOT}/include/qe/*)
# archives depend on each other circularly -> link group
file(GLOB QE_ARCHIVES ${QE_ROOT}/lib*/libqe_*.a ${QE_ROOT}/lib*/libmbd.a)

# third-party dependencies, as used by QE ======================================
# search where QE found its dependencies: prefixes of its external include dirs
foreach(_inc IN LISTS QE_INCS)
  get_filename_component(_prefix ${_inc} DIRECTORY)
  list(APPEND CMAKE_PREFIX_PATH ${_prefix})
endforeach()
set(DEP_LIBS "")
qe_uses(__MPI mpi)
if(mpi)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  list(APPEND DEP_LIBS MPI::MPI_Fortran)
endif()
qe_uses(__HDF5 hdf5)
if(hdf5)
  find_package(HDF5 REQUIRED COMPONENTS Fortran)
  list(APPEND DEP_LIBS HDF5::HDF5)
endif()
qe_uses(__LIBXC libxc)
if(libxc)
  find_library(XCF03 xcf03 REQUIRED)
  find_library(XC xc REQUIRED)
  list(APPEND DEP_LIBS ${XCF03} ${XC})
endif()
if(_qe MATCHES "OpenMP::OpenMP_Fortran")
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
  list(APPEND DEP_LIBS OpenMP::OpenMP_Fortran)
endif()
# BLAS/LAPACK(/ScaLAPACK): QE records absolute library paths when it found
# them itself (e.g. MKL ScaLAPACK); otherwise search like QE does (MKL for Intel
# compilers; it also provides DFTI).
if(NOT QE_LIBS MATCHES "mkl|lapack|openblas")
  if(NOT BLA_VENDOR)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
      set(BLA_VENDOR Intel10_64lp_seq)
    else()
      set(BLA_VENDOR All)
    endif()
  endif()
  find_package(LAPACK REQUIRED)
  list(APPEND DEP_LIBS ${LAPACK_LIBRARIES})
endif()

# PP-RMTA ======================================================================
add_compile_definitions(${QE_DEFS})
set(CMAKE_Fortran_PREPROCESS ON)

file(GLOB RMTA_SRC CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/*.f90)
list(REMOVE_ITEM RMTA_SRC ${CMAKE_SOURCE_DIR}/src/rmta.f90)

add_library(rmta STATIC ${RMTA_SRC})
target_include_directories(rmta PUBLIC ${QE_MODDIRS} ${QE_INCS})

add_executable(rmta.x src/rmta.f90)
target_link_libraries(rmta.x PRIVATE rmta
  "$<LINK_GROUP:RESCAN,${QE_ARCHIVES}>" ${QE_LIBS} ${DEP_LIBS})
set_target_properties(rmta.x PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Same helper file the make build writes; the example run.sh scripts source it.
file(GENERATE OUTPUT ${CMAKE_SOURCE_DIR}/paths.sh CONTENT
"export QE_ROOT=${QE_ROOT}
export PPRMTA_ROOT=${CMAKE_SOURCE_DIR}
")
