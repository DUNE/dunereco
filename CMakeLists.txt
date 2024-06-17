# ======================================================================
#  dunereco main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.14)

find_package(cetmodules 3.12.00 REQUIRED)
project(dunereco LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.90.02d00)

# cetbuildtools contains our cmake modules

# Project variable defaults (_INIT) and overrides.
set(${PROJECT_NAME}_SCRIPTS_DIR_INIT scripts)
set(${PROJECT_NAME}_WP_DIR wire-cell-cfg)
include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
              -Wno-deprecated-declarations
              $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
)

find_package( Boost REQUIRED COMPONENTS system )

cet_report_compiler_flags()

cet_cmake_module_directories(Modules BINARY)

add_subdirectory(Modules)

find_package( Eigen3 REQUIRED )
# testing use of rstartree ups product
#find_package( rstartree REQUIRED )
find_package( art REQUIRED EXPORT )
find_package( art_root_io REQUIRED )
find_package( canvas_root_io REQUIRED )
find_package( canvas REQUIRED EXPORT)
find_package(ROOT REQUIRED)
find_package( dunecore REQUIRED EXPORT )
find_package( dunepdlegacy REQUIRED EXPORT )
find_package( larcore REQUIRED EXPORT )
find_package( larcoreobj REQUIRED EXPORT )
find_package( larcorealg REQUIRED EXPORT )
find_package( larreco REQUIRED EXPORT )
find_package( larevt REQUIRED EXPORT )
find_package( lardataobj REQUIRED EXPORT )
find_package( lardataalg REQUIRED EXPORT )
find_package( lardata REQUIRED EXPORT )
find_package( larrecodnn REQUIRED EXPORT )
find_package( larpandora REQUIRED EXPORT )
find_package( CLHEP REQUIRED )
find_package( nufinder REQUIRED )
find_package(GENIE REQUIRED EXPORT)
find_package( GSL REQUIRED )
find_package( hep_hpc REQUIRED )
find_package( Geant4 REQUIRED )
if(DEFINED ENV{CAFFE_LIB})
  find_package(caffe REQUIRED)
endif()
if(DEFINED ENV{LAPACK_LIB})
  find_package( lapack REQUIRED )
endif()
find_package(TritonCommon EXPORT QUIET)
find_package(TritonClient EXPORT QUIET)
find_package( larrecodnn REQUIRED )
find_package(TensorFlow EXPORT QUIET)
if(DEFINED ENV{LIBTORCH_DIR})
#  find_package( libtorch REQUIRED )
  find_package(Torch REQUIRED)
endif()

find_ups_product( dunepdlegacy )

find_package(HDF5 REQUIRED EXPORT)

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunereco)

cet_cmake_config()
