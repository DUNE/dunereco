# ======================================================================
#  dunereco main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.12...3.14 FATAL_ERROR)

project(dunereco LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.58.02d00)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

# Project variable defaults (_INIT) and overrides.
set(${PROJECT_NAME}_SCRIPTS_DIR_INIT scripts)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

find_ups_boost( )

cet_report_compiler_flags()

cet_cmake_module_directories(Modules BINARY)

add_subdirectory(Modules)

find_ups_product( eigen )
include_directories( $ENV{EIGEN_INC} )
# testing use of rstartree ups product
find_ups_product( rstartree )
find_ups_product( art )
find_ups_product( art_root_io )
find_ups_product( canvas_root_io )
find_ups_root()
find_ups_product( cetbuildtools )
find_ups_product( dunecore )
find_ups_product( larcore )
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larreco )
find_ups_product( larevt )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata )
find_ups_product( larrecodnn )
find_ups_product( larpandora )
find_ups_product( clhep )
find_package(GENIE REQUIRED EXPORT)
find_ups_product( gsl )
find_ups_product( hep_hpc )
find_ups_product( hdf5 )
find_ups_geant4( )
if(DEFINED ENV{CAFFE_LIB} )
  find_ups_product(caffe)
endif()
if(DEFINED ENV{LAPACK_LIB} )
  find_ups_product( lapack )
endif()
find_package(TRITON QUIET EXPORT)
find_ups_product( larrecodnn )
find_package(TensorFlow 2.6.0 QUIET EXPORT)
find_ups_product( libtorch )
find_package(Torch REQUIRED)


# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunereco)
add_subdirectory(ups)

include(UseCPack)
