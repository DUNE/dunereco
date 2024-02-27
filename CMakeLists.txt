# ======================================================================
#  dunereco main build file
# ======================================================================


# use cmake 3.14 or later
cmake_minimum_required (VERSION 3.20 FATAL_ERROR) ### MIGRATE-CHANGED (migrate-3.22.02) - VERSION 3.12...3.14 -> 3.20
# cmake_policy(VERSION 3.14) ### MIGRATE-ACTION-COMPATIBILITY (migrate-3.22.02) - Uncomment to preserve compatibility with older CMake versions

find_package(cetmodules REQUIRED) ### MIGRATE-ADDED (migrate-3.22.02) - find_package(cetmodules) must precede project()
project(dunereco LANGUAGES CXX)
set(${PROJECT_NAME}_CMAKE_PROJECT_VERSION_STRING 09.83.01d00)

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

find_package( Boost REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)

cet_report_compiler_flags() ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - add args: REPORT_THRESHOLD VERBOSE

cet_cmake_module_directories(Modules BINARY)

add_subdirectory(Modules)

find_package( Eigen3 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(eigen ...) -> find_package(eigen REQUIRED ...)
# testing use of rstartree ups product
#find_package( rstartree REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(rstartree ...) -> find_package(rstartree REQUIRED ...)
find_package( art REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art ...) -> find_package(art REQUIRED ...)
find_package( art_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(art_root_io ...) -> find_package(art_root_io REQUIRED ...)
find_package( canvas_root_io REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(canvas_root_io ...) -> find_package(canvas_root_io REQUIRED ...)
find_package(ROOT REQUIRED) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
find_package( dunecore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(dunecore ...) -> find_package(dunecore REQUIRED ...)
find_package( larcore REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcore ...) -> find_package(larcore REQUIRED ...)
find_package( larcoreobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcoreobj ...) -> find_package(larcoreobj REQUIRED ...)
find_package( larcorealg REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larcorealg ...) -> find_package(larcorealg REQUIRED ...)
find_package( larreco REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larreco ...) -> find_package(larreco REQUIRED ...)
find_package( larevt REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larevt ...) -> find_package(larevt REQUIRED ...)
find_package( lardataobj REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardataobj ...) -> find_package(lardataobj REQUIRED ...)
find_package( lardataalg REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardataalg ...) -> find_package(lardataalg REQUIRED ...)
find_package( lardata REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lardata ...) -> find_package(lardata REQUIRED ...)
find_package( larrecodnn REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larrecodnn ...) -> find_package(larrecodnn REQUIRED ...)
find_package( larpandora REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larpandora ...) -> find_package(larpandora REQUIRED ...)
find_package( CLHEP REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(clhep ...) -> find_package(clhep REQUIRED ...)
find_package(GENIE REQUIRED EXPORT)
find_package( GSL REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(gsl ...) -> find_package(gsl REQUIRED ...)
find_package( hep_hpc REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(hep_hpc ...) -> find_package(hep_hpc REQUIRED ...)
find_package( Geant4 REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(...) -> find_package(REQUIRED ...)
if(DEFINED ENV{CAFFE_LIB})
  find_package(caffe REQUIRED) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(caffe ...) -> find_package(caffe REQUIRED ...)
endif()
if(DEFINED ENV{LAPACK_LIB})
  find_package( lapack REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(lapack ...) -> find_package(lapack REQUIRED ...)
endif()
find_package(TritonCommon EXPORT QUIET)
find_package(TritonClient EXPORT QUIET)
find_package( larrecodnn REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(larrecodnn ...) -> find_package(larrecodnn REQUIRED ...)
find_package(TensorFlow EXPORT QUIET)
if(DEFINED ENV{LIBTORCH_DIR})
  find_package( libtorch REQUIRED ) ### MIGRATE-CHANGED (migrate-3.22.02) - find_ups_product(libtorch ...) -> find_package(libtorch REQUIRED ...)
  find_package(Torch REQUIRED)
endif()

find_package(HDF5 REQUIRED EXPORT)

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(dunereco)

cet_cmake_config() ### MIGRATE-ADDED (migrate-3.22.02) - required cmd
