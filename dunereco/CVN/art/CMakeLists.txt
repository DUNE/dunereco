include_directories( $ENV{GSL_INC}/gsl )

cet_find_library( GSLLIB NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )

# Tensorflow stuff
if (DEFINED ENV{TENSORFLOW_DIR})
# For some reason we apparently need to care if larreco is checked out.
#if (larreco_not_in_ups)
#set (LARRECO_RECOALG_IMAGEPATTERNALGS_TF larrecodnn_ImagePatternAlgs_Tensorflow_TF)
#else (larreco_not_in_ups)
#cet_find_library( LARRECO_RECOALG_IMAGEPATTERNALGS_TF NAMES larrecodnn_ImagePatternAlgs_Tensorflow_TF PATHS ENV LARRECO_LIB NO_DEFAULT_PATH)
#endif (larreco_not_in_ups)
else(DEFINED ENV{TENSORFLOW_DIR})
set (EXCLUDE_TF TFNetHandler.cxx)
endif (DEFINED ENV{TENSORFLOW_DIR})

# HEP HPC stuff
if (DEFINED ENV{HEP_HPC_DIR} AND DEFINED ENV{HDF5_DIR})
include_directories( $ENV{HEP_HPC_INC} )
include_directories( $ENV{HDF5_INC} )
cet_find_library( HEPHPCLIB NAMES hep_hpc_hdf5 PATHS ENV HEP_HPC_LIB NO_DEFAULT_PATH )
string(TOUPPER ${CMAKE_BUILD_TYPE} BTYPE_UC )
if( ${BTYPE_UC} MATCHES "DEBUG" )
  cet_find_library( HDF5LIB NAMES hdf5_debug PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
else()
  cet_find_library( HDF5LIB NAMES hdf5 PATHS ENV HDF5_LIB NO_DEFAULT_PATH )
endif()
else (DEFINED ENV{HEP_HPC_DIR} AND DEFINED ENV{HDF5_DIR})
set (EXCLUDE_HEPHPC GCNH5_module.cc)
endif (DEFINED ENV{HEP_HPC_DIR} AND DEFINED ENV{HDF5_DIR})



art_make(BASENAME_ONLY
#  LIBRARY_NAME      CVNArt
  EXCLUDE ${EXCLUDE_TF}
  EXCLUDE ${EXCLUDE_HEPHPC}
  LIB_LIBRARIES dune_CVN_func
  dune_CVN_tf
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  ART_PERSISTENCY_COMMON
  canvas::canvas
  art::Persistency_Provenance
  ART_UTILITIES
  messagefacility::MF_MessageLogger
  z
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  lardataobj_RecoBase
  lardata_Utilities
  nusimdata::SimulationBase
  larsim_MCCheater_ParticleInventoryService_service
  GSLLIB
  HEPHPCLIB
  HDF5LIB
  MVAAlg
  boost_filesystem
  MODULE_LIBRARIES dune_CVN_func
  dune_CVN_tf
  dune_CVN_art
  stdc++fs
  )


install_headers()
install_fhicl()
install_source()
