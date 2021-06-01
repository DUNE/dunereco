#Caffe stuff
#include_directories( $ENV{CAFFE_INC} )
#include_directories( $ENV{CAFFE_SRC} )
include_directories( $ENV{GSL_INC}/gsl )
#include_directories( $ENV{PROTOBUF_INC} )
#include_directories( $ENV{LEVELDB_INC} )
#include_directories( $ENV{GFLAGS_INC} )
#include_directories( $ENV{GLOG_INC} )

#cet_find_library( CAFFELIB NAMES caffe-d caffe PATHS ENV CAFFE_LIB NO_DEFAULT_PATH )
cet_find_library( GSLLIB NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )
#cet_find_library( PROTOBUFLIB NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
#cet_find_library( LEVELDBLIB NAMES leveldb PATHS ENV LEVELDB_LIB NO_DEFAULT_PATH )
#cet_find_library( GFLAGSLIB NAMES gflags PATHS ENV GFLAGS_LIB NO_DEFAULT_PATH )
#cet_find_library( GLOGLIB NAMES glog PATHS ENV GLOG_LIB NO_DEFAULT_PATH )

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
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ART_PERSISTENCY_COMMON}
  canvas
  art_Persistency_Provenance
  ${ART_UTILITIES}
  ${MF_MESSAGELOGGER}
  z
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  lardataobj_RecoBase
  lardata_Utilities
  nusimdata_SimulationBase
#  larsim_MCCheater_BackTrackerService_service           
  larsim_MCCheater_ParticleInventoryService_service
  #${CAFFELIB}
  ${GSLLIB}
  ${HEPHPCLIB}
  ${HDF5LIB}
  #${PROTOBUFLIB}
  #${LEVELDBLIB}
  #${GFLAGSLIB}
  #${GLOGLIB}                             
  MVAAlg
  #larreco_RecoAlg
#  ${LARRECO_RECOALG_IMAGEPATTERNALGS_TF}
  boost_filesystem
#  stdc++fs
  MODULE_LIBRARIES dune_CVN_func
  dune_CVN_tf
  dune_CVN_art
  stdc++fs
  )


install_headers()
install_fhicl()
install_source()
