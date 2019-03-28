#for classes_def.xml
# not ready for c2 yet
if( DEFINED ENV{TENSORFLOW_DIR} )
  # find tensorflow library
  find_ups_product( protobuf )
  find_ups_product( tensorflow )
  cet_find_library(TENSORFLOW NAMES tensorflow-core PATHS ENV TENSORFLOW_LIB NO_DEFAULT_PATH )
  cet_find_library(PROTOBUF NAMES protobuf PATHS ENV PROTOBUF_LIB NO_DEFAULT_PATH )
  # This is a hack for a broken tensorflow link, which
  # requires -lrt for glibc earlier than 2.17 (i.e. SLF6)
  # We take advantage of OSTYPE defined by cetbuildtools
  if ( "${OSTYPE}" MATCHES "slf6" )
    set(FOR_TENSORFLOW "-lrt")
  endif ()
endif ()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{TENSORFLOW_INC}/eigen )

art_make( BASENAME_ONLY
  LIBRARY_NAME     RegCNNFunc
  LIB_LIBRARIES    nusimdata_SimulationBase
  lardataobj_RecoBase
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${MF_MESSAGELOGGER}
  
  cetlib cetlib_except
  ${Boost_FILESYSTEM_LIBRARY}            
  ${Boost_SYSTEM_LIBRARY}
  ${ROOT_BASIC_LIB_LIST}
  ${TENSORFLOW}
  ${PROTOBUF}
  ${FOR_TENSORFLOW}
  DICT_LIBRARIES   lardataobj_RecoBase
  RegCNNFunc
  )

install_headers()
install_source()

