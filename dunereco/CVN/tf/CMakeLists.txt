if( DEFINED ENV{TENSORFLOW_DIR} )
# This is a hack for a broken tensorflow link, which
  # requires -lrt for glibc earlier than 2.17 (i.e. SLF6)
  # We take advantage of OSTYPE defined by cetbuildtools
  if ( "${OSTYPE}" MATCHES "slf6" )
    set(FOR_TENSORFLOW "-lrt")
  endif ()
endif ()

#art_make(
#          LIB_LIBRARIES
#                        ${FHICLCPP}
#                        cetlib cetlib_except
#                        ${TENSORFLOW}
#                        ${PROTOBUF}
#                        ${FOR_TENSORFLOW}
#      #rt #workaround for larsoft issue 17960
#        )

art_make(BASENAME_ONLY
  LIBRARY_NAME      CVNTf
  LIB_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${ART_PERSISTENCY_COMMON}
  pthread
  protobuf
  canvas
  art_Persistency_Provenance
  ${ART_UTILITIES}
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${TENSORFLOW}
  ${FOR_TENSORFLOW}
  ${ROOT_BASIC_LIB_LIST}
  #lardataobj_RecoBase
  #lardata_Utilities
  nusimdata_SimulationBase
  #larsim_MCCheater_BackTrackerService_service           
  #larsim_MCCheater_ParticleInventoryService_service
  MVAAlg
  )


install_headers()
install_fhicl()
install_source()

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories ( $ENV{TENSORFLOW_INC}/eigen )
include_directories ( $ENV{TENSORFLOW_INC} )
include_directories ( $ENV{PROTOBUF_INC} )
