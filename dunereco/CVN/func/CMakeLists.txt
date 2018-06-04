#for classes_def.xml

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( BASENAME_ONLY
  LIBRARY_NAME     CVNFunc
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
  DICT_LIBRARIES   lardataobj_RecoBase
  CVNFunc
  )

install_headers()
install_source()

