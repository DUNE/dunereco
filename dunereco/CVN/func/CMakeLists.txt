#for classes_def.xml

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( BASENAME_ONLY
#  LIBRARY_NAME     CVNFunc
  LIB_LIBRARIES    nusimdata_SimulationBase
  larcorealg_Geometry
  larcorealg_GeoAlgo
  lardataobj_RecoBase
  larsim_MCCheater_BackTrackerService_service 
  larsim_MCCheater_ParticleInventoryService_service
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  art_Persistency_Provenance
  messagefacility::MF_MessageLogger
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  cetlib cetlib_except
  canvas
  ${Boost_FILESYSTEM_LIBRARY}            
  ${Boost_SYSTEM_LIBRARY}
  ROOT_BASIC_LIB_LIST
  DICT_LIBRARIES   lardataobj_RecoBase
  dune_CVN_func
  )

install_headers()
install_source()

