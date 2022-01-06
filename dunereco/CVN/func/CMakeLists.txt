#for classes_def.xml

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( BASENAME_ONLY
#  LIBRARY_NAME     CVNFunc
  LIB_LIBRARIES    nusimdata::SimulationBase
  larcorealg_Geometry
  larcorealg_GeoAlgo
  lardataobj_RecoBase
  larsim_MCCheater_BackTrackerService_service 
  larsim_MCCheater_ParticleInventoryService_service
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  art::Persistency_Provenance
  messagefacility::MF_MessageLogger
  art::Persistency_Common canvas
  art::Persistency_Provenance canvas
  art::Utilities canvas
  cetlib::cetlib cetlib_except
  canvas::canvas
  Boost::filesystem            
  
  ROOT_BASIC_LIB_LIST
  DICT_LIBRARIES   lardataobj_RecoBase
  dunereco_CVN_func
  )

install_headers()
install_source()

