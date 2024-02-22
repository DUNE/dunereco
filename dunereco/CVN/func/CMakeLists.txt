#for classes_def.xml

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make( BASENAME_ONLY
#  LIBRARY_NAME     CVNFunc
  LIB_LIBRARIES    nusimdata::SimulationBase
  larcorealg::Geometry
  larcorealg::GeoAlgo
  lardataobj::RecoBase
  larsim::MCCheater_BackTrackerService_service 
  larsim::MCCheater_ParticleInventoryService_service
  art::Persistency_Provenance
  art::Framework_Principal
  messagefacility::MF_MessageLogger
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  cetlib::cetlib
  canvas::canvas
  Boost::filesystem            
  ROOT::Hist  
  DICT_LIBRARIES   lardataobj::RecoBase
  dunereco_CVN_func
  ) ### MIGRATE ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()

