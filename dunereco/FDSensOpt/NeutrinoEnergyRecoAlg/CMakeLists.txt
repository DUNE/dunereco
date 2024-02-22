include_directories(${CMAKE_CURRENT_SOURCE_DIR}) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make( BASENAME_ONLY
  LIBRARY_NAME  NeutrinoEnergyRecoAlg
  LIB_LIBRARIES 
  larcorealg::Geometry
  dunereco_AnaUtils
  larcore::Geometry_Geometry_service
  larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
  lardataobj::RawData
  lardataobj::RecoBase
  lardata::ArtDataHelper
  larreco::RecoAlg
  larreco::RecoAlg_PMAlg
  lardata::Utilities
  larreco::Calorimetry
  IniSegAlg
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  nusimdata::SimulationBase
  ROOT::Core
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  cetlib::cetlib 
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


install_headers()
install_fhicl()
install_source()
