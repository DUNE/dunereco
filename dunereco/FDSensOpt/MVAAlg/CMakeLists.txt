include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
  LIBRARY_NAME  MVAAlg
  LIB_LIBRARIES larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
  lardataalg::DetectorInfo
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
  art::Persistency_Common
  ROOT::Core
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  art::Framework_Principal
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  cetlib::cetlib
  )


install_headers()
install_fhicl()
install_source()
