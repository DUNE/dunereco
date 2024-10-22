include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
  LIBRARY_NAME  NeutrinoAngularRecoAlg
  LIB_LIBRARIES 
  larcorealg::Geometry
  dunereco::AnaUtils
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
  art::Framework_Core
  art::Framework_Principal
  art::Persistency_Common 
  art::Framework_Services_Registry
  art_root_io::tfile_support
  ROOT::Core
  art_root_io::TFileService_service
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities canvas::canvas
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except::cetlib_except
  ROOT::TMVA
  ROOT::EG
  ROOT::Minuit
  ROOT::Minuit2
  )


install_headers()
install_fhicl()
install_source()
