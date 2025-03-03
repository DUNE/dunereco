include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
  LIBRARY_NAME  ParticleSelectionAlg
  LIB_LIBRARIES 
  larcorealg::Geometry
  dunereco_AnaUtils
  larcore::Geometry_Geometry_service
  larsim::Simulation 
  lardataobj::RawData
  lardataobj::RecoBase
  lardata::ArtDataHelper
  larreco::RecoAlg
  larpandora::LArPandoraInterface
  lardata::Utilities
  ROOT::Core
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  cetlib::cetlib 
  )


install_headers()
install_fhicl()
install_source()
