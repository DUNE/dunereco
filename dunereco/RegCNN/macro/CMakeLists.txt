
art_make(MODULE_LIBRARIES
  ROOT::Core
  dunereco_AnaUtils
  canvas::canvas
  art::Persistency_Provenance
  art_root_io::tfile_support
  ROOT::Tree
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib
  lardataobj::RecoBase
  lardata::Utilities
  nusimdata::SimulationBase
  larsim::MCCheater_BackTrackerService_service           
  larsim::MCCheater_ParticleInventoryService_service
  larreco::RecoAlg
  larreco::Calorimetry
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larpandora::LArPandoraInterface
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


install_headers()
install_fhicl()
install_source()
