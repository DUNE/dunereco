
art_make(MODULE_LIBRARIES
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  ART_PERSISTENCY_COMMON
  dune_AnaUtils
  canvas::canvas
  art::Persistency_Provenance
  ART_UTILITIES
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  lardataobj_RecoBase
  lardata_Utilities
  nusimdata::SimulationBase
  larsim_MCCheater_BackTrackerService_service           
  larsim_MCCheater_ParticleInventoryService_service
  larreco_RecoAlg
  larreco_Calorimetry
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larpandora_LArPandoraInterface
  )


install_headers()
install_fhicl()
install_source()
