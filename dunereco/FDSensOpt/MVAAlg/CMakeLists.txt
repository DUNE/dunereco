include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
  LIBRARY_NAME  MVAAlg
  LIB_LIBRARIES larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
  lardataalg_DetectorInfo
  lardataobj_RawData
  lardataobj_RecoBase
  lardata_ArtDataHelper
  larreco_RecoAlg
  larreco_RecoAlg_PMAlg
  lardata_Utilities
  larreco_Calorimetry
  IniSegAlg
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  nusimdata_SimulationBase
  ART_FRAMEWORK_CORE
  ART_FRAMEWORK_PRINCIPAL
  art_Persistency_Common canvas
  ART_FRAMEWORK_SERVICES_REGISTRY
  ART_ROOT_IO_TFILE_SUPPORT
  ROOT::Core
  ART_ROOT_IO_TFILESERVICE_SERVICE
  art_Persistency_Common canvas
  art_Persistency_Provenance canvas
  art_Utilities canvas
  FHICLCPP
  messagefacility::MF_MessageLogger
  
  cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  ${Boost_SYSTEM_LIBRARY}
  ROOT_TMVA
  )


install_headers()
install_fhicl()
install_source()
