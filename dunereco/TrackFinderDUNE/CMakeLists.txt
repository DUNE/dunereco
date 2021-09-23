art_make( BASENAME_ONLY
          MODULE_LIBRARIES larreco_RecoAlg
                           lardataalg_DetectorInfo
                           lardataobj_RecoBase
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_Utilities
                           larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
                           larsim_MCCheater_BackTrackerService_service
                           larsim_MCCheater_ParticleInventoryService_service
                           nusimdata_SimulationBase
                           ART_FRAMEWORK_CORE
                           ART_FRAMEWORK_PRINCIPAL
                           ART_FRAMEWORK_SERVICES_REGISTRY
                           ART_ROOT_IO_TFILE_SUPPORT
                           ROOT::Core
                           ART_ROOT_IO_TFILESERVICE_SERVICE
                           ART_PERSISTENCY_PROVENANCE
                           art_Persistency_Common
                           art_Utilities canvas
                           messagefacility::MF_MessageLogger
                           FHICLCPP
                           cetlib cetlib_except
                           ROOT_BASIC_LIB_LIST
         )

install_headers()
install_fhicl()
install_source()
