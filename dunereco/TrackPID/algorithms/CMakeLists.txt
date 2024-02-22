include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make( 
    LIB_LIBRARIES larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation lardataobj_Simulation
    lardataobj_RawData
    lardataobj_RecoBase
    lardataobj_AnalysisBase
    lardata_Utilities
    nusimdata::SimulationBase
    larsim_MCCheater_BackTrackerService_service           
    larsim_MCCheater_ParticleInventoryService_service
    ART_FRAMEWORK_CORE
    ART_FRAMEWORK_PRINCIPAL
    ART_FRAMEWORK_SERVICES_REGISTRY
    ART_ROOT_IO_TFILE_SUPPORT
    ROOT::Core
    ART_ROOT_IO_TFILESERVICE_SERVICE
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib cetlib_except
    dunereco_AnaUtils
    dunereco_TrackPID_tf
    dunereco_TrackPID_products
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

