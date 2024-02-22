include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make( 
    MODULE_LIBRARIES larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation lardataobj_Simulation
    lardataobj_RawData
    lardataobj_RecoBase
    lardataobj_AnalysisBase
    lardata_Utilities
    nusimdata::SimulationBase
    larsim_MCCheater_BackTrackerService_service           
    larsim_MCCheater_ParticleInventoryService_service
    ROOT::Core
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib cetlib_except
    dunereco_AnaUtils
    dunereco_TrackPID_tf
    dunereco_TrackPID_algorithms
    dunereco_TrackPID_products
#    MODULE_LIBRARIES dunereco_TrackPID
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

