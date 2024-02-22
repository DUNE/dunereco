include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

art_make( 
    LIB_LIBRARIES 
#    larcorealg_Geometry
#    larcore_Geometry_Geometry_service
#    larsim_Simulation lardataobj_Simulation
#    lardataobj_RawData
    lardataobj::RecoBase
#    lardataobj_AnalysisBase
#    lardata_Utilities
#    nusimdata::SimulationBase
#    larsim_MCCheater_BackTrackerService_service           
#    larsim_MCCheater_ParticleInventoryService_service
#    art::Framework_Core
#    art::Framework_Principal
#    art::Framework_Services_Registry
#    art_root_io::tfile_support ROOT::Core
#    art_root_io::TFileService_service
#    art::Persistency_Common
#    art::Persistency_Provenance
#    art::Utilities
#    canvas::canvas
#    messagefacility::MF_MessageLogger
#    fhiclcpp::fhiclcpp
#    cetlib::cetlib cetlib_except
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

