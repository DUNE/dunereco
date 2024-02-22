art_make( 
    LIB_LIBRARIES larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation lardataobj::Simulation
    lardataalg::DetectorInfo
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::AnalysisBase
    lardata::Utilities
    nusimdata::SimulationBase
    ROOT::Core
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    dunereco_CVN_func
    dunereco_TrackPID_products
) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

