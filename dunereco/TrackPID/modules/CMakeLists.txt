include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( 
    MODULE_LIBRARIES larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation
    lardataobj::Simulation
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::AnalysisBase
    lardata::Utilities
    nusimdata::SimulationBase
    larsim::MCCheater_BackTrackerService_service           
    larsim::MCCheater_ParticleInventoryService_service
    ROOT::Core ROOT::Tree
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities
    canvas::canvas
    messagefacility::MF_MessageLogger
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    cetlib_except::cetlib_except
    dunereco_AnaUtils
    dunereco_TrackPID_tf
    dunereco_TrackPID_algorithms
    dunereco_TrackPID_products
#    MODULE_LIBRARIES dunereco_TrackPID
)

install_headers()
install_fhicl()
install_source()

