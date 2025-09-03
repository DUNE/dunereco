art_make(
    LIB_LIBRARIES larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::Simulation lardataobj::Simulation
    lardataalg::DetectorInfo
    lardataobj::RawData
    lardataobj::RecoBase
    lardataobj::AnalysisBase
    lardata::Utilities
    larpandora::LArPandoraInterface
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
)

install_headers()
install_fhicl()
install_source()

