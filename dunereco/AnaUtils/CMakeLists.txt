art_make( 
    LIB_LIBRARIES larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation lardataobj_Simulation
    lardataalg_DetectorInfo
    lardataobj_RawData
    lardataobj_RecoBase
    lardataobj_AnalysisBase
    lardata_Utilities
    nusimdata::SimulationBase
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
    
    ROOT_BASIC_LIB_LIST
    dune_CVN_func
    dune_TrackPID_products
)

install_headers()
install_fhicl()
install_source()

