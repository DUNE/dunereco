art_make( 
    LIB_LIBRARIES larcorealg_Geometry
    larcore_Geometry_Geometry_service
    larsim_Simulation lardataobj_Simulation
    lardataalg_DetectorInfo
    lardataobj_RawData
    lardataobj_RecoBase
    lardataobj_AnalysisBase
    lardata_Utilities
    nusimdata_SimulationBase
    ART_FRAMEWORK_CORE
    ART_FRAMEWORK_PRINCIPAL
    ART_FRAMEWORK_SERVICES_REGISTRY
    ART_ROOT_IO_TFILE_SUPPORT
    ROOT::Core
    ART_ROOT_IO_TFILESERVICE_SERVICE
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    messagefacility::MF_MessageLogger
    FHICLCPP
    cetlib cetlib_except
    ${Boost_SYSTEM_LIBRARY}
    ROOT_BASIC_LIB_LIST
    dune_CVN_func
    dune_TrackPID_products
)

install_headers()
install_fhicl()
install_source()

