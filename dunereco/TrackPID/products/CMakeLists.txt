include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make( 
    LIB_LIBRARIES 
#    larcorealg_Geometry
#    larcore_Geometry_Geometry_service
#    larsim_Simulation lardataobj_Simulation
#    lardataobj_RawData
#    lardataobj_RecoBase
#    lardataobj_AnalysisBase
#    lardata_Utilities
#    nusimdata_SimulationBase
#    larsim_MCCheater_BackTrackerService_service           
#    larsim_MCCheater_ParticleInventoryService_service
#    ${ART_FRAMEWORK_CORE}
#    ${ART_FRAMEWORK_PRINCIPAL}
#    ${ART_FRAMEWORK_SERVICES_REGISTRY}
#    ${ART_FRAMEWORK_SERVICES_BASIC}
#    ${ART_ROOT_IO_TFILE_SUPPORT} ROOT::Core
#    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
#    art_Persistency_Common
#    art_Persistency_Provenance
#    art_Utilities
#    canvas
#    ${messagefacility::MF_MessageLogger}
#    ${FHICLCPP}
#    cetlib cetlib_except
#    ${Boost_SYSTEM_LIBRARY}
#    ROOT_BASIC_LIB_LIST
)

install_headers()
install_fhicl()
install_source()

