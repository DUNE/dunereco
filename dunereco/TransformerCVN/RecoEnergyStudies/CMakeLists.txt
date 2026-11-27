# dunereco::CVN_art (PixelMapProducer, used only for the vertical-drift global wire mapping) is built only with
# TensorFlow and Triton (dunereco/CVN/CMakeLists.txt). Without them the module is built for horizontal drift only.
simple_plugin(RecoEnergyS "module"
    lardataobj::RawData
    lardataobj::RecoBase
    lardata::RecoObjects
    larcorealg::Geometry
    larreco::Calorimetry
    larreco::RecoAlg
    lardata::ArtDataHelper
    lardataalg::DetectorInfo
    larcore::Geometry_Geometry_service
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
    lardata::Utilities
    larpandora::LArPandoraInterface
    nusimdata::SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
    ${ART_ROOT_IO_TFILESERVICE_SERVICE}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    art::Persistency_Common canvas::canvas
    art::Persistency_Provenance canvas::canvas
    art::Utilities canvas::canvas
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    cetlib::cetlib cetlib_except::cetlib_except
    ${ROOT_BASIC_LIB_LIST}
    BASENAME_ONLY
    dunereco::CVN_func
    )
if (( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  target_link_libraries(RecoEnergyS_module PRIVATE dunereco::CVN_art)
  target_compile_definitions(RecoEnergyS_module PRIVATE RECOENERGYS_WITH_CVN_ART)
endif()
install_headers()
install_fhicl()
install_source()
install_scripts()
