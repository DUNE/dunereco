# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(BASENAME_ONLY ALLOW_UNDERSCORES
  LIBRARY_NAME LowEUtils
  LIB_LIBRARIES
    fhiclcpp::fhiclcpp
    lardata::Utilities
    lardataobj::RecoBase
    larreco::RecoAlg_Cluster3DAlgs
    larreco::RecoAlg_ClusterRecoUtil
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    messagefacility::MF_MessageLogger
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    ROOT::Core
    art_root_io::TFileService_service
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Common
    art::Utilities 
    canvas::canvas
    cetlib::cetlib 
    cetlib_except::cetlib_except 
    Boost::filesystem
  SERVICE_LIBRARIES
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    lardata::headers
    lardataalg::DetectorInfo
    duneopdet::OpticalDetector
    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    ROOT::Core
    art_root_io::TFileService_service
    art::Persistency_Common
    art::Utilities 
    canvas::canvas
    cetlib::cetlib cetlib_except::cetlib_except
    CLHEP::CLHEP
    Boost::filesystem
  MODULE_LIBRARIES LowEUtils
    dunecore::DuneObj
    lardataalg::DetectorInfo
    lardataobj::RecoBase
    larreco::RecoAlg_Cluster3DAlgs
    larreco::RecoAlg_ClusterRecoUtil
    lardata::headers
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::MCCheater_ParticleInventoryService_service
    nusimdata::SimulationBase
    nurandom::RandomUtils_NuRandomService_service
    fhiclcpp::fhiclcpp
    messagefacility::MF_MessageLogger
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    ROOT::Core
    ROOT::Tree
    ROOT::Graf
    art_root_io::TFileService_service
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Persistency_Provenance
    art::Persistency_Common
    art::Utilities 
    canvas::canvas
    cetlib::cetlib 
    cetlib_except::cetlib_except
    CLHEP::CLHEP
    Boost::filesystem
    FFTW3::FFTW3
  DICT_LIBRARIES
    larreco::RecoAlg_Cluster3DAlgs
)

install_headers()
install_fhicl()
install_scripts()
install_source()