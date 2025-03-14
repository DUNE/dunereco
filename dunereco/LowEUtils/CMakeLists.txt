# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

art_make(BASENAME_ONLY ALLOW_UNDERSCORES
  LIBRARY_NAME LowEUtils
  LIB_LIBRARIES
    fhiclcpp::fhiclcpp
    lardata::Utilities
    lardataobj::RecoBase
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_PhotonBackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
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
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art::Persistency_Common
    art::Utilities 
    art_root_io::tfile_support
    art_root_io::TFileService_service
    Boost::filesystem
    canvas::canvas
    cetlib::cetlib cetlib_except::cetlib_except
    CLHEP::CLHEP
    dunecore::ProducerUtils
    duneopdet::OpticalDetector
    fhiclcpp::fhiclcpp
    larcore::Geometry_Geometry_service
    larcorealg::Geometry
    lardata::headers
    lardataalg::DetectorInfo
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_PhotonBackTrackerService_service
    larsim::MCCheater_ParticleInventoryService_service
    messagefacility::MF_MessageLogger
    ROOT::Core
  MODULE_LIBRARIES LowEUtils
    dunecore::DuneObj
    dunecore::ProducerUtils
    lardataalg::DetectorInfo
    lardataobj::RecoBase
    lardata::headers
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    larsim::MCCheater_BackTrackerService_service
    larsim::MCCheater_PhotonBackTrackerService_service
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
)

add_subdirectory(fcl)

install_headers()
install_fhicl()
install_scripts()
install_source()