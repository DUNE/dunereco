
art_make( BASENAME_ONLY
          LIBRARY_NAME     HitFinderDUNE
          LIB_LIBRARIES    lardataobj::RecoBase
                           lardataalg::DetectorInfo
                           larreco::RecoAlg
                           larreco::RecoAlg_PMAlg
                           larevt::Filters
                           larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           lardata::RecoObjects
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           nusimdata::SimulationBase
                           ROOT::Core
			   ROOT::Hist
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities
			   art::Framework_Services_Registry
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
                           CLHEP::CLHEP
         MODULE_LIBRARIES  HitFinderDUNE
                           lardataobj::RecoBase
                           lardata::ArtDataHelper
                           larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           lardata::Utilities
                           nurandom::RandomUtils_NuRandomService_service
                           nusimdata::SimulationBase
                           ROOT::Core
                           art::Persistency_Common 
                           art::Utilities 
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
         )

install_headers()
install_fhicl()
install_source()
