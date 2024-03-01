art_make( BASENAME_ONLY
          MODULE_LIBRARIES larreco::RecoAlg
                           lardataalg::DetectorInfo
                           lardataobj::RecoBase
                           larcorealg::Geometry
                           larcore::Geometry_Geometry_service
                           lardata::Utilities
                           larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
                           larsim::MCCheater_BackTrackerService_service
                           larsim::MCCheater_ParticleInventoryService_service
                           nusimdata::SimulationBase
                           ROOT::Core
			   ROOT::Hist
                           art::Persistency_Common
                           art::Utilities
			   art_root_io::tfile_support
                           messagefacility::MF_MessageLogger
                           fhiclcpp::fhiclcpp
                           cetlib::cetlib
         )

install_headers()
install_fhicl()
install_source()
