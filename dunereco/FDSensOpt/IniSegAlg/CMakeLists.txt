include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  IniSegAlg
          LIB_LIBRARIES larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
                        lardataobj::RawData
                        lardataobj::RecoBase
			larreco::RecoAlg
			larreco::RecoAlg_PMAlg
                        lardata::Utilities
                        nusimdata::SimulationBase
                        art::Persistency_Common
                        ROOT::Core
                        art::Persistency_Provenance
			art_root_io::tfile_support
                        art::Utilities
                        fhiclcpp::fhiclcpp
			messagefacility::MF_MessageLogger
                        cetlib::cetlib
			
          )


install_headers()
install_fhicl()
install_source()
