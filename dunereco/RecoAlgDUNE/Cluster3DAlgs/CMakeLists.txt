
art_make( 
	  LIBRARY_NAME  dunereco_RecoAlgDUNE_Cluster3DAlgs
          LIB_LIBRARIES larreco_VertexFinder_HarrisVertexFinder_module
                        lardataobj_RecoBase
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        larevt_Filters
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
                        lardataobj_AnalysisBase
                        lardata_Utilities
			ROOT::Core
			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
			CLHEP::CLHEP
                        Boost::filesystem
                        
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
