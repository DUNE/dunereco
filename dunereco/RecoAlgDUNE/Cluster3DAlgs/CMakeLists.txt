
art_make( 
	  LIBRARY_NAME  dune_RecoAlgDUNE_Cluster3DAlgs
          LIB_LIBRARIES larreco_VertexFinder_HarrisVertexFinder_module
                        lardataobj_RecoBase
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        larevt_Filters
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects
                        lardataobj_AnalysisBase
                        lardata_Utilities
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
                        ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE
			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        CLHEP
			ROOT_BASIC_LIB_LIST
                        ROOT_MINUIT
			ROOT_MINUIT2	
                        Boost::filesystem
                        
        )

install_headers()
install_fhicl()
install_source()
