include_directories(${CMAKE_CURRENT_SOURCE_DIR})

art_make( BASENAME_ONLY
          LIBRARY_NAME  IniSegAlg
          LIB_LIBRARIES larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
                        lardataobj_RawData
                        lardataobj_RecoBase
			larreco_RecoAlg
			larreco_RecoAlg_PMAlg
                        lardata_Utilities
                        nusimdata::SimulationBase
                        ART_FRAMEWORK_CORE
                        ART_FRAMEWORK_PRINCIPAL
                        art::Persistency_Common canvas
                        ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        fhiclcpp::fhiclcpp
			messagefacility::MF_MessageLogger
			
                        cetlib::cetlib cetlib_except
                        ROOT_BASIC_LIB_LIST
			
          )


install_headers()
install_fhicl()
install_source()
