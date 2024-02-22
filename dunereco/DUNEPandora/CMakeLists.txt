# where should the scripts/..xml file be installed?  Perhaps in bin?

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible


art_make( BASENAME_ONLY
          LIBRARY_NAME  DUNEPandora
          LIB_LIBRARIES larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			larpandoracontent::LArPandoraContent
                        PandoraPFA::PandoraSDK
                        PandoraPFA::PandoraMonitoring
                        nusimdata::SimulationBase
                        ROOT::Core
                        art::Persistency_Common canvas
                        art::Persistency_Provenance canvas
                        art::Utilities canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
                        ROOT::Geom
			
          MODULE_LIBRARIES DUNEPandora
          ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

add_subdirectory(scripts)
