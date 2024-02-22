art_make( LIBRARY_NAME PandrizzleAlg
          LIB_LIBRARIES larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
                        larsim::MCCheater_BackTrackerService_service
                        larsim::MCCheater_ParticleInventoryService_service
                        lardata::Utilities
                        larevt::Filters
                        lardataobj::RawData
                        lardataobj::RecoBase
                        lardata::RecoObjects 
                        ROOT::Tree
                        larreco::RecoAlg
                        larreco::Calorimetry
                        larpandora::LArPandoraInterface
                        nusimdata::SimulationBase
                        larsim::Utils
                        NeutrinoEnergyRecoAlg

       MODULE_LIBRARIES PandrizzleAlg
                        dune_TrackPID_algorithms
                        dune_TrackPID_products
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
install_scripts()
