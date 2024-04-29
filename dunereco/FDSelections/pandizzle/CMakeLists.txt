art_make(LIBRARY_NAME     PandizzleAlg
         LIB_LIBRARIES    larcorealg::Geometry
                          larcore::Geometry_Geometry_service
                          larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
                          larsim::MCCheater_BackTrackerService_service
                          larsim::MCCheater_ParticleInventoryService_service
                          lardata::Utilities
                          larevt::Filters
                          lardataobj::RawData
                          lardataobj::RecoBase
                          lardata::RecoObjects 
                          larreco::RecoAlg
                          larreco::Calorimetry
                          larpandora::LArPandoraInterface
                          nusimdata::SimulationBase
                          larsim::Utils
                          NeutrinoEnergyRecoAlg
         MODULE_LIBRARIES PandizzleAlg
                          #dune_TrackPID
                          dune_TrackPID_algorithms
                          dune_TrackPID_products
                          art::Persistency_Common
                          art::Persistency_Provenance
                          art::Utilities
        )

install_headers()
install_fhicl()
install_source()
install_scripts()
