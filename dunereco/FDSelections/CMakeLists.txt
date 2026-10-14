add_subdirectory(pandrizzle)
add_subdirectory(pandizzle)
add_subdirectory(FDSelectionData)

if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  
art_make( LIBRARY_NAME FDSelection
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
#larreco_Deprecated
	  	        larreco::RecoAlg
         			larreco::Calorimetry
	  	        larpandora::LArPandoraInterface
          	        nusimdata::SimulationBase
            larsim::Utils
            NeutrinoEnergyRecoAlg
            PandrizzleAlg
            PandizzleAlg

MODULE_LIBRARIES FDSelection
	  	           art::Persistency_Common
	  	           art::Persistency_Provenance 
	  	           art::Utilities canvas::canvas
        )

    endif ()
    
add_subdirectory(tools)
add_subdirectory(weights)

install_headers()
install_fhicl()
install_source()
install_scripts()
