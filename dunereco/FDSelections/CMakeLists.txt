add_subdirectory(pandrizzle)
add_subdirectory(pandizzle)
add_subdirectory(FDSelectionData)

if(( DEFINED ENV{TENSORFLOW_DIR} ) AND ( DEFINED ENV{TRITON_DIR} ))
  
art_make( LIBRARY_NAME FDSelection
	  LIB_LIBRARIES larcorealg_Geometry
	  	        larcore_Geometry_Geometry_service
		          larsim_Simulation nug4_ParticleNavigation lardataobj_Simulation
 	  	       larsim_MCCheater_BackTrackerService_service
         	  larsim_MCCheater_ParticleInventoryService_service
	  	        lardata_Utilities
	  	        larevt_Filters
	  	        lardataobj_RawData
	  	        lardataobj_RecoBase
	  	        lardata_RecoObjects 
#larreco_Deprecated
	  	        larreco_RecoAlg
         			larreco_Calorimetry
	  	        larpandora_LArPandoraInterface
          	        nusimdata_SimulationBase
          	        ${LARRECO_LIB}
          	        ${LARDATA_LIB}
	  	        ${MF_MESSAGELOGGER}
	  	        ${MF_UTILITIES}
	  	        ${FHICLCPP}
	  	        cetlib cetlib_except
          	        ${ROOT_GEOM}
          	        ${ROOT_XMLIO}
          	        ${ROOT_GDML}
	  	        ${ROOT_BASIC_LIB_LIST}
            ${ROOT_TMVA}
            larsim_Utils
            NeutrinoEnergyRecoAlg
            PandrizzleAlg
            PandizzleAlg

MODULE_LIBRARIES FDSelection
          	           ${ART_FRAMEWORK_CORE}
	  	           ${ART_FRAMEWORK_PRINCIPAL}
		       	   ${ART_FRAMEWORK_SERVICES_REGISTRY}
	  	       	   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	  	           art_Persistency_Common canvas
	  	           art_Persistency_Provenance canvas
	  	           art_Utilities canvas
                   dunereco_TrackPID_algorithms
        )

    endif ()
    
add_subdirectory(tools)
add_subdirectory(weights)

install_headers()
install_fhicl()
install_source()
install_scripts()
