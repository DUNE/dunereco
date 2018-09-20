
art_make( BASENAME_ONLY
	  MODULE_LIBRARIES larreco_RecoAlg
	  		   lardataobj_RecoBase
                           larcorealg_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_Utilities
			   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			   larsim_MCCheater_BackTrackerService_service
         larsim_MCCheater_ParticleInventoryService_service
			   nusimdata_SimulationBase
 			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   art_Persistency_Common canvas
			   art_Utilities canvas
			   ${MF_MESSAGELOGGER}
			   
	                   ${FHICLCPP}
			   cetlib cetlib_except
                           ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
