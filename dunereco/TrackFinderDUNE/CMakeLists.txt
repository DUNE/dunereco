
art_make( BASENAME_ONLY
	  MODULE_LIBRARIES larreco_RecoAlg
	  		   lardataobj_RecoBase
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_Utilities
			   larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
			   larsim_MCCheater_BackTracker_service
			   nusimdata_SimulationBase
 			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   art_Persistency_Common canvas_Persistency_Common
			   art_Utilities canvas_Utilities
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	                   ${FHICLCPP}
			   ${CETLIB}
                           ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
