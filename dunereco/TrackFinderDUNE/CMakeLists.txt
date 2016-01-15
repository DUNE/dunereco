
art_make( BASENAME_ONLY
	  MODULE_LIBRARIES larreco_RecoAlg
	  		   lardata_RecoBase
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
                           lardata_Utilities
			   larsim_Simulation
			   larsim_MCCheater_BackTracker_service
			   ${SIMULATIONBASE}
 			   ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			   ${ART_PERSISTENCY_COMMON}
			   ${ART_UTILITIES}
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	                   ${FHICLCPP}
			   ${CETLIB}
                           ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
