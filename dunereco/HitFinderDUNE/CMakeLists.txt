
art_make( BASENAME_ONLY
          LIBRARY_NAME     HitFinderDUNE
          LIB_LIBRARIES    lardata_RecoBase
	  		   larreco_RecoAlg
			   larreco_RecoAlg_PMAlg
			   larevt_Filters
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_RecoObjects
                           lardata_Utilities_LArProperties_service
                           lardata_Utilities_DetectorProperties_service
			   larsim_MCCheater_BackTracker_service
                           lardata_Utilities
			   ${SIMULATIONBASE}
                           ${ART_FRAMEWORK_CORE}
   			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			   ${ART_PERSISTENCY_COMMON}
		       	   ${ART_PERSISTENCY_PROVENANCE}
			   ${ART_UTILITIES}
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	                   ${FHICLCPP}
			   ${CETLIB}
                           ${CLHEP}
			   ${Boost_SYSTEM_LIBRARY}
			   ${ROOT_BASIC_LIB_LIST}
         MODULE_LIBRARIES HitFinderDUNE
                           lardata_RecoBase
                           lardata_RecoBaseArt
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_Utilities_LArProperties_service
                           lardata_Utilities_DetectorProperties_service
			   larsim_MCCheater_BackTracker_service
                           lardata_Utilities
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
