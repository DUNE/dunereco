
art_make( BASENAME_ONLY
          LIBRARY_NAME     HitFinderDUNE
          LIB_LIBRARIES    RecoBase
	  		   RecoAlg
			   Filters
                           Geometry
                           Geometry_service
                           RecoObjects
                           LArProperties_service
                           DetectorProperties_service
			   BackTracker_service
                           Utilities
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
                           RecoBase
                           RecoBaseArt
                           Geometry
                           Geometry_service
                           LArProperties_service
                           DetectorProperties_service
			   BackTracker_service
                           Utilities
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
