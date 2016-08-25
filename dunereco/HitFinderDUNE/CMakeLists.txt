
art_make( BASENAME_ONLY
          LIBRARY_NAME     HitFinderDUNE
          LIB_LIBRARIES    lardataobj_RecoBase
	  		   larreco_RecoAlg
			   larreco_RecoAlg_PMAlg
			   larevt_Filters
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           lardata_RecoObjects
			   larsim_MCCheater_BackTracker_service
                           lardata_Utilities
			   lbne_daqinput35t
			   nusimdata_SimulationBase
                           ${ART_FRAMEWORK_CORE}
   			   ${ART_FRAMEWORK_PRINCIPAL}
			   ${ART_FRAMEWORK_SERVICES_REGISTRY}
			   ${ART_FRAMEWORK_SERVICES_BASIC}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			   ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                           ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			   art_Persistency_Common canvas_Persistency_Common
		       	   art_Persistency_Provenance canvas_Persistency_Provenance
			   art_Utilities canvas_Utilities
			   ${MF_MESSAGELOGGER}
			   ${MF_UTILITIES}
	                   ${FHICLCPP}
			   ${CETLIB}
                           ${CLHEP}
			   ${Boost_SYSTEM_LIBRARY}
			   ${ROOT_BASIC_LIB_LIST}
         MODULE_LIBRARIES HitFinderDUNE
                           lardataobj_RecoBase
                           lardata_RecoBaseArt
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
			   larsim_MCCheater_BackTracker_service
                           lardata_Utilities
			   larsim_RandomUtils_LArSeedService_service
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
