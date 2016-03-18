
art_make( 
	  LIBRARY_NAME  dune_ClusterFinderDUNE          
          LIB_LIBRARIES larcore_Geometry
	  		dune_RecoAlgDUNE
                        lardata_RecoBase
                        larreco_RecoAlg_ClusterRecoUtil
          MODULE_LIBRARIES dune_RecoAlgDUNE_Cluster3DAlgs
                           larreco_RecoAlg_ClusterRecoUtil			   
                           larreco_ClusterFinder
                           larsim_MCCheater_BackTracker_service
                           larevt_Filters
                           lardata_RecoBase
                           larcore_Geometry
                           larcore_Geometry_Geometry_service
                           SeedService_service # artextensions
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
         )

install_headers()
install_fhicl()
install_source()


