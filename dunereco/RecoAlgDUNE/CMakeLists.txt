
art_make(
	  LIBRARY_NAME  dune_RecoAlgDUNE
          LIB_LIBRARIES larsim_MCCheater_BackTrackerService_service
           larsim_MCCheater_ParticleInventoryService_service
                        lardataobj_RecoBase
                        lardata_ArtDataHelper
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        larevt_Filters
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_RecoObjects larreco_Deprecated
                        lardataobj_AnalysisBase
			larreco_Calorimetry
                        lardata_Utilities
                        larreco_RecoAlg_ClusterRecoUtil
			larreco_RecoAlg_CMTool_CMToolBase
			larreco_RecoAlg_CMTool_CMTAlgMerge
			larreco_RecoAlg_CMTool_CMTAlgPriority
			larreco_RecoAlg_CMTool_CMTAlgMatch
			larreco_RecoAlg_CMTool_CMToolApp
			larreco_RecoAlg_PMAlg
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			${MF_MESSAGELOGGER}
			
	                ${FHICLCPP}
			cetlib cetlib_except
                        ${CLHEP}
			${ROOT_BASIC_LIB_LIST}
                        ${ROOT_MINUIT}
			${ROOT_MINUIT2}	
                        ${Boost_FILESYSTEM_LIBRARY}
                        ${Boost_SYSTEM_LIBRARY}
        )

#add_subdirectory(Cluster3DAlgs)

install_headers()
install_fhicl()
install_source()

