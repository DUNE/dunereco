add_subdirectory(IniSegAlg)
add_subdirectory(MVAAlg)

include_directories ( $ENV{GENIE_INC}/GENIE )

simple_plugin(NueAna   "module"
		        lardataobj_RawData
		        lardataobj_RecoBase
		        lardata_RecoBaseArt
			lardata_AnalysisAlg
		        larreco_RecoAlg_PMAlg
		        IniSegAlg
		        larcore_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation larsimobj_Simulation
			larsim_MCCheater_BackTracker_service
			nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(MVASelect   "module"
                        MVAAlg
                        nutools_NuReweight
                        nutools_NuReweight_art
			nusimdata_SimulationBase
                        ${GHEP}
                        ${GEVGCORE}
			${GPDG}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(IniSegReco "module"
		      lardataobj_RawData
		      lardataobj_RecoBase
		      larreco_RecoAlg
		      larreco_RecoAlg_PMAlg
		      IniSegAlg
		      larcore_Geometry
		      larcore_Geometry_Geometry_service
		      lardata_Utilities
		      larsim_Simulation larsimobj_Simulation
		      larsim_MCCheater_BackTracker_service
		      nusimdata_SimulationBase
      		      ${ART_FRAMEWORK_CORE}
		      ${ART_FRAMEWORK_PRINCIPAL}
		      ${ART_FRAMEWORK_SERVICES_REGISTRY}
      		      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
      		      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 		      art_Persistency_Common canvas_Persistency_Common
		      art_Persistency_Provenance canvas_Persistency_Provenance
		      art_Utilities canvas_Utilities
		      ${MF_MESSAGELOGGER}
		      ${MF_UTILITIES}
		      ${CETLIB}
		      ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(ShSeg   "module"
			lardataobj_RawData
			lardataobj_RecoBase
			larreco_RecoAlg
			lardata_AnalysisAlg
			larreco_RecoAlg_PMAlg
			IniSegAlg
			larcore_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			larsim_Simulation larsimobj_Simulation
			
			
			larsim_MCCheater_BackTracker_service
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
			BASENAME_ONLY
)



install_headers()
install_fhicl()
install_source()
install_scripts()
