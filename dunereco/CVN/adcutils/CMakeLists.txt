simple_plugin(SPMultiTpcDump "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_PointIdAlg
			larreco_RecoAlg_PMAlg
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
 			art_Persistency_Common
			art_Persistency_Provenance
			art_Utilities
			canvas
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

