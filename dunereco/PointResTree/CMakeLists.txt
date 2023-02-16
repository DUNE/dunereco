
simple_plugin(PointResTree "module"
		larcore_Geometry_Geometry_service
		larcorealg_Geometry
		larsim_Simulation
		lardataobj_RawData
		lardataobj_RecoBase
		lardataalg_DetectorInfo
		larreco_Calorimetry
		larsim_MCCheater_BackTrackerService_service
		larsim_MCCheater_ParticleInventoryService_service
		nusimdata_SimulationBase
		${ART_FRAMEWORK_CORE}
		${ART_FRAMEWORK_PRINCIPAL}
		${ART_FRAMEWORK_SERVICES_REGISTRY}
		${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
		${ART_ROOT_IO_TFILESERVICE_SERVICE}
		art_Persistency_Common canvas
		art_Persistency_Provenance canvas
		art_Utilities canvas
		${MF_MESSAGELOGGER}
		${MF_UTILITIES}
		${CETLIB}
		${ROOT_BASIC_LIB_LIST}
        dune_AnaUtils
)

install_headers()
install_fhicl()
install_source()
install_scripts()CMakeLists.cc 
