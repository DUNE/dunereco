art_make_library( 
  LIBRARY_NAME     GlobImage
  SOURCE EventImageData.h EventImageData.cxx
  LIBRARIES    CVNFunc
  LIBRARY_NAME_VAR GLOBIMAGE
#  lardataobj_RecoBase
#  larsim_MCCheater_BackTrackerService_service
#  larsim_MCCheater_ParticleInventoryService_service
#  ${ART_FRAMEWORK_SERVICES_REGISTRY}
#  cetlib cetlib_except
#  ${Boost_FILESYSTEM_LIBRARY}            
#  ${Boost_SYSTEM_LIBRARY}
#  ${ROOT_BASIC_LIB_LIST}
#  DICT_LIBRARIES   lardataobj_RecoBase
  )

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
      ${GLOBIMAGE}
)

install_headers()
install_fhicl()
install_source()

