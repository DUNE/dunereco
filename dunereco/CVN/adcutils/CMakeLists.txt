
if( DEFINED ENV{TENSORFLOW_DIR} )

art_make_library( 
  LIBRARY_NAME     GlobImage
  SOURCE EventImageData.h EventImageData.cxx
  LIBRARIES    CVNFunc
  LIBRARY_NAME_VAR GLOBIMAGE
  )

simple_plugin(SPMultiTpcDump "module"
			larcore_Geometry_Geometry_service
			larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larreco_RecoAlg_ImagePatternAlgs_Tensorflow_PointIdAlg
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
			
			cetlib cetlib_except
			${ROOT_BASIC_LIB_LIST}
      ${GLOBIMAGE}
)
endif()

install_headers()
install_fhicl()
install_source()

