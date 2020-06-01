
if( DEFINED ENV{TENSORFLOW_DIR} )

art_make_library( 
  LIBRARY_NAME     GlobImage
  SOURCE EventImageData.h EventImageData.cxx
  LIBRARIES    dune_CVN_func
  LIBRARY_NAME_VAR GLOBIMAGE
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  )

simple_plugin(SPMultiTpcDump "module"
			larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larrecodnn_ImagePatternAlgs_Tensorflow_PointIdAlg
			larreco_RecoAlg_PMAlg
			nusimdata_SimulationBase
			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
			${ART_ROOT_IO_TFILESERVICE_SERVICE}
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

