
if( DEFINED ENV{TENSORFLOW_DIR} )

art_make_library( 
  LIBRARY_NAME     GlobImage
  LIBRARY_NAME_VAR GLOBIMAGE
  SOURCE EventImageData.h EventImageData.cxx
  LIBRARIES    dunereco_CVN_func
               ART_ROOT_IO_TFILESERVICE_SERVICE
               ART_ROOT_IO_TFILE_SUPPORT
  )

cet_build_plugin(SPMultiTpcDump art::module
			larcore_Geometry_Geometry_service
                        lardataalg_DetectorInfo
			larcorealg_Geometry
			lardataobj_RecoBase
			larreco_Calorimetry
			larrecodnn_ImagePatternAlgs_Tensorflow_PointIdAlg
			larreco_RecoAlg_PMAlg
			nusimdata::SimulationBase
			ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
 			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			
			cetlib::cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
      GLOBIMAGE
)
endif()

install_headers()
install_fhicl()
install_source()

