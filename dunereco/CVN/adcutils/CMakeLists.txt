
if( DEFINED ENV{TENSORFLOW_DIR} )

art_make_library( 
  LIBRARY_NAME     GlobImage
  LIBRARY_NAME_VAR GLOBIMAGE
  SOURCE EventImageData.h EventImageData.cxx
  LIBRARIES    dunereco_CVN_func
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
			ROOT::Core
 			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except
      GLOBIMAGE
)
endif()

install_headers()
install_fhicl()
install_source()

