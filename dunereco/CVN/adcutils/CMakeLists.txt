
if( DEFINED ENV{TENSORFLOW_DIR} )

art_make_library( 
  LIBRARY_NAME     GlobImage
  LIBRARY_NAME_VAR GLOBIMAGE
  SOURCE EventImageData.h EventImageData.cxx
  LIBRARIES    dunereco_CVN_func
               larrecodnn::ImagePatternAlgs_Tensorflow_PointIdAlg
                        lardataalg::DetectorInfo
                        lardata::headers
                        art_root_io::tfile_support
  )

cet_build_plugin(SPMultiTpcDump art::module LIBRARIES
			larcore::Geometry_Geometry_service
                        lardataalg::DetectorInfo
			larcorealg::Geometry
			lardataobj::RecoBase
                        lardata::headers
			larreco::Calorimetry
			larrecodnn::ImagePatternAlgs_Tensorflow_PointIdAlg
			larreco::RecoAlg_PMAlg
			nusimdata::SimulationBase
			ROOT::Core ROOT::Tree
 			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities
			canvas::canvas
			messagefacility::MF_MessageLogger
			cetlib::cetlib cetlib_except::cetlib_except
      GLOBIMAGE
)
endif()

install_headers()
install_fhicl()
install_source()

