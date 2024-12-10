# Tensorflow stuff
if (DEFINED ENV{TENSORFLOW_DIR})
# For some reason we apparently need to care if larreco is checked out.
#if (larreco_not_in_ups)
#set (LARRECO_RECOALG_IMAGEPATTERNALGS_TF larrecodnn_ImagePatternAlgs_Tensorflow_TF)
#else (larreco_not_in_ups)
#cet_find_library( LARRECO_RECOALG_IMAGEPATTERNALGS_TF NAMES larrecodnn_ImagePatternAlgs_Tensorflow_TF PATHS ENV LARRECO_LIB NO_DEFAULT_PATH)
#endif (larreco_not_in_ups)
else(DEFINED ENV{TENSORFLOW_DIR})
set (EXCLUDE_TF TFNetHandler.cxx)
endif (DEFINED ENV{TENSORFLOW_DIR})


include_directories(${HEP_HPC_INCLUDE_DIRS})


# TODO -- replace when Onnx is ups'd
include_directories ( $ENV{ONNXRUNTIME_ROOTDIR}/include )
link_directories($ENV{ONNXRUNTIME_ROOTDIR}/lib)

set(cvn_art_sources 
  PixelMapWireProducer.cxx
  PixelMapSimProducer.cxx
  PixelMapProducer.cxx
)
if (DEFINED ENV{TENSORFLOW_DIR})
  list(APPEND cvn_art_sources TFNetHandler.cxx)
endif()

cet_make_library(
  LIBRARY_NAME
    CVN_art
  SOURCE
    ${cvn_art_sources}
  LIBRARIES
    dunereco::CVN_func
    dunereco::CVN_tf
    dunereco::ONNXObjs
    art::Framework_Core
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::tfile_support
    art_root_io::TFileService_service
    ROOT::Core
    art::Persistency_Common
    art::Persistency_Provenance
    art::Utilities 
    canvas::canvas
    messagefacility::MF_MessageLogger
    z
    fhiclcpp::fhiclcpp
    cetlib::cetlib cetlib_except::cetlib_except
    ROOT::Core
    lardataobj::RecoBase
    lardata::Utilities
    nusimdata::SimulationBase
    larsim::MCCheater_ParticleInventoryService_service
    GSL::gsl
    hep_hpc_hdf5
    hep_hpc_Utilities
    HDF5::HDF5
    MVAAlg
    Boost::filesystem
)

set(cvn_libraries dunereco::CVN_func dunereco::CVN_tf dunereco::CVN_art)

set(cvn_art_producers 
  CVNEvaluator CVNMapper CVNMapperProtoDUNE CVNMapperWire CVNSparseMapper3D
  CVNSparseMapper GCNGraphMaker2D GCNGraphMaker GCNGraphMakerProtoDUNE
  GCNHitGraphMaker
)
foreach(module IN LISTS cvn_art_producers)
  cet_build_plugin(
    ${module} art::EDProducer
    LIBRARIES PRIVATE
      ${cvn_libraries}
  )
endforeach()

set(cvn_art_anaylzers
  CVNEventDump CVNEventDumpProtoDUNE CVNSparseROOT CVNValidation CVNZlibMaker
  CVNZlibMakerProtoDUNE GCNH5 GCNZlibMaker GCNZlibMakerProtoDUNE
)
foreach(module IN LISTS cvn_art_anaylzers)
  cet_build_plugin(
    ${module} art::EDAnalyzer
    LIBRARIES PRIVATE
      ${cvn_libraries}
  )
endforeach()


install_headers()
install_fhicl()
install_source()
