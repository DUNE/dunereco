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

art_make(BASENAME_ONLY
#  LIBRARY_NAME      CVNArt
  EXCLUDE ${EXCLUDE_TF}
  LIB_LIBRARIES 
  dunereco::CVN_func
  dunereco::CVN_tf
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
  cetlib::cetlib cetlib_except
  ROOT_BASIC_LIB_LIST
  lardataobj::RecoBase
  lardata::Utilities
  nusimdata::SimulationBase
  larsim_MCCheater_ParticleInventoryService_service
  gsl
  hep_hpc_hdf5
  hep_hpc_Utilities
  ${HDF5_LIBRARIES}
  MVAAlg
  Boost::filesystem
  MODULE_LIBRARIES 
  dunereco::CVN_func
  dunereco::CVN_tf
  dunereco::CVN_art
  stdc++fs
  )


install_headers()
install_fhicl()
install_source()
