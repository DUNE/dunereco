#for classes_def.xml
# not ready for c2 yet

include_directories( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make( BASENAME_ONLY
  LIBRARY_NAME     RegCNNFunc
  LIB_LIBRARIES    nusimdata::SimulationBase
  lardataobj::RecoBase
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except::cetlib_except
  Boost::filesystem
  TensorFlow::framework
  TensorFlow::cc
  ROOT::Hist
  DICT_LIBRARIES   lardataobj::RecoBase
  RegCNNFunc
  )

install_headers()
install_source()
