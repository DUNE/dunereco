#for classes_def.xml
# not ready for c2 yet

include_directories( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make( BASENAME_ONLY
  LIBRARY_NAME     RegCNNFunc
  LIB_LIBRARIES    nusimdata::SimulationBase
  lardataobj_RecoBase
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  ART_FRAMEWORK_SERVICES_REGISTRY
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  Boost::filesystem
  ROOT_BASIC_LIB_LIST
  TensorFlow::framework
  TensorFlow::cc
  DICT_LIBRARIES   lardataobj_RecoBase
  RegCNNFunc
  )

install_headers()
install_source()
