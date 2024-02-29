#for classes_def.xml
# not ready for c2 yet

include_directories( $ENV{TENSORFLOW_INC}/absl ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

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
  DICT_LIBRARIES   lardataobj::RecoBase
  RegCNNFunc
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_source()
