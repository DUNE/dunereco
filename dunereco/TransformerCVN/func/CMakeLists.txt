#for classes_def.xml
# not ready for c2 yet

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make( BASENAME_ONLY
  LIBRARY_NAME     TransformerCVNFunc
  LIB_LIBRARIES
      nusimdata::SimulationBase
      RegCNNFunc
  lardataobj::RecoBase
  larsim::MCCheater_BackTrackerService_service
  larsim::MCCheater_ParticleInventoryService_service
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except::cetlib_except
  Boost::filesystem
  DICT_LIBRARIES   lardataobj::RecoBase
  TransformerCVNFunc
  )

install_headers()
install_source()
