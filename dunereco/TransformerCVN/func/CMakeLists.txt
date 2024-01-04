#for classes_def.xml
# not ready for c2 yet

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make( BASENAME_ONLY
  LIBRARY_NAME     TransformerCVNFunc
  LIB_LIBRARIES
      nusimdata::SimulationBase
      RegCNNFunc
  lardataobj_RecoBase
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  ART_FRAMEWORK_SERVICES_REGISTRY
  messagefacility::MF_MessageLogger
  cetlib::cetlib cetlib_except
  Boost::filesystem
  ROOT_BASIC_LIB_LIST
  DICT_LIBRARIES   lardataobj_RecoBase
  TransformerCVNFunc
  )

install_headers()
install_source()
