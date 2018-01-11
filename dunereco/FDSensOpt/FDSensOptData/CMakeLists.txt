# this will create both the library and the dictionary
art_make(
  LIB_LIBRARIES
     ${ROOT_BASIC_LIB_LIST} #TLorentzVector
#    lardataobj_RecoBase
#    ${ROOT_EG} # TDatabasePDG
  DICT_LIBRARIES
#     ${ROOT_BASIC_LIB_LIST}
#    nusimdata_SimulationBase
  )
install_source()
install_headers()
