# this will create both the library and the dictionary
art_make(
  LIB_LIBRARIES
     ${ROOT_GENVECTOR} #LorentzVector
  DICT_LIBRARIES
     lardataobj_RecoBase #just for associations
  )
install_source()
install_headers()
