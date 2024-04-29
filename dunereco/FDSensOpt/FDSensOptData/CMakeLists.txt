# this will create both the library and the dictionary
art_make(
  LIB_LIBRARIES
     ROOT::GenVector #LorentzVector
  DICT_LIBRARIES
     lardataobj::RecoBase #just for associations
  )
install_source()
install_headers()
