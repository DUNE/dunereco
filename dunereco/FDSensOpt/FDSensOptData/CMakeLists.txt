# this will create both the library and the dictionary
art_make(
  LIB_LIBRARIES
     ${ROOT_GENVECTOR} #LorentzVector
  DICT_LIBRARIES
  )
install_source()
install_headers()
