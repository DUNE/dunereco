# this will create both the library and the dictionary
art_make(
  LIB_LIBRARIES
     ROOT::GenVector #LorentzVector
  DICT_LIBRARIES
     lardataobj::RecoBase #just for associations
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types
install_source()
install_headers()
