
include_directories ( $ENV{TENSORFLOW_INC}/absl ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make(BASENAME_ONLY
  LIB_LIBRARIES
  pthread
  TensorFlow::cc
  TensorFlow::framework
  ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()

