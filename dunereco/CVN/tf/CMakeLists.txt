include_directories ( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make(BASENAME_ONLY
  LIB_LIBRARIES
  pthread
  TensorFlow::cc
  TensorFlow::framework
  )

install_headers()
install_fhicl()
install_source()

