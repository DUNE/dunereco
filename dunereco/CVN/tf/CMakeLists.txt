include_directories ( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make(BASENAME_ONLY
  LIB_LIBRARIES
  pthread
  PROTOBUF
  TENSORFLOW_CC
  TENSORFLOW_FRAMEWORK
  )

install_headers()
install_fhicl()
install_source()

