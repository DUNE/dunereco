
include_directories ( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make(BASENAME_ONLY
  LIB_LIBRARIES
  pthread
  protobuf
  tensorflow_cc
  tensorflow_framework
  )

install_headers()
install_fhicl()
install_source()

