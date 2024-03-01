
include_directories( $ENV{TENSORFLOW_INC}/absl )

cet_add_compiler_flags(CXX -Wno-pedantic)

art_make_library(
    LIBRARY_NAME VLNModels
    SOURCE
        tf_model/ModelConfig.cxx
        tf_model/TFModel.cxx
        zoo/VLNEnergyModel.cxx
    LIBRARIES
        TensorFlow::cc
        TensorFlow::framework
)

install_headers(SUBDIRS tf_model zoo)
install_source(SUBDIRS tf_model zoo)

