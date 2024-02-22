
include_directories( $ENV{TENSORFLOW_INC}/absl ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible

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

