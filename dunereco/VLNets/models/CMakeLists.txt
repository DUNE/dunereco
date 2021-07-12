
# Block Copy&Paste from TrackPID/tf/CMakeLists.txt
if (DEFINED ENV{TENSORFLOW_DIR})
    if ("${OSTYPE}" MATCHES "slf6")
        set(FOR_TENSORFLOW "-lrt")
    endif()
endif()

include_directories( $ENV{TENSORFLOW_INC}/third_party/eigen3 )
include_directories( $ENV{TENSORFLOW_INC}/absl )

art_make_library(
    LIBRARY_NAME VLNModels
    SOURCE
        tf_model/ModelConfig.cxx
        tf_model/TFModel.cxx
        zoo/VLNEnergyModel.cxx
    LIBRARIES
        ${TENSORFLOW}
        ${FOR_TENSORFLOW}
)

cet_remove_compiler_flags(CXX -pedantic)

install_headers(SUBDIRS tf_model zoo)
install_source(SUBDIRS tf_model zoo)

