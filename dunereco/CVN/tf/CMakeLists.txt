include_directories ( $ENV{TENSORFLOW_INC}/absl )

# TODO -- replace when Onnx is ups'd
include_directories ( $ENV{ONNXRUNTIME_ROOTDIR}/include )

cet_add_compiler_flags(CXX -Wno-pedantic)

cet_make_library(
  LIBRARY_NAME
    CVN_tf
  SOURCE
    tf_graph.cc
    tf_bundle.cc
  LIBRARIES
    pthread
    TensorFlow::cc
    TensorFlow::framework
)

link_directories($ENV{ONNXRUNTIME_ROOTDIR}/lib)

cet_make_library(
  LIBRARY_NAME
    ONNXObjs
  SOURCE
    onnx_graph.cxx
  LIBRARIES
    onnxruntime
)

cet_make_exec(
  NAME
    test_onnx
  SOURCE
    test_onnx.cxx
  LIBRARIES
    ONNXObjs
)

install_headers()
install_fhicl()
install_source()
