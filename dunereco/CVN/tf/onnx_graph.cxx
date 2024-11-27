#include "onnx_graph.h"
onnx::Model::Model(const std::string& model_file_name, 
      const std::vector<std::string>& outputs,
      bool& success, int ninputs, int noutputs)
: n_inputs(ninputs), n_outputs(noutputs), fSession(nullptr) {
    // Load the model using ONNX Runtime Inference Session
    // onnxruntime setup
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ONNXModel");
    Ort::SessionOptions session_options;
    fSession = Ort::Session(env, model_file_name.c_str(), session_options);


    // print name/shape of inputs
    Ort::AllocatorWithDefaultOptions allocator;
    
    std::cout << "Input Node Name/Shape (" << fInputNames.size() << "):" << std::endl;
    for (std::size_t i = 0; i < fSession.GetInputCount(); i++) {
        fInputNames.emplace_back(fSession.GetInputNameAllocated(i, allocator).get());
        fInputShapes = fSession.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
        std::cout << "\t" << fInputNames.at(i) << " : " /*<< print_shape(fInputShapes)*/ << std::endl;
    }
    // some models might have negative shape values to indicate dynamic shape, e.g., for variable batch size.
    for (auto& s : fInputShapes) {
        if (s < 0) {
            s = 1;
        }
    }

    // print name/shape of outputs
    std::cout << "Output Node Name/Shape (" << fOutputNames.size() << "):" << std::endl;
    for (std::size_t i = 0; i < fSession.GetOutputCount(); i++) {
        fOutputNames.emplace_back(fSession.GetOutputNameAllocated(i, allocator).get());
        auto output_shapes = fSession.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
        std::cout << "\t" << fOutputNames.at(i) << " : " << /*print_shape(output_shapes) <<*/ std::endl;
    }

    // TODO -- CONSIDER CHECKING THIS -- add in input and output shapes
    // assert(fInputNames.size() == 1 && fOutputNames.size() == 1);

}
