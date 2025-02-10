#include "onnx_graph.h"
#include <algorithm>
onnx::Model::Model(const std::string& model_file_name, 
                   bool& success)
: /*n_inputs(ninputs), n_outputs(noutputs),*/ fSession(nullptr), fEnv(ORT_LOGGING_LEVEL_WARNING, "ONNXModel") {
    // Load the model using ONNX Runtime Inference Session
    // onnxruntime setup
    Ort::SessionOptions session_options;
    fSession = Ort::Session(fEnv, model_file_name.c_str(), session_options);


    // print name/shape of inputs
    Ort::AllocatorWithDefaultOptions allocator;
    
    //TODO make the char * vecs easier

    std::cout << "Input Node Name/Shape (" << fInputNames.size() << "):" << std::endl;
    for (std::size_t i = 0; i < fSession.GetInputCount(); i++) {
        fInputNames.emplace_back(fSession.GetInputNameAllocated(i, allocator).get());
        fInputShapes.emplace_back(fSession.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape());
        std::cout << "\t" << fInputNames.at(i) << " : " << std::endl;
        for (const auto & sv : fInputShapes)
            for (const auto & s : sv)
                std::cout << "\t\t" << s << std::endl;
    }
    // some models might have negative shape values to indicate dynamic shape, e.g., for variable batch size.
    for (auto& sv : fInputShapes) {
        for (auto & s : sv ) {
            if (s < 0) {
                s = 1;
            }
        }
    }

    // print name/shape of outputs
    std::cout << "Output Node Name/Shape (" << fOutputNames.size() << "):" << std::endl;
    for (std::size_t i = 0; i < fSession.GetOutputCount(); i++) {
        fOutputNames.emplace_back(fSession.GetOutputNameAllocated(i, allocator).get());
        auto output_shapes = fSession.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
        std::cout << "\t" << fOutputNames.at(i) << " : " << /*print_shape(output_shapes) <<*/ std::endl;
    }

    // std::transform(std::begin(fInputNames), std::end(fInputNames), std::begin(fInputNames_char),
    //                 [&](const std::string& str) { return str.c_str(); });
    
    // std::transform(std::begin(fOutputNames), std::end(fOutputNames), std::begin(fOutputNames_char),
    //                 [&](const std::string& str) { return str.c_str(); });


    success = true;
}

template <typename T>
Ort::Value onnx::Model::vec_to_tensor(
    std::vector<T>& data, const std::vector<std::int64_t>& shape) {
  Ort::MemoryInfo mem_info =
      Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator,
                                 OrtMemType::OrtMemTypeDefault);
  auto tensor = Ort::Value::CreateTensor<T>(
    mem_info, data.data(), data.size(), shape.data(), shape.size());
  return tensor;
}

onnx::Model::ModelOutput onnx::Model::run(const onnx::Model::ModelInput & x, long long int samples) {
    ModelOutput result;


    std::vector<const char*> fInputNames_char(fInputNames.size(), nullptr);
    std::transform(std::begin(fInputNames), std::end(fInputNames), std::begin(fInputNames_char),
                    [&](const std::string& str) { return str.c_str(); });
    std::vector<const char*> fOutputNames_char(fOutputNames.size(), nullptr);
    std::transform(std::begin(fOutputNames), std::end(fOutputNames), std::begin(fOutputNames_char),
                    [&](const std::string& str) { return str.c_str(); });

    for (auto & cstr : fOutputNames_char) {
        std::cout << cstr << std::endl;
    }

    size_t rows = x.front().size(),
           cols = x.front().front().size(),
           nviews = x.front().front().front().size();
    std::vector<std::vector<float>> view_vals(nviews);
    
    size_t temp_samples = 1; //One sample for now
    for(size_t view = 0; view < nviews; ++view) {
        auto & these_vals = view_vals[view];
        for (size_t s = 0; s < temp_samples; ++s) {
            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    these_vals.push_back(x[s][r][c][view]);
                }
            }
        }
    }
    
    std::vector<Ort::Value> input_tensors;
    for (size_t i = 0; i < nviews; ++i) {
        input_tensors.emplace_back(vec_to_tensor<float>(view_vals[i], fInputShapes[i]));
    }
    auto onnx_output = fSession.Run(
        Ort::RunOptions{nullptr}, fInputNames_char.data(), input_tensors.data(),
        fInputNames_char.size(), fOutputNames_char.data(), fOutputNames_char.size());


    //Assume only one sample -- so make a new output vector the same size of the output
    // then we'll have possibly several output vals per output bucket
    result.emplace_back(
        std::vector<std::vector<float>>()
    );
    for (size_t i = 0; i < onnx_output.size(); ++i) {
        auto & output = onnx_output[i];
        int nelements = 1;
        auto this_shape = output.GetTensorTypeAndShapeInfo().GetShape();
        for (auto & s : this_shape) {
            nelements *= s;
        }
        auto * ptr = output.GetTensorMutableData<float>();
        result.back().push_back(
            std::vector<float>(ptr, ptr+nelements)
        );
    }


    return result;
}