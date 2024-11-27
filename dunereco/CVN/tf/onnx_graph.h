#include <onnxruntime_cxx_api.h>
#include <iostream>

namespace onnx
{

    class Model {
    public:
        int n_inputs = 1;
        int n_outputs = 1;

        static std::unique_ptr<Model> create(const std::string& model_file_name, const std::vector<std::string>& outputs = {}, int ninputs = 1, int noutputs = 1) {
            bool success;
            std::unique_ptr<Model> ptr(new Model(model_file_name, outputs, success, ninputs, noutputs));
            if (success) { return ptr; }
            else { return nullptr; }
        }

        ~Model() {}

        // Run with a vector of 2D inputs (batch of inputs)
        std::vector<float> run(const std::vector<std::vector<float>>& x) {
            std::vector<float> result;
            // Convert the input vector to a tensor and run the model
            // Use ONNXRuntime API to run
            return result;
        }

        // Process vector of 3D inputs
        std::vector<std::vector<std::vector<float>>> run(const std::vector<std::vector<std::vector<std::vector<float>>>>& x, long long int samples = -1) {
            std::vector<std::vector<std::vector<float>>> result;
            // Handle multi-dimensional inputs and run ONNX inference
            return result;
        }

    private:
        /// Constructor, initializes ONNX model
        Model(const std::string& model_file_name, const std::vector<std::string>& outputs, bool& success, int ninputs, int noutputs);
        // : n_inputs(ninputs), n_outputs(noutputs), fSession(nullptr) {
        //     // Load the model using ONNX Runtime Inference Session
        //     // onnxruntime setup
        //     Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ONNXModel");
        //     Ort::SessionOptions session_options;
        //     fSession = Ort::Session(env, model_file_name.c_str(), session_options);


        //     // print name/shape of inputs
        //     Ort::AllocatorWithDefaultOptions allocator;
            
        //     std::cout << "Input Node Name/Shape (" << fInputNames.size() << "):" << std::endl;
        //     for (std::size_t i = 0; i < fSession.GetInputCount(); i++) {
        //         fInputNames.emplace_back(fSession.GetInputNameAllocated(i, allocator).get());
        //         fInputShapes = fSession.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
        //         std::cout << "\t" << fInputNames.at(i) << " : " << print_shape(fInputShapes) << std::endl;
        //     }
        //     // some models might have negative shape values to indicate dynamic shape, e.g., for variable batch size.
        //     for (auto& s : fInputShapes) {
        //         if (s < 0) {
        //             s = 1;
        //         }
        //     }

        //     // print name/shape of outputs
        //     std::cout << "Output Node Name/Shape (" << fOutputNames.size() << "):" << std::endl;
        //     for (std::size_t i = 0; i < fSession.GetOutputCount(); i++) {
        //         fOutputNames.emplace_back(fSession.GetOutputNameAllocated(i, allocator).get());
        //         auto output_shapes = fSession.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
        //         std::cout << "\t" << fOutputNames.at(i) << " : " << print_shape(output_shapes) << std::endl;
        //     }

        //     // TODO -- CONSIDER CHECKING THIS -- add in input and output shapes
        //     // assert(fInputNames.size() == 1 && fOutputNames.size() == 1);

        // }

        // ONNX Runtime session for inference
        Ort::Session fSession;

        // Input and output names
        std::vector<std::string> fInputNames;
        std::vector<std::string> fOutputNames;
        std::vector<std::int64_t> fInputShapes;
    };
}
