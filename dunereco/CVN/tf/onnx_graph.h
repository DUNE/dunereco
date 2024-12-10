#include <onnxruntime_cxx_api.h>
#include <iostream>

namespace onnx
{

    class Model {
    public:
        // int n_inputs = 1;
        // int n_outputs = 1;

        static std::unique_ptr<Model> create(const std::string& model_file_name) {
            bool success;
            std::unique_ptr<Model> ptr(new Model(model_file_name, success));
            if (success) { return ptr; }
            else { return nullptr; }
        }

        ~Model() {}

        // // Run with a vector of 2D inputs (batch of inputs)
        // std::vector<float> run(const std::vector<std::vector<float>>& x) {
        //     std::vector<float> result;
        //     // Convert the input vector to a tensor and run the model
        //     // Use ONNXRuntime API to run
        //     return result;
        // }

        // Process vector of 3D inputs
        std::vector<std::vector<std::vector<float>>> run(const std::vector<std::vector<std::vector<std::vector<float>>>>& x, long long int samples = -1) {
            std::vector<std::vector<std::vector<float>>> result;
            // Handle multi-dimensional inputs and run ONNX inference
            return result;
        }

    private:
        /// Constructor, initializes ONNX model
        Model(const std::string& model_file_name, bool& success);
        // ONNX Runtime session for inference
        Ort::Session fSession;

        // Input and output names
        std::vector<std::string> fInputNames;
        std::vector<std::string> fOutputNames;
        std::vector<std::int64_t> fInputShapes;
    };
}

