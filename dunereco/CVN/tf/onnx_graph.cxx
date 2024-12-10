#include "onnx_graph.h"
onnx::Model::Model(const std::string& model_file_name, 
                   bool& success)
: /*n_inputs(ninputs), n_outputs(noutputs),*/ fSession(nullptr) {
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

    success = true;
}

// std::vector< std::vector< std::vector<float> > > onnx::Model::run(
// 	const std::vector<  std::vector<  std::vector< std::vector<float> > > > & x,
// 	long long int samples)
// {
//     if ((samples == 0) || x.empty() || x.front().empty() || x.front().front().empty() || x.front().front().front().empty())
//         return std::vector< std::vector< std::vector<float> > >();

//     if ((samples == -1) || (samples > (long long int)x.size())) { samples = x.size(); }

//     long long int
//               rows = x.front().size(),
//               cols = x.front().front().size(),
//               depth = x.front().front().front().size();

//     //Change this -- might not need it for Onnx
//     std::vector< tensorflow::Tensor > _x;
//     // Single-input network
//     if (n_inputs == 1)
//     {
//         _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, depth })));
//         auto input_map = _x[0].tensor<float, 4>();
//         for (long long int s = 0; s < samples; ++s) {
//             const auto & sample = x[s];
//             for (long long int r = 0; r < rows; ++r) {
//                 const auto & row = sample[r];
//                 for (long long int c = 0; c < cols; ++c) {
//                     const auto & col = row[c];
//                     for (long long int d = 0; d < depth; ++d) {
//                         input_map(s, r, c, d) = col[d];
//                     }
//                 }
//             }
//         }
//     }
//     // Multi-input network
//     else
//     {
//         for(int i=0; i<depth; ++i){
//             _x.push_back(tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, 1 })));
//         }

//         //tensorflow::Tensor _x(tensorflow::DT_FLOAT, tensorflow::TensorShape({ samples, rows, cols, depth }));

//         for(int view=0; view<depth; ++view){
//             auto input_map = _x[view].tensor<float, 4>();
//             for (long long int s = 0; s < samples; ++s) {
//                 const auto & sample = x[s];
//                 for (long long int r = 0; r < rows; ++r) {
//                     const auto & row = sample[r];
//                     for (long long int c = 0; c < cols; ++c) {
//                         const auto & col = row[c];
//                         long long int d = view;
//                         input_map(s, r, c, 0) = col[d];
//                     }
//                 }
//             }
//         }
//     }

//     return run(_x);
// }
