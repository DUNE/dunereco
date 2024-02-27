#include "tf2_graph.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"

void printHelloTF() 
{

    std::cout << "Hello TensorFlow!" << std::endl;
}

tensorflow::Tensor to_tensor( std::vector< std::vector< std::vector< std::vector<float> > > > x){
    long long int
              rows = x.front().size(),
              cols = x.front().front().size(),
              depth = x.front().front().front().size();
    std::cout<<"reading map.."<<std::endl;
    std::cout<<"rows "<<rows<<std::endl;
    std::cout<<"cols "<<cols<<std::endl;
    std::cout<<"depth "<<depth<<std::endl;
	
    auto input_tensor = tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({ 1, rows, cols, depth })); 
    auto input_map = input_tensor.tensor<float, 4>();
    for (long long int s = 0; s < 1; ++s) {
        const auto & sample = x[s];
        for (long long int r = 0; r < rows; ++r) {
            const auto & row = sample[r];
            for (long long int c = 0; c < cols; ++c) {
                const auto & col = row[c];
                for (long long int d = 0; d < depth; ++d) {
                     input_map(s, r, c, d) = col[d];
                }
            }
        }
    }
    std::cout<<"done with to_tensor"<<std::endl;
    return input_tensor;
}
//std::vector< std::vector< std::vector< float > > > run(const std::vector< tensorflow::Tensor > & x)
void my_prediction(std::vector<std::vector<std::vector<std::vector<float> > > > & x)
{


    std::cout << "loading and predicting..."<<std::endl;
    const std::string model_dir = "/exp/dune/data/users/higuera/CNN/atmo_model_tf_2_11/ResNet_20240204/";
    tensorflow::SavedModelBundle bundle;
    tensorflow::SessionOptions session_options;
    tensorflow::RunOptions run_options;
    tensorflow::Status status = tensorflow::LoadSavedModel(session_options, run_options, model_dir, {tensorflow::kSavedModelTagServe}, &bundle);
    if(!status.ok()) {
        std::cout << "Failed to load saved model: " << status.ToString() << std::endl;
    }

    std::cout<<"here"<<std::endl;
   
    tensorflow::Tensor input = to_tensor(x);
    /*
    for(int b = 0; b < 1; ++b) {
    for(int h = 0; h < 200; ++h) {
        for(int w = 0; w < 200; ++w) {
            for(int c = 0; c < 3; ++c) {
                // Accessing the tensor value at [b, h, w, c]
                //                 // Tensor access in TensorFlow C++ API requires flat indexing
                     int index = b * 200 * 200 * 3 + h * 200 * 3 + w * 3 + c;
                     float value = input.flat<float>()(index);
			std::cout<<"value "<<value<<std::endl;
                    }
              }
         }
     }
    */

    //tensor is fine
    //
    //check model
    auto sig_map = bundle.GetSignatures();
    if (!sig_map.contains("serving_default")) {
	std::cout << "Could not find serving_default in model signatures." << std::endl;
    }
 
    auto model_def = sig_map.at("serving_default");
    std::cout<<"ww"<<std::endl;
    for (auto const &p : model_def.inputs())
    {
        std::cout << "key: " << p.first.c_str() << " " << p.second.name().c_str() << std::endl;
    }
    for (auto const &p : model_def.outputs())
    {
        std::cout << "key: " << p.first.c_str() << " " << p.second.name().c_str() << std::endl;
    }

    auto input_name = model_def.inputs().at("input_1").name();

   
    std::cout<<"hey 111"<<std::endl;
    //The above works
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::Status tf_status;
    tf_status = bundle.session->Run({{input_name, input}},
    	                              {"StatefulPartitionedCall:0","StatefulPartitionedCall:1","StatefulPartitionedCall:2"}, // Output operation names
    	                               {}, &outputs);
    std::cout<<"hey 222"<<std::endl;
    if(!tf_status.ok()) {
        std::cout << "Failed to run model: " << status.ToString() << std::endl;
    }

    
    std::cout << "out size " << outputs.size() << std::endl;
    const tensorflow::Tensor& flavour_output = outputs[0]; // Access the 'flavour' output
    const tensorflow::Tensor& protons_output = outputs[1]; // Access the 'protons' output
    const tensorflow::Tensor& pions_output = outputs[2];   // Access the 'pions' output

      
    // TODO: Now process these tensors as needed
    //     // Each line is a placeholder for our eventual inference. 
    std::cout << "Flavour output shape: " << flavour_output.shape().DebugString() << std::endl;
    std::cout << "Protons output shape: " << protons_output.shape().DebugString() << std::endl;
    std::cout << "Pions output shape: " << pions_output.shape().DebugString() << std::endl;
    //std::vector<std::vector<std::vector<float> > > test;
    //return test;
   
}


