////////////////////////////////////////////////////////////////////////////////////////////////////
//// Class:       Bundle
//// Authors:     A. Higuera, 			from DUNE, Rice U, Feb, 2024
//// 
//// Iterface to run Tensorflow saved model bundle to a file.
////
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include "tf_bundle.h"
#include "tensorflow/core/public/session_options.h"
#include "tensorflow/cc/saved_model/loader.h"
#include "tensorflow/cc/saved_model/tag_constants.h"


Bundle::Bundle(const char* bundle_file_name, const std::vector<std::string> & outputs, bool & success, int ninputs, int noutputs)
{

    success = false;
    tensorflow::SavedModelBundle bundle;
    tensorflow::SessionOptions session_options;
    tensorflow::RunOptions run_options;
    tensorflow::Status status = tensorflow::LoadSavedModel(session_options, run_options, bundle_file_name, {tensorflow::kSavedModelTagServe}, &bundle);
    if(!status.ok()) {
        std::cout << "Failed to load saved model: " << status.ToString() << std::endl;
    }
    success = true;

}
tensorflow::Tensor Bundle::to_tensor( std::vector< std::vector< std::vector< std::vector<float> > > > x){

    long long int
              rows = x.front().size(),
              cols = x.front().front().size(),
              depth = x.front().front().front().size();

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
    return input_tensor;
}
std::vector<std::vector<std::vector<float> > >  Bundle::run(const char* bundle_file_name, std::vector<std::vector<std::vector<std::vector<float> > > >  x)
{
    tensorflow::SavedModelBundle bundle;
    tensorflow::SessionOptions session_options;
    tensorflow::RunOptions run_options;
    tensorflow::Status status = tensorflow::LoadSavedModel(session_options, run_options, bundle_file_name, {tensorflow::kSavedModelTagServe}, &bundle);
    std::vector< std::string > fInputNames;
    auto sig_map = bundle.GetSignatures();
    if (!sig_map.contains("serving_default")) {
	std::cout << "Could not find serving_default in model signatures." << std::endl;
    }
 
    auto model_def = sig_map.at("serving_default");
    std::vector<std::string> input_names;
    std::vector< std::string > output_names;  
    for (auto const &p : model_def.inputs()){
        input_names.push_back(p.second.name().c_str());
    }
    for (auto const &p : model_def.outputs()){
        output_names.push_back(p.second.name().c_str());
    }

    tensorflow::Tensor tensor = Bundle::to_tensor(x); 
    std::vector< std::pair<std::string, tensorflow::Tensor> > inputs;

    //To do, what if input > 1?
    if (input_names.size() == 1 ){
       inputs.push_back({input_names[0], tensor});
    }

    std::vector<tensorflow::Tensor> outputs;
    tensorflow::Status tf_status;
 
    tf_status = bundle.session->Run({inputs}, output_names, {}, &outputs);
    if(!tf_status.ok()) {
        std::cout << "Failed to run model: " << status.ToString() << std::endl;
    }
    //NOTE:
    //For some reason, TensorFlow decided to change the order of the outputs
    //Therefore, we need to rearrange them accordingly
    tensorflow::Tensor flavour;
    tensorflow::Tensor protons;
    tensorflow::Tensor pions;
    for (size_t i = 0; i < outputs.size(); ++i) {
        if (outputs[i].shape().dims() == 2 && 
            outputs[i].shape().dim_size(0) == 1 && 
            outputs[i].shape().dim_size(1) == 3) {
            flavour = outputs[i];
        } 
        else if (outputs[i].shape().dims() == 2 && 
            outputs[i].shape().dim_size(0) == 1 && 
            outputs[i].shape().dim_size(1) == 4) {
            protons = outputs[i];
        } 
        else if (outputs[i].shape().dims() == 2 && 
            outputs[i].shape().dim_size(0) == 1 && 
            outputs[i].shape().dim_size(1) == 2) { 
            pions = outputs[i];
        } 
        else{
            std::cout << "Output with wrong shape!" << std::endl;
        }
    }
    
    std::vector<tensorflow::Tensor> final_outputs;
    final_outputs.push_back(flavour); 
    final_outputs.push_back(protons); 
    final_outputs.push_back(pions); 

    std::vector< std::vector< std::vector< float > > > result;
    result.resize(1, std::vector< std::vector< float > >(final_outputs.size()));
    for (size_t s = 0; s < 1; ++s){
            for (size_t o = 0; o < final_outputs.size(); ++o){
                size_t n = final_outputs[o].dim_size(1);                
                auto output_map = final_outputs[o].tensor<float, 2>();
                result[s][o].resize(final_outputs[o].dim_size(1));
                std::vector< float > & vs = result[s][o];
                for (size_t i = 0; i < n; ++i){
                    vs[i] = output_map(s, i);
                }            
            }
        }

    return result;
}
Bundle::~Bundle()
{

}


