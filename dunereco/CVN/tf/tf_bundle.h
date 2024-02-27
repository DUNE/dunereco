#ifndef TF_BUNDLE_H
#define TF_BUNDLE_H

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

#include <memory>
#include <vector>
#include <string>

class Bundle{


public:
    static std::unique_ptr<Bundle> create(const char* bundle_file_name, const std::vector<std::string> & outputs = {}, int ninputs = 1, int noutputs = 1){
        bool success;
        std::unique_ptr<Bundle> ptr(new Bundle(bundle_file_name, outputs, success, ninputs, noutputs));
        if (success){
             return ptr;
        }
        else {
             return nullptr;   
        }
    }


    std::vector<std::vector<std::vector<float> > > run(const char* bundle_file_name, std::vector<std::vector<std::vector<std::vector<float> > > >  x);
    tensorflow::Tensor to_tensor(std::vector<std::vector<std::vector<std::vector<float> > > > x); 
    ~Bundle();
private:

    Bundle(const char* bundle_file_name, const std::vector<std::string> & outputs, bool & success, int ninputs, int noutputs);

};
#endif  
