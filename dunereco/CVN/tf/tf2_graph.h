#ifndef TF2_GRAPH_H
#define TF2_GRAPH_H

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

#include <memory>
#include <vector>
#include <string>

void printHelloTF();
//std::vector< std::vector < std::vector< float > > > run(const std::vector< tensorflow::Tensor > & x);
void my_prediction(std::vector<std::vector<std::vector<std::vector<float> > > > & x);
tensorflow::Tensor to_tensor(std::vector<std::vector<std::vector<std::vector<float> > > > x);
//
#endif  
