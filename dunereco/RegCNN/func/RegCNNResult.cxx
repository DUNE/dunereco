////////////////////////////////////////////////////////////////////////
/// \file    RegCNNResult.h
/// \brief   RegCNNResult for RegCNN modified from Result.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "dune/RegCNN/func/RegCNNResult.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cnn
{

  RegCNNResult::RegCNNResult(const float* output, unsigned int& nOutputs):
  fOutput(nOutputs)
  {
    for(size_t i = 0; i < nOutputs; ++i) fOutput[i] = output[i];
  }

  RegCNNResult::RegCNNResult(const std::vector<float> output){
    fOutput = output; 
  }

  RegCNNResult::RegCNNResult():
  fOutput()
  {}

}
