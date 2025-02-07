////////////////////////////////////////////////////////////////////////
/// \file    TransformerCVNResult.h
/// \brief   TransformerCVNResult for TransformerCVN modified from Result.h
/// \author  Alejandro Yankelevich - ayankele@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "dunereco/TransformerCVN/func/TransformerCVNResult.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cnn
{

  TransformerCVNResult::TransformerCVNResult(const float* output, unsigned int& nOutputs):
  fOutput(nOutputs)
  {
    for(size_t i = 0; i < nOutputs; ++i) fOutput[i] = output[i];
  }

  TransformerCVNResult::TransformerCVNResult(const std::vector<float> output){
    fOutput = output; 
  }

  TransformerCVNResult::TransformerCVNResult():
  fOutput()
  {}

}
