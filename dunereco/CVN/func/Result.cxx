////////////////////////////////////////////////////////////////////////
/// \file    Result.h
/// \brief   Result for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "dune/CVN/func/Result.h"

namespace cvn
{

  Result::Result(const float* output, unsigned int& nOutputs):
  fOutput(nOutputs)
  {
    for(size_t i = 0; i < nOutputs; ++i) fOutput[i] = output[i];
  }

  Result::Result(const std::vector<float> output){
    fOutput = output; 
  }

  Result::Result():
  fOutput()
  {}

  unsigned int Result::ArgMax(){
    // Get the max element iterator and convert to vector index
    return std::max_element(fOutput.begin(),fOutput.end()) - fOutput.begin();
  }

  float Result::Max(){
    // Get the maximum value by dereferencing the iterator
    return *std::max_element(fOutput.begin(),fOutput.end());
  }

  unsigned int Result::NOutput(){
    return fOutput.size();
  }

  /// Return the predicted interaction type
  TFResultType Result::PredictedInteractionType(){
    return static_cast<TFResultType>(this->ArgMax());
  }

  /// Return the numu flavour probability
  float Result::GetNumuProbability(){
    return fOutput[TFResultType::kTFNumuQE] + fOutput[TFResultType::kTFNumuRes]
         + fOutput[TFResultType::kTFNumuDIS] + fOutput[TFResultType::kTFNumuOther];
  }
  
  /// Return the nue flavour probability
  float Result::GetNueProbability(){
    return fOutput[TFResultType::kTFNueQE] + fOutput[TFResultType::kTFNueRes]
         + fOutput[TFResultType::kTFNueDIS] + fOutput[TFResultType::kTFNueOther];
  }

  /// Return the nutau flavour probability
  float Result::GetNutauProbability(){
    return fOutput[TFResultType::kTFNutauQE] + fOutput[TFResultType::kTFNutauRes]
         + fOutput[TFResultType::kTFNutauDIS] + fOutput[TFResultType::kTFNutauOther];
  }

  /// Return the NC probability
  float Result::GetNCProbability(){
    return fOutput[TFResultType::kTFNC];
  }

}
