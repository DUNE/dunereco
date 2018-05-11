////////////////////////////////////////////////////////////////////////
/// \file    Result.h
/// \brief   Result for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CVN_RESULT_H
#define CVN_RESULT_H

#include <vector>
#include "dune/CVN/func/InteractionType.h"

namespace cvn
{
  /// Result, basic output of CVN neural net
  class Result
  {
  public:
    Result(const float* output, unsigned int& nOutputs);
    // Vector version of the constructor
    Result(const std::vector<float> output);
    //Result(const float* output, unsigned int& nOutputs, const float* features, unsigned int& nFeatures);
    Result();

    /// Index of maximum value in vector
    unsigned int ArgMax();

    /// Maximum value in vector
    float Max();

    /// Return the predicted interaction type
    TFResultType PredictedInteractionType();

    /// Return the numu flavour probability
    float GetNumuProbability();

    /// Return the nue flavour probability
    float GetNueProbability();

    /// Return the nutau flavour probability
    float GetNutauProbability();

    /// Return the NC probability
    float GetNCProbability();

    /// Number of outputs, i.e. size of vector
    unsigned int NOutput();

    std::vector<float> fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // CVN_RESULT_H

