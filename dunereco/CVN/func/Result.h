////////////////////////////////////////////////////////////////////////
/// \file    Result.h
/// \brief   Result for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_RESULT_H
#define CVN_RESULT_H

#include <vector>

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
    unsigned int Max();

    /// Number of outputs, i.e. size of vector
    unsigned int NOutput(){
      return fOutput.size();
    };

    std::vector<float> fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // CVN_RESULT_H

