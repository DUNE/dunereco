#ifndef TRANSFORMERCVN_RESULT_H
#define TRANSFORMERCVN_RESULT_H

#include <vector>

namespace cnn
{
  /// TransformerCVNResult, basic output of CNN neural net
  class TransformerCVNResult
  {
  public:
    TransformerCVNResult(const float* output, unsigned int& nOutputs);
    // Vector version of the constructor
    TransformerCVNResult(const std::vector<float> output);
    TransformerCVNResult();

    std::vector<float> fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // TRANSFORMERCVN_RESULT_H

