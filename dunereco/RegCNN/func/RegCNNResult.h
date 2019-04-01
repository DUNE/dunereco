////////////////////////////////////////////////////////////////////////
/// \file    RegCNNResult.h
/// \brief   RegCNNResult for RegCNN modified from Result.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_RESULT_H
#define REGCNN_RESULT_H

#include <vector>

namespace cnn
{
  /// RegCNNResult, basic output of CNN neural net
  class RegCNNResult
  {
  public:
    RegCNNResult(const float* output, unsigned int& nOutputs);
    // Vector version of the constructor
    RegCNNResult(const std::vector<float> output);
    RegCNNResult();

    std::vector<float> fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // REGCNN_RESULT_H

