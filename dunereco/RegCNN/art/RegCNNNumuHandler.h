////////////////////////////////////////////////////////////////////////
/// \file    RegCNNNumuHandler.h
/// \brief   RegCNNNumuHandler for numu energy estimation
/// \author  Wenjie Wu - wenjieww@uci.edu 
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_NUMUHANDLER_H
#define REGCNN_NUMUHANDLER_H

#include <vector>
#include <memory>

#include "fhiclcpp/ParameterSet.h"
#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/art/TFRegNetHandler.h"

namespace cnn
{

  /// Wrapper for caffe::Net which handles construction and prediction
  class RegCNNNumuHandler
  {
    public:

      /// Constructor which takes a pset with DeployProto and ModelFile fields
      RegCNNNumuHandler(const fhicl::ParameterSet& pset);

      /// Return prediction arrays for RegPixelMap
      std::vector<float> Predict(const RegPixelMap& pm, bool fLongestTrackContained);

    private:

      cnn::TFRegNetHandler fTFHandlerContained;
      cnn::TFRegNetHandler fTFHandlerExiting;

    };
}
#endif  // REGCNN_NUMUHANDLER_H
