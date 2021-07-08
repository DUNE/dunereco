////////////////////////////////////////////////////////////////////////
/// \file    TFRegNetHandler.h
/// \brief   TFRegNetHandler for RegCNN modified from TFNetHandler.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_TFNETHANDLER_H
#define REGCNN_TFNETHANDLER_H

#include <vector>
#include <memory>

#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/func/RegPixelMap3D.h"
#include "fhiclcpp/ParameterSet.h"
#include "dune/RegCNN/func/RegCNN_TF_Graph.h"
//#include "larreco/RecoAlg/ImagePatternAlgs/Tensorflow/TF/tf_graph.h"

namespace cnn
{

  /// Wrapper for caffe::Net which handles construction and prediction
  class TFRegNetHandler
  {
  public:

    /// Constructor which takes a pset with DeployProto and ModelFile fields
    TFRegNetHandler(const fhicl::ParameterSet& pset);

    /// Return prediction arrays for RegPixelMap
    std::vector<float> Predict(const RegPixelMap& pm);
    std::vector<float> Predict(const RegPixelMap& pm, const std::vector<float> cm_list);

    std::vector<float> PredictNuEEnergy(const RegPixelMap& pm);

  private:

    std::string  fLibPath;  ///< Library path (typically dune_pardata...)
    std::string  fTFProtoBuf;  ///< location of the tf .pb file in the above path
    unsigned int fInputs;   ///< Number of tdcs for the network to classify
    std::vector<std::string> fOutputName;
    std::vector<bool> fReverseViews; ///< Do we need to reverse any views?
    std::unique_ptr<tf::RegCNNGraph> fTFGraph; ///< Tensorflow graph

  };

}

#endif  // REGCNN_TFNETHANDLER_H
