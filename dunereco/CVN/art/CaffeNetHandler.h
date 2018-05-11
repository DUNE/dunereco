////////////////////////////////////////////////////////////////////////
/// \file    CaffeNetHandler.h
/// \brief   CaffeNetHandler for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_CAFFENETHANDLER_H
#define CVN_CAFFENETHANDLER_H


#include <array>
#include <vector>
#include <map>

#define CPU_ONLY
// Suppress warnings originating in Caffe that we can't do anything about
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#include "caffe/caffe.hpp"
#include "caffe/layers/memory_data_layer.hpp"
#pragma GCC diagnostic pop

#include <glog/logging.h>

#include "dune/CVN/func/PixelMap.h"
#include "fhiclcpp/ParameterSet.h"


namespace cvn
{
  /// Wrapper for caffe::Net which handles construction and prediction
  class CaffeNetHandler
  {
  public:
    /// Basic constructor, takes path to deploy prototxt and model weights
//    CaffeNetHandler(const std::string& libPath,
//                    const std::string& deployProto,
//                    const std::string& modelFile,
//                    const std::string& featureMap = "");

    /// Constructor which takes a pset with DeployProto and ModelFile fields
    CaffeNetHandler(const fhicl::ParameterSet& pset);

    /// Number of outputs in neural net
    int NOutput() const;

    /// Number of outputs in neural net
    int NFeatures() const;

    /// Return prediction array for PixelMap
    //const float* Predict(const PixelMap& pm);

    /// Return prediction arrays for PixelMap
    std::pair<const float*, const float*> Predict(const PixelMap& pm);

    std::string  fLibPath;  ///< location of deploy prototxt
    std::string  fDeployProto;  ///< location of deploy prototxt
    std::string  fModelFile;    ///< location of model weights
    std::string  fFeatureMap;    ///< which feature map to extract if any
    int          fLogLevel;     ///< Stashed version of the log level in glog
    bool         fUseLogChargeScale;  ///< Is the charge using a log scale?
    unsigned int fImageWires;  ///< Number of wires for the network to classify
    unsigned int fImageTDCs;   ///< Number of tdcs for the network to classify
    std::vector<bool> fReverseViews; ///< Do we need to reverse any views?
    caffe::Net<float> fNet;    ///< Network object

    caffe::Datum PixelMapToDatum(const PixelMap& pm);

  };


  /// This function just returns the thing you give it, but in the meanwhile
  /// it sets the glog level to what you pass in.
  /// It's a hokey thing to do, but I wanted to call it in the initializer
  /// list for the class before the fNet was constructed.
  int LogLevel(int level);

}

#endif  // CVN_CAFFENETHANDLER_H
