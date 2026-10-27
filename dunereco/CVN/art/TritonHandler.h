////////////////////////////////////////////////////////////////////////
/// \file    TritonHandler.h
/// \brief   TritonHandler for CVN
////////////////////////////////////////////////////////////////////////
#ifndef CVN_TRITONHANDLER_H
#define CVN_TRITONHANDLER_H

#include <string>
#include <vector>
#include <memory>

#include "fhiclcpp/ParameterSet.h"

#include "dunereco/CVN/func/PixelMap.h"
#include "dunereco/CVN/func/InteractionType.h"
#include "dunereco/CVN/func/CVNImageUtils.h"

#include "larrecodnn/ImagePatternAlgs/NuSonic/Triton/TritonClient.h"

namespace cvn
{
  /// Wrapper for a Triton inference client which handles construction and prediction
  class TritonHandler
  {
  public:
    /// Constructor which takes a pset with TritonConfig and image/view fields
    TritonHandler(const fhicl::ParameterSet& pset);

    /// Number of outputs in neural net (i.e. number of output heads)
    int NOutput() const;

    /// Number of input views/features fed to the neural net
    int NFeatures() const;

    /// Return prediction arrays for PixelMap, one vector<float> per output head:
    /// [0] output_is_antineutrino (size 1)
    /// [1] output_flavour         (size 4)
    /// [2] output_interaction     (size 4)
    /// [3] output_protons         (size 4)
    /// [4] output_pions           (size 4)
    /// [5] output_pizeros         (size 4)
    /// [6] output_neutrons        (size 4)
    std::vector<std::vector<float>> Predict(const PixelMap& pm);

    /// Return four element vector with summed numu, nue, nutau and NC elements
    std::vector<float> PredictFlavour(const PixelMap& pm);

  private:
    /// Split a flattened multi-view image into per-view buffers for Triton inputs
    void SplitViews(const ImageVectorF& fullImage,
                     std::vector<std::vector<float>>& viewData) const;

    std::unique_ptr<lartriton::TritonClient> fTritonClient; ///< Triton inference client

    bool         fUseLogChargeScale; ///< Is the charge using a log scale?
    unsigned int fImageWires;        ///< Number of wires for the network to classify
    unsigned int fImageTDCs;         ///< Number of tdcs for the network to classify
    unsigned int fNViews;            ///< Number of views (input tensors)
    std::vector<bool> fReverseViews; ///< Do we need to reverse any views?

    std::vector<std::string> fInputNames;  ///< Names of the Triton input tensors (per view)
    std::vector<std::string> fOutputNames; ///< Names of the Triton output tensors (per head)
  };
}
#endif  // CVN_TRITONHANDLER_H
