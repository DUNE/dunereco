////////////////////////////////////////////////////////////////////////
/// \file    RegCNNImageUtils.h
/// \brief   Utilities for producing images for the RegCNN
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef RegCNN_IMAGE_UTILS_H
#define RegCNN_IMAGE_UTILS_H

#include <vector>

#include "dune/RegCNN/func/RegPixelMap.h"

namespace cnn
{

  /// Useful typedefs
  typedef std::vector<std::vector<float> > ViewVectorF;
  typedef std::vector<ViewVectorF> ImageVectorF;

  /// Class containing some utility functions for all things RegCNN
  class RegCNNImageUtils
  {
  public:
    RegCNNImageUtils();
    RegCNNImageUtils(unsigned int nWires, unsigned int nTDCs, unsigned int nViews);
    ~RegCNNImageUtils();


    /// Function to set any views that need reversing
    void SetViewReversal(bool reverseX, bool reverseY, bool reverseZ);
    void SetViewReversal(std::vector<bool> reverseViews);

    /// Set the input pixel map size
    void SetPixelMapSize(unsigned int nWires, unsigned int nTDCs);


    // Scale Charges
    float ConvertToScaledCharge(float charge);

    /// Convert a pixel map into an image vector (float version)
    void ConvertPixelMapToImageVectorF(const RegPixelMap &pm, ImageVectorF &imageVec);

    /// Float version of conversion for convenience of TF interface
    void ConvertChargeVectorsToImageVectorF(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                           std::vector<float> &v2pe, ImageVectorF &imageVec);  


  private:

    /// Base function for conversion of the Pixel Map to our required output format
    void ConvertChargeVectorsToViewVectors(std::vector<float> &v0pe, std::vector<float> &v1pe, std::vector<float> &v2pe,
                                  ViewVectorF& view0, ViewVectorF& view1, ViewVectorF& view2);

    /// Make the image vector from the view vectors
    ImageVectorF BuildImageVectorF(ViewVectorF v0, ViewVectorF v1, ViewVectorF v2);


    /// Funtion to actually reverse the view
    void ReverseView(std::vector<float> &peVec);

    /// Number of views of each event
    unsigned int fNViews;

    /// Input pixel map sizes
    unsigned int fPixelMapWires;
    unsigned int fPixelMapTDCs; 

    /// Vector of bools to decide if any views need to be reversed
    std::vector<bool> fViewReverse;

  };

}

#endif  // RegCNN_IMAGE_UTILS_H
