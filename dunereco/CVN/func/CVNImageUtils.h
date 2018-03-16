////////////////////////////////////////////////////////////////////////
/// \file    CVNImageUtils.h
/// \brief   Utilities for producing images for the CVN
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CVN_IMAGE_UTILS_H
#define CVN_IMAGE_UTILS_H

#include <vector>

#include "dune/CVN/func/PixelMap.h"

namespace cvn
{

  /// Class containing some utility functions for all things CVN
  class CVNImageUtils
  {
  public:
    CVNImageUtils();
    CVNImageUtils(unsigned int nWires, unsigned int nTDCs, unsigned int nViews);
    ~CVNImageUtils();

    /// Convert the hit charge into the range 0 to 255 required by the CVN
    unsigned char ConvertChargeToChar(float charge);

    /// Set up the image size that we want to have
    void SetImageSize(unsigned int nWires, unsigned int nTDCs, unsigned int nViews);
   
    /// Function to set any views that need reversing
    void SetViewReversal(bool reverseX, bool reverseY, bool reverseZ);
    void SetViewReversal(std::vector<bool> reverseViews);

    /// Set the log scale for charge
    void SetLogScale(bool setLog);

    /// Set the input pixel map size
    void SetPixelMapSize(unsigned int nWires, unsigned int nTDCs);

    /// Convert a Pixel Map object into a single pixel array with an image size nWire x nTDC
    void ConvertPixelMapToPixelArray(const PixelMap &pm, std::vector<unsigned char> &pix);

    /// Convert (up to) three vectors (sorted in the same way as the vectors in the PixelMap object) 
    /// into a single pixel array with an image size nWire x nTDC
    void ConvertChargeVectorsToPixelArray(std::vector<float> &v0pe, std::vector<float> &v1pe,
                                          std::vector<float> &v2pe, std::vector<unsigned char> &pix);  
  private:

    /// Get the minimum and maximum wires from the pixel map needed to make the image
    void GetMinMaxWires(std::vector<float> &wireCharges, unsigned int &minWire, unsigned int &maxWire); 

    /// Get the minimum and maximum tdcs from the pixel map needed to make the image
    void GetMinMaxTDCs(std::vector<float> &tdcCharges, unsigned int &minTDC, unsigned int &maxTDC); 

    /// Funtion to actually reverse the view
    void ReverseView(std::vector<float> &peVec);

    /// Number of views of each event
    unsigned int fNViews;

    /// Number of wires to use for the image width
    unsigned int fNWires;

    /// Number of TDCs to use for the image height
    unsigned int fNTDCs;

    /// Input pixel map sizes
    unsigned int fPixelMapWires;
    unsigned int fPixelMapTDCs; 

    /// Vector of bools to decide if any views need to be reversed
    std::vector<bool> fViewReverse;

    /// Use a log scale for charge?
    bool fUseLogScale;

  };

}

#endif  // CVN_IMAGE_UTILS_H
