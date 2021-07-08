////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap3D.h
/// \brief   RegPixelMap3D for RegCNN modified from PixelMap.h
/// \author  Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_REGPIXELMAP3D_H
#define REGCNN_REGPIXELMAP3D_H

#include <ostream>
#include <vector>

#include "dune/RegCNN/func/RegCNNBoundary3D.h"
#include "dune/RegCNN/func/HitType.h"
#include "TH3F.h"
#include "TAxis.h"

namespace cnn
{

  /// RegPixelMap3D, input to 3D CNN neural net
  class RegPixelMap3D
  {
  public:
    RegPixelMap3D(const RegCNNBoundary3D& bound, const bool& Cropped, const bool& prongOnly);
    RegPixelMap3D() {};
    
    void AddHit(float rel_x, float rel_y, float rel_z, float charge, int hit_prong_tag);
    bool IsCroppedPM() const {return fCropped;};
    std::vector<float> GetPM() const {return fPE;};
    std::vector<float> GetCroppedPM() const {return fPECropped;};

    // Add Finish method in order to determine whether to produce prong only/cropped
    // pixel maps or the full event/uncropped pixel maps
    void Finish();

    unsigned int LocalToIndex(const unsigned int& bin_x,
                              const unsigned int& bin_y,
                              const unsigned int& bin_z)
       const; 
    TH3F* ToTH3() const;
    TH3F* ToCroppedTH3() const;

    RegCNNBoundary3D fBound;
    bool fCropped;
    bool fProngOnly;     //< whether to use prong only pixel map
    unsigned int fInPM;
    std::vector<float> fPE; //< charges of all pixels
    std::vector<float> fPECropped;
    std::vector<int> fProngTag;

  private:
    TAxis x_axis;
    TAxis y_axis;
    TAxis z_axis;

  }; // class RegPixelMap3D

  std::ostream& operator<<(std::ostream& os, const RegPixelMap3D& m);

} // namespace cnn

#endif  // REGCNN_PIXELMAP3D_H
