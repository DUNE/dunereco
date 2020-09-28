////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap3D.cxx
/// \brief   RegPixelMap3D for RegCNN modifed from PixelMap.cxx
/// \author  Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <iomanip>
#include "dune/RegCNN/func/RegPixelMap3D.h"

namespace cnn
{

  RegPixelMap3D::RegPixelMap3D(const RegCNNBoundary3D& bound, bool prongOnly):
      fBound(bound),
      fProngOnly(prongOnly),
      fInPM(0),
      fPE(fBound.NBins(0)*fBound.NBins(1)*fBound.NBins(2))
  {
      x_axis.Set(fBound.NBins(0), fBound.StartPos(0), fBound.StopPos(0));
      y_axis.Set(fBound.NBins(1), fBound.StartPos(1), fBound.StopPos(1));
      z_axis.Set(fBound.NBins(2), fBound.StartPos(2), fBound.StopPos(2));
  }

  void RegPixelMap3D::AddHit(float rel_x, float rel_y, float rel_z, float charge) 
  {
      if (fBound.IsWithin(rel_x, rel_y, rel_z)) {
          fInPM = 1;
          int xbin = x_axis.FindBin(rel_x);
          int ybin = y_axis.FindBin(rel_y);
          int zbin = z_axis.FindBin(rel_z);
          fPE[LocalToIndex(xbin-1, ybin-1, zbin-1)] += charge;
      }
  }

  unsigned int RegPixelMap3D::LocalToIndex(const unsigned int& bin_x, 
          const unsigned int& bin_y, const unsigned int& bin_z) const 
  {
      unsigned int index = bin_x*fBound.NBins(1)*fBound.NBins(2) + bin_y*fBound.NBins(2) + bin_z%fBound.NBins(2);
      assert(index < fPE.size());
      return index;
  }

  std::ostream& operator<<(std::ostream& os, const RegPixelMap3D& m)
  {
    os << "RegPixelMap3D...... ";
    return os;
  }
}
