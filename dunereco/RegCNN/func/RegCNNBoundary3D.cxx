////////////////////////////////////////////////////////////////////////
/// \file    RegCNNBoundary3D.cxx
/// \brief   RegCNNBoundary3D for 3D RegCNN PixelMap modified from CVNBoundary.cxx
/// \author  Wenjie Wu   - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <utility>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dune/RegCNN/func/RegCNNBoundary3D.h"

namespace cnn
{

  RegCNNBoundary3D::RegCNNBoundary3D(const int& nbinsX, // # of wires of pixel map
                     const int& nbinsY,              // # of TDC of pixel map
                     const int& nbinsZ,              // # of TDC of pixel map
                     const double& x_length,              // # of Wire merging
		             const double& y_length,              // # of TDC merging
		             const double& z_length,              // # of TDC merging
                     const float& center_X,
                     const float& center_Y,
                     const float& center_Z):
    fStartBin{0, 0, 0},
    fStopBin{nbinsX, nbinsY, nbinsZ},
    fStart{0, 0, 0},
    fStop{x_length, y_length, z_length},
    fCenter{center_X, center_Y, center_Z}
  {
  }

  bool RegCNNBoundary3D::IsWithin(const float& rel_x, const float& rel_y, const float& rel_z)
  {
    bool inX = (rel_x >= fStart[0]) && (rel_x < fStop[0]);
    bool inY = (rel_y >= fStart[1]) && (rel_y < fStop[1]);
    bool inZ = (rel_z >= fStart[2]) && (rel_z < fStop[2]);

    return inX && inY && inZ;
  }

  std::ostream& operator<<(std::ostream& os, const RegCNNBoundary3D& b)
  {
    os<<"RegCNNBoundary3D actual pixel map coverage (in cm)\n"
      <<b.Length(0)<<", "<<b.Length(1)<<b.Length(2)<<"\n"
      <<"Actual pixel size\n"
      <<b.NBins(0)<<", "<<b.NBins(1)<<", "<<b.NBins(2)<<"\n";

    return os;
  }
}
