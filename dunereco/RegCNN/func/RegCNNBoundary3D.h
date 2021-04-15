////////////////////////////////////////////////////////////////////////
/// \file    RegCNNBoundary3D.h
/// \brief   RegCNNBoundary3D for 3D RegCNN PixelMap modified from CVNBoundary.h
/// \author  Wenjie Wu   - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_BOUNDARY3D_H
#define REGCNN_BOUNDARY3D_H

#include <ostream>
#include <vector>


namespace cnn
{

  /// RegCNNBoundary object intended for use with cnn::RegPixelMap.  
  /// Stores mean of wire and TPCs
  class RegCNNBoundary3D
  {

  public:
    /// Create new RegCNNBoundary3D object based on space points
    RegCNNBoundary3D(const int& nbinsX, const int& nbinsY, const int& nbinsZ, 
             const double& x_length, const double& y_length, const double& z_length,
             const float& center_X, const float& center_Y, const float& center_Z);

    RegCNNBoundary3D(){};

    bool IsWithin(const float& rel_x, const float& rel_y, const float& rel_z);

    double Length(const unsigned int& axis) const {return fStop[axis]-fStart[axis];};
    int NBins(const unsigned int& axis) const {return fStopBin[axis]-fStartBin[axis];};
    float StartPos(const unsigned int& axis) const {return fStart[axis];};
    float StopPos(const unsigned int& axis) const {return fStop[axis];};
    float Center(const unsigned int& axis) const {return fCenter[axis];};

    int fStartBin[3]; ///< Minimum TDC in each view, inclusive
    int fStopBin[3];  ///< Maximum TDC in each view, inclusive
    double fStart[3];  ///< Minimum wire, inclusive
    double fStop[3];   ///< Maximum wire, inclusive
    float fCenter[3];

  };

  std::ostream& operator<<(std::ostream& os, const RegCNNBoundary3D& b);
}

#endif  // REGCNN_BOUNDARY3D_H
