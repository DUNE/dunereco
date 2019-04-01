////////////////////////////////////////////////////////////////////////
/// \file    RegCNNBoundary.h
/// \brief   RegCNNBoundary for RegCNN PixelMap modified from CVNBoundary.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_BOUNDARY_H
#define REGCNN_BOUNDARY_H

#include <ostream>
#include <vector>


namespace cnn
{


  /// RegCNNBoundary object intended for use with cnn::RegPixelMap.  
  /// Stores mean of wire and TPCs
  class RegCNNBoundary
  {

  public:
    /// Create new RegCNNBoundary object based on number of wires, number of tdcs,
    /// minumum wire and mean tdc in odd and even view.
    RegCNNBoundary(const int& nWire, const int& nTDC, const int& tRes,
             const int& WireMeanX,
             const int& WireMeanY,
             const int& WireMeanZ,
             const int& TDCMeanX,
             const int& TDCMeanY,
             const int& TDCMeanZ);

    RegCNNBoundary(){};

    bool IsWithin(const int& wire, const int& tdc, const unsigned int& view);

    int FirstWire(const unsigned int& view) const {return fFirstWire[view];};
    int LastWire(const unsigned int& view) const {return fLastWire[view];};
    int FirstTDC(const unsigned int& view) const {return fFirstTDC[view];};
    int LastTDC (const unsigned int& view) const {return fLastTDC[view];};



    int fFirstWire[3];  ///< Minimum wire, inclusive
    int fLastWire[3];   ///< Maximum wire, inclusive
    int fFirstTDC[3]; ///< Minimum TDC in each view, inclusive
    int fLastTDC[3];  ///< Maximum TDC in each view, inclusive


  };

  std::ostream& operator<<(std::ostream& os, const RegCNNBoundary& b);
}

#endif  // REGCNN_BOUNDARY_H
