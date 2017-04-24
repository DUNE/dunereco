////////////////////////////////////////////////////////////////////////
/// \file    Boundary.cxx
/// \brief   Boundary for CVN PixelMap
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <utility>
#include  <cassert>


#include "dune/CVN/func/Boundary.h"

namespace cvn
{

  Boundary::Boundary(const int& nWire,
                     const double& tRes,
                     const int& minWireX,
                     const int& minWireY,
                     const int& minWireZ,
                     const double& centerTDCX,
                     const double& centerTDCY,
                     const double& centerTDCZ):
    fFirstWire{minWireX,
      minWireY,
      minWireZ},
    fLastWire{minWireX + nWire - 1,
        minWireY + nWire - 1,
        minWireZ + nWire - 1},
    fFirstTDC{centerTDCX - tRes,   // For odd nTDC, we will truncate 0.5,
        centerTDCY  - tRes,
        centerTDCZ  - tRes},  // but get it back in LastTDC
    fLastTDC{centerTDCX + tRes,// Recover the trucated 0.5
        centerTDCY + tRes,
        centerTDCZ + tRes}// with nTDC % 2
  
  {
    assert(fLastWire[0] - fFirstWire[0] == nWire - 1);
    assert(fLastWire[1] - fFirstWire[1] == nWire - 1);
    assert(fLastWire[2] - fFirstWire[2] == nWire - 1);
  }

  bool Boundary::IsWithin(const unsigned int& wire, const double& cell, const unsigned int& view)
  {
    bool inWireRcvne = (int) wire >= fFirstWire[view] && (int) wire <= fLastWire[view];
    bool inTDCRcvne = (double) cell >= fFirstTDC[view] &&
                       (double) cell <= fLastTDC[view];
    return inWireRcvne && inTDCRcvne;
  }


  std::ostream& operator<<(std::ostream& os, const Boundary& b)
  {
    os<<"Boundary with "
      <<"(first,last) wire X: (" << b.FirstWire(0)<<", "<< b.LastWire(0)
      <<"(first,last) wire Y: (" << b.FirstWire(1)<<", "<< b.LastWire(1)
      <<"(first,last) wire Z: (" << b.FirstWire(2)<<", "<< b.LastWire(2)
      <<"), (first,last) tdc X: ("<<b.FirstTDC(0)<<", "<<b.LastTDC(0)<<")"
      <<"), (first,last) tdc Y: ("<<b.FirstTDC(1)<<", "<<b.LastTDC(1)<<")"
      <<"), (first,last) tdc Z: ("<<b.FirstTDC(2)<<", "<<b.LastTDC(2)<<")";

    return os;
  }
}
