////////////////////////////////////////////////////////////////////////
/// \file    Boundary.h
/// \brief   Boundary for CVN PixelMap
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_BOUNDARY_H
#define CVN_BOUNDARY_H

#include <ostream>
#include <vector>






namespace cvn
{


  /// Boundary object intended for use with cvn::PixelMap.  Stores first and
  /// last wires, as well as first and last cell for even and odd view.
  /// CVN doesn't carefully define X/Y view, but instead simply uses
  /// odd/even wire (wire%2) as a proxy.
  class Boundary
  {

  public:
    /// Create new Boundary object based on number of wires, number of cells,
    /// minumum wire and mean cell in odd and even view.
    Boundary(const int& nWire, const double& tRes,
             const int& minWireX,
             const int& minWireY,
             const int& minWireZ,
             const double& centerTDCX,
             const double& centerTDCY,
             const double& centerTDCZ);

    Boundary(){};

    bool IsWithin(const unsigned int& wire, const double& cell, const unsigned int& view);

    int FirstWire(const unsigned int& view) const {return fFirstWire[view];};
    int LastWire(const unsigned int& view) const {return fLastWire[view];};
    double FirstTDC(const unsigned int& view) const {return fFirstTDC[view];};
    double LastTDC (const unsigned int& view) const {return fLastTDC[view];};



  private:
    int fFirstWire[3];  ///< Minimum wire, inclusive
    int fLastWire[3];   ///< Maximum wire, inclusive
    double fFirstTDC[3]; ///< Minimum cell in each view, inclusive
    double fLastTDC[3];  ///< Maximum cell in each view, inclusive


  };

  std::ostream& operator<<(std::ostream& os, const Boundary& b);
}

#endif  // CVN_BOUNDARY_H
