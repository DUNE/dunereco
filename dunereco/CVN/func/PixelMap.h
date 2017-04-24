////////////////////////////////////////////////////////////////////////
/// \file    PixelMap.h
/// \brief   PixelMap for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_PIXELMAP_H
#define CVN_PIXELMAP_H

#include <ostream>
#include <vector>

#include "dune/CVN/func/Boundary.h"
#include "dune/CVN/func/HitType.h"
#include "TH2F.h"

namespace cvn
{


  /// PixelMap, basic input to CVN neural net
  class PixelMap
  {
  public:
    PixelMap(unsigned int nWire, unsigned int nTdc, const Boundary& bound);
    PixelMap(){};

    /// Length in wires
    unsigned int NWire() const {return fNWire;};

    /// Width in tdcs
    unsigned int NTdc() const {return fNTdc;};

    /// Total number of pixels in map
    unsigned int NPixel() const {return fPE.size();};

    /// Map boundary
    Boundary Bound() const {return fBound;};

    /// Number of inputs for the neural net
    unsigned int NInput() const {return NPixel();};

    void FillInputVector(float* input) const;

    /// Add a hit to the map if it is contained within the wire, tdc rcvne
    /// Could be expanded later to add to overflow accordingly.
    void Add(const unsigned int& wire, const double& tdc,  const unsigned int& view, const double& pe);


    /// Take global wire, tdc (detector) and return index in fPE vector
    unsigned int GlobalToIndex(const unsigned int& wire,
                               const double& tdc,
                               const unsigned int& view)
      ;

    /// Take local wire, tdc (within map) and return index in fPE vector
    unsigned int LocalToIndex(const unsigned int& wire,
                              const unsigned int& tdc)
      const;

  /// Take global wire, tdc (detector) and return index in fPE vector
    unsigned int GlobalToIndexSingle(const unsigned int& wire,
                                     const double& tdc,
                                     const unsigned int& view)
;

    /// Draw pixel map to the screen.  This is pretty hokey and the aspect ratio
    /// is totally unrealistic.
    void Print();

    /// Return the pixel map as a 2D histogram for visualization.
    TH2F* ToTH2() const;
    TH2F* ToLabTH2() const;
    TH2F* SingleViewToTH2(const unsigned int& view) const;

    unsigned int      fNWire;  ///< Number of wires, length of pixel map
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    std::vector<float>   fPE;   ///< Vector of PE measurements for pixels
    std::vector<float>   fPEX;  ///< Vector of X PE measurements for pixels
    std::vector<float>   fPEY;  ///< Vector of Y PE measurements for pixels
    std::vector<float>   fPEZ;  ///< Vector of Y PE measurements for pixels
    std::vector<double>  fPur;  ///< Vector of purity for pixels
    std::vector<double>  fPurX; ///< Vector of X purity for pixels
    std::vector<double>  fPurY; ///< Vector of Y purity for pixels
    std::vector<double>  fPurZ; ///< Vector of Y purity for pixels
    std::vector<HitType> fLab;  ///< Vector of Truth labels for pixels
    std::vector<HitType> fLabX; ///< Vector of X Truth labels for pixels
    std::vector<HitType> fLabY; ///< Vector of Y Truth labels for pixels
    std::vector<HitType> fLabZ; ///< Vector of Y Truth labels for pixels

    Boundary          fBound;    //< Boundary of pixel map

  };

  std::ostream& operator<<(std::ostream& os, const PixelMap& m);

}

#endif  // CVN_PIXELMAP_H
