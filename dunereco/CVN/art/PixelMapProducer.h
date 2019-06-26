////////////////////////////////////////////////////////////////////////
/// \file    PixelMapProducer.h
/// \brief   PixelMapProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_PIXELMAPPRODUCER_H
#define CVN_PIXELMAPPRODUCER_H


#include <array>
#include <vector>

// Framework includes
#include "art/Framework/Principal/Handle.h"

#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/SparsePixelMap.h"
#include "dune/CVN/func/Boundary.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace cvn
{
  /// Producer algorithm for PixelMap, input to CVN neural net
  class PixelMapProducer
  {
  public:
    PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes);
    PixelMapProducer();

    void SetUnwrapped(unsigned short unwrap){fUnwrapped = unwrap;};
    void SetProtoDUNE(){fProtoDUNE = true;};

    /// Get boundaries for pixel map representation of cluster
    Boundary DefineBoundary(const std::vector< const recob::Hit* >& cluster);

    /// Function to convert to a global unwrapped wire number
    void GetDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const; 
    void GetDUNEGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                              unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const;
    void GetDUNE10ktGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                  unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const;
    void GetProtoDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const; 

    unsigned int NWire() const {return fNWire;};
    unsigned int NTdc() const {return fNTdc;};
    double TRes() const {return fTRes;};

    PixelMap CreateMap(const std::vector< art::Ptr< recob::Hit > >& slice);
    PixelMap CreateMap(const std::vector< const recob::Hit* >& slice);

    PixelMap CreateMapGivenBoundary(const std::vector< const recob::Hit* >& cluster,
                                    const Boundary& bound);

    /// Create sparse pixel map for SCN applications
    SparsePixelMap CreateSparseMap(std::vector< art::Ptr< recob::Hit> >& cluster);

  private:
    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    double            fTRes;   ///< Timing resolution for pixel map
    unsigned short    fUnwrapped; ///< Use unwrapped pixel maps?
    bool              fProtoDUNE; ///< Do we want to use this for particle extraction from protoDUNE?

    geo::GeometryCore const* fGeometry;
    detinfo::DetectorProperties const* fDetProp;
  };

}

#endif  // CVN_PIXELMAPPRODUCER_H
