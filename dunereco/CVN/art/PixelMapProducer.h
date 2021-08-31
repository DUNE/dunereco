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
#include "lardataobj/RecoBase/SpacePoint.h"

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
    Boundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                            const std::vector< const recob::Hit* >& cluster);

    /// Function to convert to a global unwrapped wire number
    void GetDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const; 
    void GetDUNEGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                              unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                              unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const;

    void GetDUNE10ktGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                                  unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                  unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const;
    void GetProtoDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const; 
    void GetProtoDUNEGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                   unsigned int& globalWire, double& globalTDC, unsigned int& globalPlane) const;
    // preliminary vert drift 3 view studies 
    void GetDUNEVertDrift3ViewGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const; 


    unsigned int NWire() const {return fNWire;};
    unsigned int NTdc() const {return fNTdc;};
    double TRes() const {return fTRes;};

    PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                       const std::vector< art::Ptr< recob::Hit > >& slice);
    PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                       const std::vector< const recob::Hit* >& slice);

    PixelMap CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,
                                    const std::vector< const recob::Hit* >& cluster,
                                    const Boundary& bound);

    /// Create sparse pixel map for SCN applications
    void GetHitTruth(detinfo::DetectorClocksData const& clockData,
                     art::Ptr<recob::Hit>& hit, std::vector<int>& pdgs, std::vector<int>& tracks,
      std::vector<float>& energies, std::vector<std::string>& processes);
    SparsePixelMap CreateSparseMap2D(detinfo::DetectorClocksData const& clockData,
                                     detinfo::DetectorPropertiesData const& detProp,
                                     std::vector< art::Ptr< recob::Hit> >& cluster, bool usePixelTruth=false);
    SparsePixelMap CreateSparseMap3D(detinfo::DetectorClocksData const& clockData,
                                     std::vector< art::Ptr< recob::SpacePoint> >& sp, std::vector<std::vector<art::Ptr<recob::Hit>>>& hit);

  private:
    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    double            fTRes;   ///< Timing resolution for pixel map
    unsigned short    fUnwrapped; ///< Use unwrapped pixel maps?
    bool              fProtoDUNE; ///< Do we want to use this for particle extraction from protoDUNE?

    geo::GeometryCore const* fGeometry;
    std::vector<double> fVDPlane0;
    std::vector<double> fVDPlane1;
    // std::vector<int> fPlane0GapWires;
    // std::vector<int> fPlane1GapWires;

    double _getIntercept(geo::WireID wireid) const;
    void _cacheIntercepts();
  };

}

#endif  // CVN_PIXELMAPPRODUCER_H
