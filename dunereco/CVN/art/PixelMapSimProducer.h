////////////////////////////////////////////////////////////////////////
/// \file    PixelMapSimProducer.h
/// \brief   PixelMapSimProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
////////////////////////////////////////////////////////////////////////

#ifndef CVN_PIXELMAPSIMPRODUCER_H
#define CVN_PIXELMAPSIMPRODUCER_H


#include <array>
#include <vector>

// Framework includes
#include "art/Framework/Principal/Handle.h"

#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/SparsePixelMap.h"
#include "dune/CVN/func/Boundary.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace cvn
{
  /// Producer algorithm for PixelMap, input to CVN neural net
  class PixelMapSimProducer
  {
  public:
    PixelMapSimProducer(unsigned int nWire, unsigned int nTdc, double tRes, double threshold = 0.);
    PixelMapSimProducer();

    void SetUnwrapped(unsigned short unwrap){fUnwrapped = unwrap;};
    void SetProtoDUNE(){fProtoDUNE = true;};

    /// Get boundaries for pixel map representation of cluster
    Boundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                            const std::vector< const sim::SimChannel* >& cluster);

    unsigned int NROI(){return fTotHits;};

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
                       const std::vector< art::Ptr< sim::SimChannel > >& slice);
    PixelMap CreateMap(detinfo::DetectorPropertiesData const& detProp,
                       const std::vector< const sim::SimChannel* >& slice);

    PixelMap CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,
                                    const std::vector< const sim::SimChannel* >& cluster,
                                    const Boundary& bound);

  private:
    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    double            fTRes;   ///< Timing resolution for pixel map
    double            fThreshold; ///< charge threshold for each time tick, below which isn't added to pixel map
    unsigned short    fUnwrapped; ///< Use unwrapped pixel maps?
    bool              fProtoDUNE; ///< Do we want to use this for particle extraction from protoDUNE?

    unsigned int fTotHits;  ///<How many ROIs above threshold?

    geo::GeometryCore const* fGeometry;
    std::vector<double> fVDPlane0;
    std::vector<double> fVDPlane1;
    // std::vector<int> fPlane0GapWires;
    // std::vector<int> fPlane1GapWires;

    double _getIntercept(geo::WireID wireid) const;
    void _cacheIntercepts();
  };

}

#endif  // CVN_PixelMapSimProducer_H
