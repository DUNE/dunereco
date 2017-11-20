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
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/Boundary.h"
#include "lardataobj/RecoBase/Hit.h"

//#include "art/Framework/Services/Registry/ServiceHandle.h"
//#include "Geometry/Geometry.h"

namespace cvn
{
  /// Producer algorithm for PixelMap, input to CVN neural net
  class PixelMapProducer
  {
  public:
    PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes);

    void SetUnwrapped(bool unwrap){fUnwrapped = unwrap;};

    /// Get boundaries for pixel map representation of cluster
    Boundary DefineBoundary(std::vector< art::Ptr< recob::Hit > >& cluster);

    /// Function to convert to a global unwrapped wire number
    void GetGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const; 

    unsigned int NWire() const {return fNWire;};
    unsigned int NTdc() const {return fNTdc;};
    double TRes() const {return fTRes;};

    PixelMap CreateMap(std::vector< art::Ptr< recob::Hit > >& slice);

    PixelMap CreateMapGivenBoundary(std::vector< art::Ptr< recob::Hit > >& cluster,
                                    const Boundary& bound);

   private:
    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    double            fTRes;
    bool              fUnwrapped; ///< Use unwrapped pixel maps?
  };

}

#endif  // CVN_PIXELMAPPRODUCER_H
