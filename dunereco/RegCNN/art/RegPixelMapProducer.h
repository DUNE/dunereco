////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMapProducer.h
/// \brief   RegPixelMapProducer for RegCNN modified from PixelMapProducer.h
/// \author  Ilsoo Seong - iseong@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef REGCNN_PIXELMAPPRODUCER_H
#define REGCNN_PIXELMAPPRODUCER_H


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
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/func/RegCNNBoundary.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"

//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

namespace cnn
{
  /// Producer algorithm for RegPixelMap, input to CVN neural net
  class RegPixelMapProducer
  {
  public:
    RegPixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes, int Global);

    /// Get boundaries for pixel map representation of cluster
    RegCNNBoundary DefineBoundary(std::vector< art::Ptr< recob::Hit > >& cluster);
    RegCNNBoundary DefineBoundary(std::vector< art::Ptr< recob::Hit > >& cluster, const float *vtx);

    /// Function to convert to a global unwrapped wire number
    double GetGlobalWire(const geo::WireID& wireID);
    void GetDUNEGlobalWireTDC(const geo::WireID& wireID, float localTDC,
                              double& globalWire, unsigned int& globalPlane, 
                              float& globalTDC);


    unsigned int NWire() const {return fNWire;};
    unsigned int NTdc() const {return fNTdc;};
    double TRes() const {return fTRes;};

    RegPixelMap CreateMap(std::vector< art::Ptr< recob::Hit > >& cluster, art::FindManyP<recob::Wire> fmwire);
    RegPixelMap CreateMap(std::vector< art::Ptr< recob::Hit > >& cluster, 
                          art::FindManyP<recob::Wire> fmwire,
                          const float *vtx);

    RegPixelMap CreateMapGivenBoundary(std::vector< art::Ptr< recob::Hit > >& cluster,
                                    const RegCNNBoundary& bound,
                                    art::FindManyP<recob::Wire> fmwire);

    void ShiftGlobalWire(std::vector< art::Ptr< recob::Hit > >& cluster);

   private:

    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    unsigned int      fTRes;
    int               fGlobalWireMethod;
    double            fOffset[2];
    std::vector<int> hitwireidx; // collect hit wire
    std::vector<int> tmin_each_wire;
    std::vector<int> tmax_each_wire;
    std::vector<float> trms_max_each_wire;

    detinfo::DetectorProperties const* detprop;
    art::ServiceHandle<geo::Geometry> geom;
  };

}

#endif  // REGCNN_PIXELMAPPRODUCER_H
