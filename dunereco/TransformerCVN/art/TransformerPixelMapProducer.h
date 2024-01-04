////////////////////////////////////////////////////////////////////////
/// \file    TransformerPixelMapProducer.h
/// \brief   TransformerPixelMapProducer for TransformerCVN modified from RegPixelMapProducer.h
/// \author  Alejandro Yankelevich - ayankele@uci.edu
////////////////////////////////////////////////////////////////////////

#ifndef TRANSFORMERCVN_PIXELMAPPRODUCER_H
#define TRANSFORMERCVN_PIXELMAPPRODUCER_H


#include <array>
#include <vector>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "dunereco/TransformerCVN/func/TransformerPixelMap.h"
#include "dunereco/RegCNN/func/RegCNNBoundary.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}

//#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

namespace cnn
{
  /// Producer algorithm for TransformerPixelMap, input to CVN neural net
  class TransformerPixelMapProducer
  {
  public:
    TransformerPixelMapProducer(unsigned int nWire, unsigned int wRes, unsigned int nTdc, double tRes, int Global,
            bool ProngOnly, bool ByHit);

    /// Get boundaries for pixel map representation of cluster
    RegCNNBoundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                  std::vector< art::Ptr< recob::Hit > > const& cluster);
    RegCNNBoundary DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                  std::vector< art::Ptr< recob::Hit > > const& cluster,
                                  const std::vector<float> &vtx);

    /// Function to convert to a global unwrapped wire number
    double GetGlobalWire(const geo::WireID& wireID);
    void GetDUNEGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                              const geo::WireID& wireID, double localTDC,
                              unsigned int& globalWire, unsigned int& globalPlane, 
                              double& globalTDC);


    unsigned int NWire() const {return fNWire;};
    unsigned int NTdc() const {return fNTdc;};
    double TRes() const {return fTRes;};
    double WRes() const {return fWRes;};

    TransformerPixelMap CreateMap(detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp,
                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                          art::FindManyP<recob::Wire> const& fmwire);

    TransformerPixelMap CreateMap(detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp,
                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                          art::FindManyP<recob::Wire> const& fmwire,
                          const std::vector<float> &vtx);

    TransformerPixelMap CreateMapGivenBoundary(detinfo::DetectorClocksData const& clockData,
                                       detinfo::DetectorPropertiesData const& detProp,
                                       std::vector< art::Ptr< recob::Hit > > const& cluster,
                                       const RegCNNBoundary& bound,
                                       art::FindManyP<recob::Wire> const& fmwire);

    // tag prong by track
    TransformerPixelMap CreateMap(detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp,
                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                          art::FindManyP<recob::Wire> const& fmwire,
                          art::FindManyP<recob::Track> const& fmtrkhit,
                          const std::vector<float> &vtx);
    TransformerPixelMap CreateMapGivenBoundaryByHit(detinfo::DetectorClocksData const& clockData,
                                       detinfo::DetectorPropertiesData const& detProp,
                                       std::vector< art::Ptr< recob::Hit > > const& cluster,
                                       const RegCNNBoundary& bound,
                                       art::FindManyP<recob::Wire> const& fmwire,
                                       art::FindManyP<recob::Track> const& fmtrkhit,
                                       const bool& ProngOnly);

    // tag prong by shower
    TransformerPixelMap CreateMap(detinfo::DetectorClocksData const& clockData,
                          detinfo::DetectorPropertiesData const& detProp,
                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                          art::FindManyP<recob::Wire> const& fmwire,
                          art::FindManyP<recob::Shower> const& fmshwhit,
                          const std::vector<float> &vtx);
    TransformerPixelMap CreateMapGivenBoundaryByHit(detinfo::DetectorClocksData const& clockData,
                                       detinfo::DetectorPropertiesData const& detProp,
                                       std::vector< art::Ptr< recob::Hit > > const& cluster,
                                       const RegCNNBoundary& bound,
                                       art::FindManyP<recob::Wire> const& fmwire,
                                       art::FindManyP<recob::Shower> const& fmshwhit,
                                       const bool& ProngOnly);

    void ShiftGlobalWire(std::vector< art::Ptr< recob::Hit > > const& cluster);

   private:

    unsigned int      fNWire;  ///< Number of wires, length for pixel maps
    unsigned int      fWRes;
    unsigned int      fNTdc;   ///< Number of tdcs, width of pixel map
    unsigned int      fTRes;
    int               fGlobalWireMethod;
    bool              fProngOnly;
    bool              fByHit;
    double            fOffset[2];
    std::vector<int> hitwireidx; // collect hit wire
    std::vector<int> tmin_each_wire;
    std::vector<int> tmax_each_wire;
    std::vector<float> trms_max_each_wire;

    art::ServiceHandle<geo::Geometry> geom;
  };

}

#endif  // TRANSFORMERCVN_PIXELMAPPRODUCER_H
