////////////////////////////////////////////////////////////////////////
// \file    CVNSparseMapper_module.cc
// \brief   Producer module for creating CVN SparsePixelMap objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

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

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"

#include "dune/CVN/art/PixelMapProducer.h"
#include "dune/CVN/func/SparsePixelMap.h"

namespace cvn {

  class CVNSparseMapper : public art::EDProducer {
  public:
    explicit CVNSparseMapper(fhicl::ParameterSet const& pset);
    ~CVNSparseMapper();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    std::string    fHitsModuleLabel; ///< Module label for input clusters
    std::string    fClusterPMLabel;  ///< Instance label for cluster pixelmaps
    unsigned short fMinClusterHits;  ///< Minimum number of hits for cluster to be converted to pixel map
    bool fIncludePixelTruth;         ///< Whether to include per-pixel ground truth for segmentation
    PixelMapProducer fProducer;      ///< PixelMapProducer does the heavy lifting

  };



  //.......................................................................
  CVNSparseMapper::CVNSparseMapper(fhicl::ParameterSet const& pset): EDProducer{pset},
  fHitsModuleLabel(pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>     ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short>  ("MinClusterHits")),
  fIncludePixelTruth(pset.get<bool>         ("IncludePixelTruth")),
  fProducer()
  {

    produces< std::vector<cvn::SparsePixelMap> >(fClusterPMLabel);

  }

  //......................................................................
  CVNSparseMapper::~CVNSparseMapper()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNSparseMapper::beginJob()
  {}

  //......................................................................
  void CVNSparseMapper::endJob()
  {}

  //......................................................................
  void CVNSparseMapper::produce(art::Event& evt)
  {
    std::vector< art::Ptr< recob::Hit > > hitlist;
    auto hitListHandle = evt.getHandle< std::vector< recob::Hit > >(fHitsModuleLabel);
    if (hitListHandle)
      art::fill_ptr_vector(hitlist, hitListHandle);
    unsigned short nhits = hitlist.size();

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::SparsePixelMap> >
      pmCol(new std::vector<cvn::SparsePixelMap>);

    if (nhits > fMinClusterHits) {
      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockData);
      SparsePixelMap map = fProducer.CreateSparseMap2D(clockData, detProp, hitlist, fIncludePixelTruth);
      mf::LogInfo("CVNSparseMapper") << "Created sparse pixel map with "
        << map.GetNPixels(0) << ", " << map.GetNPixels(1) << ", "
        << map.GetNPixels(2) << " pixels.";
      pmCol->push_back(map);
    }

    evt.put(std::move(pmCol), fClusterPMLabel);
  }

  //----------------------------------------------------------------------

DEFINE_ART_MODULE(cvn::CVNSparseMapper)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
