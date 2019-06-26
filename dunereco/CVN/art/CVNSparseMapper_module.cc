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
    /// Module lablel for input clusters
    std::string    fHitsModuleLabel;

    /// Instance lablel for cluster pixelmaps
    std::string    fClusterPMLabel;

    /// Minimum number of hits for cluster to be converted to pixel map
    unsigned short fMinClusterHits;

    /// PixelMapProducer does the work for us
    PixelMapProducer fProducer;

  };



  //.......................................................................
  CVNSparseMapper::CVNSparseMapper(fhicl::ParameterSet const& pset): EDProducer{pset},
  fHitsModuleLabel(pset.get<std::string>    ("HitsModuleLabel")),
  fClusterPMLabel(pset.get<std::string>     ("ClusterPMLabel")),
  fMinClusterHits(pset.get<unsigned short>  ("MinClusterHits")),
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
    art::Handle< std::vector< recob::Hit > > hitListHandle;
    std::vector< art::Ptr< recob::Hit > > hitlist;
    if (evt.getByLabel(fHitsModuleLabel, hitListHandle))
      art::fill_ptr_vector(hitlist, hitListHandle);
    unsigned short nhits = hitlist.size();

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::SparsePixelMap> >
      pmCol(new std::vector<cvn::SparsePixelMap>);

    if (nhits > fMinClusterHits) {
      pmCol->push_back(fProducer.CreateSparseMap(hitlist));
    }

    evt.put(std::move(pmCol), fClusterPMLabel);
  }

  //----------------------------------------------------------------------

DEFINE_ART_MODULE(cvn::CVNSparseMapper)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







