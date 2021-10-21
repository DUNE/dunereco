////////////////////////////////////////////////////////////////////////
// \file    CVNSparseMapper3D_module.cc
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
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/SpacePoint.h"

#include "dune/CVN/art/PixelMapProducer.h"
#include "dune/CVN/func/SparsePixelMap.h"

namespace cvn {

  class CVNSparseMapper3D : public art::EDProducer {
  public:
    explicit CVNSparseMapper3D(fhicl::ParameterSet const& pset);
    ~CVNSparseMapper3D();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    std::string    fSPModuleLabel;   ///< Module label for reconstructed spacepoints
    std::string    fPixelMapLabel;   ///< Instance label for spacepoint pixelmaps
    unsigned short fMinSP;           ///< Minimum number of spacepoints to be converted to pixel map
    PixelMapProducer fProducer;      ///< PixelMapProducer does the heavy lifting

  };



  //.......................................................................
  CVNSparseMapper3D::CVNSparseMapper3D(fhicl::ParameterSet const& pset): EDProducer{pset},
  fSPModuleLabel(pset.get<std::string>      ("SpacePointModuleLabel")),
  fPixelMapLabel(pset.get<std::string>      ("PixelMapLabel")),
  fMinSP(pset.get<unsigned short>           ("MinSP")),
  fProducer()
  {
    produces< std::vector<cvn::SparsePixelMap> >(fPixelMapLabel);
  }

  //......................................................................
  CVNSparseMapper3D::~CVNSparseMapper3D()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNSparseMapper3D::beginJob()
  {}

  //......................................................................
  void CVNSparseMapper3D::endJob()
  {}

  //......................................................................
  void CVNSparseMapper3D::produce(art::Event& evt)
  {
    // Get spacepoints from art event record
    std::vector< art::Ptr< recob::SpacePoint > > splist;
    auto spListHandle = evt.getHandle< std::vector< recob::SpacePoint > >(fSPModuleLabel);
    if (spListHandle)
      art::fill_ptr_vector(splist, spListHandle);
    unsigned short nsp = splist.size();

    // Get assocations from spacepoints to hits
    art::FindManyP<recob::Hit> fmp(spListHandle, evt, fSPModuleLabel);
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit(nsp);
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint

    //Declaring containers for things to be stored in event
    std::unique_ptr< std::vector<cvn::SparsePixelMap> >
      pmCol(new std::vector<cvn::SparsePixelMap>);

    if (nsp > fMinSP) {
      auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
      SparsePixelMap map = fProducer.CreateSparseMap3D(clockData, splist, sp2Hit);
      pmCol->push_back(map);
      mf::LogInfo("CVNSparseMapper3D") << "Created sparse pixel map from "
        << nsp << " spacepoints and " << map.GetNPixels(0) << " hits.";
    }

    evt.put(std::move(pmCol), fPixelMapLabel);
  }

  //----------------------------------------------------------------------

DEFINE_ART_MODULE(cvn::CVNSparseMapper3D)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
