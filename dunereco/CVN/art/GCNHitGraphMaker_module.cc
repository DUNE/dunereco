////////////////////////////////////////////////////////////////////////
// Class:       GCNHitGraphMaker
// Plugin Type: producer (art v3_05_00)
// File:        GCNHitGraphMaker_module.cc
//
// Generated at Mon Apr 13 11:51:42 2020 by Jeremy Hewes using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/art/PixelMapProducer.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include <memory>

namespace cvn {


  class GCNHitGraphMaker : public art::EDProducer {
  public:
    explicit GCNHitGraphMaker(fhicl::ParameterSet const& p);

    GCNHitGraphMaker(GCNHitGraphMaker const&) = delete;
    GCNHitGraphMaker(GCNHitGraphMaker&&) = delete;
    GCNHitGraphMaker& operator=(GCNHitGraphMaker const&) = delete;
    GCNHitGraphMaker& operator=(GCNHitGraphMaker&&) = delete;

    void configure(fhicl::ParameterSet const& p);

    void produce(art::Event& e) override;

  private:

    std::string fHitModuleLabel;

  };


  GCNHitGraphMaker::GCNHitGraphMaker(fhicl::ParameterSet const& p)
    : EDProducer{p}
  {
    configure(p);
    produces<std::vector<cvn::GCNGraph>>();
    consumes<std::vector<recob::Hit>>(fHitModuleLabel);
  }

  void GCNHitGraphMaker::configure(fhicl::ParameterSet const& p)
  {
    fHitModuleLabel = p.get<std::string>("HitModuleLabel");
  }

  void GCNHitGraphMaker::produce(art::Event& e)
  {
    // Retrieve hits
    std::vector<art::Ptr<recob::Hit>> hits;
    auto hitHandle = e.getHandle<std::vector<recob::Hit>>(fHitModuleLabel);
    if (!hitHandle) {
      throw art::Exception(art::errors::LogicError)
        << "Could not find hits with module label " << fHitModuleLabel;
    }
    art::fill_ptr_vector(hits, hitHandle);

    // Create output
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>(1));

    // Helper classes and services
    cvn::PixelMapProducer pixelUtil;
    art::ServiceHandle<cheat::BackTrackerService> bt;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);

    // Loop over hits
    for (art::Ptr<recob::Hit> hit : hits) {

      // Get global coordinates
      geo::WireID wireid = hit->WireID();
      if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue; // skip dummy TPCs
      unsigned int wire, plane;
      double time;
      pixelUtil.GetDUNE10ktGlobalWireTDC(detProp, wireid.Wire, hit->PeakTime(),
        wireid.Plane, wireid.TPC, wire, plane, time);
      std::vector<float> pos = {(float)plane, (float)wire, (float)time};
      std::vector<float> feats = {hit->Integral(), hit->RMS(), (float)wireid.Wire,
        hit->PeakTime(), (float)wireid.Plane, (float)wireid.TPC};

      // Get true particle
      std::vector<sim::TrackIDE> ides = bt->HitToTrackIDEs(clockData, hit);
      if (ides.size() == 0) continue; // skip hits with no truth
      float trueID = abs(std::max_element(ides.begin(), ides.end(), [](const sim::TrackIDE &lhs,
        const sim::TrackIDE &rhs) { return lhs.energy < rhs.energy; })->trackID);
      std::vector<float> truth = {trueID};

      // Get 
      graphs->at(0).AddNode(pos, feats, truth);
      // std::cout << "Adding node with plane " << plane << ", wire " << wire << ", time " << time
      //   << ", integral " << hit->Integral() << ", RMS " << hit->RMS() << ", raw wire " << wireid.Wire
      //   << ", raw time " << hit->PeakTime() << ", raw plane " << wireid.Plane << ", TPC " << wireid.TPC
      //   << " and true ID " << trueID << std::endl;
    } // for hit

    e.put(std::move(graphs));

  }

  DEFINE_ART_MODULE(GCNHitGraphMaker)

} // namespace cvn
