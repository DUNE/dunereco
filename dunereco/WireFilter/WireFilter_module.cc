////////////////////////////////////////////////////////////////////////
// Class:       WireFilter
// Plugin Type: producer (Unknown Unknown)
// File:        WireFilter_module.cc
//
// Generated at Wed May 21 22:29:59 2025 by root using cetskelgen
// from cetlib version 3.18.02.
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
#include "lardataobj/RecoBase/Wire.h"
#include "larcore/Geometry/WireReadout.h"
#include <memory>

class WireFilter;


class WireFilter : public art::EDProducer {
public:
  explicit WireFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireFilter(WireFilter const&) = delete;
  WireFilter(WireFilter&&) = delete;
  WireFilter& operator=(WireFilter const&) = delete;
  WireFilter& operator=(WireFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::vector<int> fValidTPCs;
  std::string fInputWireLabel;
};


WireFilter::WireFilter(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fValidTPCs(p.get<std::vector<int>>("ValidTPCs", {})),
    fInputWireLabel(p.get<std::string>("InputWireLabel", "wclsdatahd"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Wire>>();
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void WireFilter::produce(art::Event& e)
{
  
  auto const& wireReadoutGeom = art::ServiceHandle<geo::WireReadout const>()->Get();

  // Implementation of required member function here.
  auto output_wires = std::make_unique<std::vector<recob::Wire>>();
  auto input_wires = e.getValidHandle<std::vector<recob::Wire>>(fInputWireLabel);
  for (const auto & wire : (*input_wires)) {
    raw::ChannelID_t channel = wire.Channel();

    std::vector<geo::WireID> wids = wireReadoutGeom.ChannelToWire(channel);
    bool found = false;
    for (auto & wid : wids) {
      // std::cout << wids.size() << " " << wid.TPC << std::endl;
      if (std::find(fValidTPCs.begin(), fValidTPCs.end(), wid.TPC) != fValidTPCs.end()) {
        found = true;
        break;
      }
    }
    if (!found) continue;
    output_wires->push_back(wire);
  }

  e.put(std::move(output_wires));
}

DEFINE_ART_MODULE(WireFilter)
