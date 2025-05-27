////////////////////////////////////////////////////////////////////////
// Class:       HitFilter
// Plugin Type: producer (Unknown Unknown)
// File:        HitFilter_module.cc
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
#include "lardataobj/RecoBase/Hit.h"

#include <memory>

class HitFilter;


class HitFilter : public art::EDProducer {
public:
  explicit HitFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitFilter(HitFilter const&) = delete;
  HitFilter(HitFilter&&) = delete;
  HitFilter& operator=(HitFilter const&) = delete;
  HitFilter& operator=(HitFilter&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  // Declare member data here.
  std::vector<int> fValidTPCs;
  std::string fInputHitLabel;
};


HitFilter::HitFilter(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fValidTPCs(p.get<std::vector<int>>("ValidTPCs", {})),
    fInputHitLabel(p.get<std::string>("InputHitLabel", "hitpdune"))
  // More initializers here.
{
  // Call appropriate produces<>() functions here.
  produces<std::vector<recob::Hit>>();
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void HitFilter::produce(art::Event& e)
{
  // Implementation of required member function here.
  auto output_hits = std::make_unique<std::vector<recob::Hit>>();
  auto input_hits = e.getValidHandle<std::vector<recob::Hit>>(fInputHitLabel);
  for (const auto & hit : (*input_hits)) {
    if (std::find(fValidTPCs.begin(), fValidTPCs.end(), hit.WireID().TPC)
        == fValidTPCs.end()) continue;
    output_hits->push_back(hit);
  }

  e.put(std::move(output_hits));
}

DEFINE_ART_MODULE(HitFilter)
