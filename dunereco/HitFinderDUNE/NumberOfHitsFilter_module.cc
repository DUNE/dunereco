////////////////////////////////////////////////////////////////////////
//
//  Very simple module to filter out events with too many hits
//  - Designed to prevent jobs hanging on very large events
//
//  Leigh Whitehead - leigh.howard.whitehead@cern.ch
//
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Hit.h"

namespace hit{
  class NumberOfHitsFilter;
}

class hit::NumberOfHitsFilter : public art::EDFilter{
public:
  
  explicit NumberOfHitsFilter(fhicl::ParameterSet const& pset);
  virtual ~NumberOfHitsFilter();
  
  void beginJob();
  bool filter(art::Event& evt);
  void endJob();
  
private:

  bool fLimitPerTPC;
  unsigned int fHitLimit;
  std::string fHitModule;  
};
  
//-----------------------------------------------------------------------
hit::NumberOfHitsFilter::NumberOfHitsFilter(fhicl::ParameterSet const& pset):
  EDFilter(pset)
{
  fLimitPerTPC = pset.get<bool>("LimitPerTPC");
  fHitLimit = pset.get<unsigned int>("HitLimit");
  fHitModule = pset.get<std::string>("HitModule");
}

//-----------------------------------------------------------------------
hit::NumberOfHitsFilter::~NumberOfHitsFilter(){}

//-----------------------------------------------------------------------
void hit::NumberOfHitsFilter::beginJob() {}

//-----------------------------------------------------------------------
bool hit::NumberOfHitsFilter::filter(art::Event& evt){

  // Get the hit collection from the event 
  auto allHits = evt.getValidHandle<std::vector<recob::Hit> >(fHitModule);  
 
  bool result = true;

  if(fLimitPerTPC){
    // Find the number of hits per TPC and then filter based on a large value
    std::map<unsigned int,unsigned int> hitsPerTPC;

    for(auto const hit : *allHits){
      hitsPerTPC[hit.WireID().TPC]++;
    }
    
    for(auto const m:  hitsPerTPC){
      if(m.second > fHitLimit){
        result = false;
        break;
      }
    }
  }
  else{
    // This is the simplest thing we can do, just cut on the total number of hits
    if (allHits->size() > fHitLimit){
      result = false;
    }
  }


  return result;
}

void hit::NumberOfHitsFilter::endJob() {}
 
DEFINE_ART_MODULE(hit::NumberOfHitsFilter)
