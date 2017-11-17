#ifndef DISAMBIGFROMSP__CC
#define DISAMBIGFROMSP__CC

////////////////////////////////////////////////////////////////////////
//
// DisambigFromSpacePoints module
//
// R.Sulej
//
// Just look at SpacePoints created with SpacePointSolver and put hits
// in the right wires.
//
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <utility> 
#include <memory>  
#include <iostream>

// Framework includes
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


namespace dune {
  class DisambigFromSpacePoints : public art::EDProducer {
  public:
	struct Config {
		using Name = fhicl::Name;
		using Comment = fhicl::Comment;

		fhicl::Atom<art::InputTag> SpModuleLabel { Name("SpModuleLabel"), Comment("SpacePointSolver label.") };
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit DisambigFromSpacePoints(Parameters const& config);

    void produce(art::Event& evt);

  private:
    geo::GeometryCore const* fGeom;
    art::InputTag fSpModuleLabel;
  };

  // implementation

  DisambigFromSpacePoints::DisambigFromSpacePoints(DisambigFromSpacePoints::Parameters const& config) :
    fSpModuleLabel(config().ChanHitLabel()
  {
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
    fGeom = &*(art::ServiceHandle<geo::Geometry>());
  }

  void DisambigFromSpacePoints::produce(art::Event& evt) {
    art::Handle< std::vector<recob::Hit> > ChannelHits;
    evt.getByLabel(fChanHitLabel, ChannelHits);

    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> ChannelHitRawDigits(ChannelHits, evt, fChanHitLabel);
    art::FindOneP<recob::Wire>   ChannelHitWires    (ChannelHits, evt, fChanHitLabel);

    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this, evt,
       ChannelHitWires.isValid(), // doWireAssns
       ChannelHitRawDigits.isValid() // doRawDigitAssns
      );
    
    // Make hits collection
    std::vector< art::Ptr<recob::Hit> >  ChHits;
    art::fill_ptr_vector(ChHits, ChannelHits);
      
    for( size_t h = 0; h < ChHits.size(); h++ ) {
      std::vector<geo::WireID> cwids = geom->ChannelToWire(ChHits[h]->Channel());
      for(size_t w = 0; w < cwids.size(); w++) {
        art::Ptr<recob::Hit>  hit = ChHits[h];
	      geo::WireID           wid = cwids[w];
        recob::HitCreator repeated_hit(*hit, wid);
        art::Ptr<recob::Hit>::key_type hit_index = hit.key();
	      art::Ptr<recob::Wire> wire = ChannelHitWires.at(hit_index);
	      art::Ptr<raw::RawDigit> rawdigits = ChannelHitRawDigits.at(hit_index);

	      hcol.emplace_back(repeated_hit.move(), wire, rawdigits);
      }
    }

    // put the hit collection and associations into the event
    hcol.put_into(evt);    
  }

  DEFINE_ART_MODULE(DisambigFromSpacePoints)

} // end of dune namespace
#endif 

