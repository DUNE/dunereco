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

		fhicl::Atom<art::InputTag> HitModuleLabel { Name("HitModuleLabel"), Comment("HitPointSolver label.") };
		fhicl::Atom<art::InputTag> SpModuleLabel { Name("SpModuleLabel"), Comment("SpacePointSolver label.") };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit DisambigFromSpacePoints(Parameters const& config);

    void produce(art::Event& evt);

  private:
    geo::GeometryCore const* fGeom;
	art::InputTag fHitModuleLabel;
    art::InputTag fSpModuleLabel;
  };

  // implementation

  DisambigFromSpacePoints::DisambigFromSpacePoints(DisambigFromSpacePoints::Parameters const& config) :
    fHitModuleLabel(config().HitModuleLabel()),
    fSpModuleLabel(config().SpModuleLabel())
  {
    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(*this);
    fGeom = &*(art::ServiceHandle<geo::Geometry>());
  }

  void DisambigFromSpacePoints::produce(art::Event& evt)
{
	auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
	auto spHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpModuleLabel);

    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    art::FindOneP<raw::RawDigit> channelHitRawDigits(hitsHandle, evt, fHitModuleLabel);
    art::FindOneP<recob::Wire>   channelHitWires    (hitsHandle, evt, fHitModuleLabel);

    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(*this, evt,
       channelHitWires.isValid(), // doWireAssns
       channelHitRawDigits.isValid() // doRawDigitAssns
      );

    // make hits collection
    std::vector< art::Ptr<recob::Hit> > hits;
    art::fill_ptr_vector(hits, hitsHandle);

	art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSpModuleLabel);
	size_t totz = 0, totu = 0, totv = 0;
	for (size_t i = 0; i < hits.size(); ++i)
	{
		art::Ptr<recob::Hit> hit = hits[i];
		//if (hit->SignalType() == geo::kCollection)
		if (hit->View() == geo::kZ)
		{
			totz += spFromHit.at(i).size();
		}
		else if (hit->View() == geo::kU)
		{
			totu += spFromHit.at(i).size();
		}
		else if (hit->View() == geo::kV)
		{
			totv += spFromHit.at(i).size();
		}
	}

	art::FindManyP< recob::Hit > hitFromSp(spHandle, evt, fSpModuleLabel);
	size_t nz = 0, nu = 0, nv = 0;
	for (size_t i = 0; i < spHandle->size(); ++i)
	{
		for (auto const & hit : hitFromSp.at(i))
		{
			//if (hit->SignalType() == geo::kCollection) { ++nc; }
			//else { ++ni; }
			if (hit->View() == geo::kZ)
			{
				nz += spFromHit.at(i).size();
			}
			else if (hit->View() == geo::kU)
			{
				nu += spFromHit.at(i).size();
			}
			else if (hit->View() == geo::kV)
			{
				nv += spFromHit.at(i).size();
			}
		}
	}
	std::cout << "totz:" << totz << " totu:" << totu << " totv:" << totv << std::endl;
	std::cout << "nz:" << nz << " nu:" << nu << " nv:" << nv << std::endl;

/*    
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
*/
    // put the hit collection and associations into the event
    hcol.put_into(evt);
  }

  DEFINE_ART_MODULE(DisambigFromSpacePoints)

} // end of dune namespace
#endif 

