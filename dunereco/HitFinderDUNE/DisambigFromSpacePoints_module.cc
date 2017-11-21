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
  // these types to be replaced with use of feature proposed in redmine #12602
  typedef std::map< unsigned int, std::vector< size_t > > plane_keymap;
  typedef std::map< unsigned int, plane_keymap > tpc_plane_keymap;
  typedef std::map< unsigned int, tpc_plane_keymap > cryo_tpc_plane_keymap;

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
    int runOnSpacePoints(
            const std::vector< art::Ptr<recob::Hit> > & eventHits,
            const art::FindManyP< recob::SpacePoint > & spFromHit,
            std::unordered_map< size_t, geo::WireID > & assignments,
            cryo_tpc_plane_keymap & indHits, cryo_tpc_plane_keymap & unassigned) const;

    int resolveUnassigned(
            std::unordered_map< size_t, geo::WireID > & assignments,
            const std::vector< art::Ptr<recob::Hit> > & eventHits,
            const cryo_tpc_plane_keymap & indHits, const cryo_tpc_plane_keymap & unassigned) const;

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

    // all hits in the collection
    std::vector< art::Ptr<recob::Hit> > eventHits;
    art::fill_ptr_vector(eventHits, hitsHandle);
    art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSpModuleLabel);
    for (size_t i = 0; i < eventHits.size(); ++i)
	{
		art::Ptr<recob::Hit> hit = eventHits[i];
		//if (hit->SignalType() == geo::kCollection)
		if (hit->View() == geo::kZ)
		{
			//std::cout << "nc:" << spFromHit.at(i).size() << std::endl;
			continue;
		}
		else if (hit->View() == geo::kU)
		{
			std::cout << "ni1:";
		
		}
		else if (hit->View() == geo::kV)
		{
			std::cout << "ni2:";
		}
		std::cout << spFromHit.at(i).size() << ", amp:" << hit->PeakAmplitude();
		std::cout << std::endl;
	}

	art::FindManyP< recob::Hit > hitFromSp(spHandle, evt, fSpModuleLabel);
	size_t nz = 0, nu = 0, nv = 0;
	for (size_t i = 0; i < spHandle->size(); ++i)
	{
		for (auto const & hit : hitFromSp.at(i))
		{
			if (hit->View() == geo::kZ) { ++nz; }
			else if (hit->View() == geo::kU) { ++nu; }
			else if (hit->View() == geo::kV) { ++nv; }
		}
	}

    cryo_tpc_plane_keymap indHits, unassignedHits;
    std::unordered_map< size_t, geo::WireID > hitToWire;

    runOnSpacePoints(eventHits, spFromHit, hitToWire, indHits, unassignedHits);

    resolveUnassigned(hitToWire, eventHits, indHits, unassignedHits);

    for (auto const & hw : hitToWire)
    {
        size_t key = hw.first;
        geo::WireID wid = hw.second;

        recob::HitCreator new_hit(*(eventHits[key]), wid);

        art::Ptr<recob::Wire> wire = channelHitWires.at(key);
        art::Ptr<raw::RawDigit> rawdigits = channelHitRawDigits.at(key);

        hcol.emplace_back(new_hit.move(), wire, rawdigits);
    }

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

  int DisambigFromSpacePoints::runOnSpacePoints(
    const std::vector< art::Ptr<recob::Hit> > & eventHits,
    const art::FindManyP< recob::SpacePoint > & spFromHit,
    std::unordered_map< size_t, geo::WireID > & assignments,
    cryo_tpc_plane_keymap & indHits, cryo_tpc_plane_keymap & unassigned) const
  {
    int nUnassigned = 0;
    for (size_t i = 0; i < eventHits.size(); ++i)
    {
        const art::Ptr<recob::Hit> & hit = eventHits[i];
        std::vector<geo::WireID> cwids = fGeom->ChannelToWire(hit->Channel());
        if (cwids.empty()) { continue; } // add warning here
        if (hit->SignalType() == geo::kCollection)
        {
            assignments[hit.key()] = cwids.front();
        }
        else
        {
            geo::WireID id = hit->WireID();
            size_t cryo = id.Cryostat, tpc = id.TPC, plane = id.Plane;

            if (spFromHit.size() == 0)
            {
                unassigned[cryo][tpc][plane].push_back(hit.key());
                ++nUnassigned;
            }
            else
            {
                std::unordered_map< size_t, geo::WireID > tpcBestWire;
                std::unordered_map< size_t, size_t > tpcScore;

                for (const auto & sp : spFromHit.at(hit.key()))
                {
                    geo::TPCID sp_tpcid = fGeom->FindTPCAtPosition(sp->XYZ());
                    if (!sp_tpcid.isValid) { continue; }

                    const float max_dw = 1.0; // max dist to wire [wire pitch]
                    for (size_t w = 0; w < cwids.size(); w++)
                    {
                        if (cwids[w].TPC != sp_tpcid.TPC) { continue; } // not that side of APA

                        float sp_wire = fGeom->WireCoordinate(sp->XYZ()[1], sp->XYZ()[2], plane, sp_tpcid.TPC, cryo);

                        float dw = std::fabs(sp_wire - cwids[w].Wire);
                        if (dw < max_dw)
                        {
                            tpcBestWire[sp_tpcid.TPC] = cwids[w];
                            tpcScore[sp_tpcid.TPC]++;
                        }
                    }
                }
                if (!tpcScore.empty())
                {
                    geo::WireID bestId;
                    size_t maxScore = 0;
                    for (const auto & score : tpcScore)
                    {
                        if (score.second > maxScore)
                        {
                            maxScore = score.second;
                            bestId = tpcBestWire[score.first];
                        }
                    }
                    indHits[cryo][bestId.TPC][plane].push_back(hit.key());
                }
                else
                {
                    std::cout << "***** Did not find matching wire. *****" << std::endl;
                    unassigned[cryo][tpc][plane].push_back(hit.key());
                    ++nUnassigned;
                }
            }
        }
    }

    return nUnassigned;
  }

  int DisambigFromSpacePoints::resolveUnassigned(
    std::unordered_map< size_t, geo::WireID > & assignments,
    const std::vector< art::Ptr<recob::Hit> > & eventHits,
    const cryo_tpc_plane_keymap & allIndHits, const cryo_tpc_plane_keymap & unassigned) const
  {
    int nLeftInPlace = 0;

    return nLeftInPlace;
  }

  DEFINE_ART_MODULE(DisambigFromSpacePoints)

} // end of dune namespace
#endif 

