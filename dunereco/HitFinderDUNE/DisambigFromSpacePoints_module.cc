#ifndef DISAMBIGFROMSP__CC
#define DISAMBIGFROMSP__CC

////////////////////////////////////////////////////////////////////////
//
// DisambigFromSpacePoints module
//
// R.Sulej
//
// Just look at SpacePoints created with SpacePointSolver and put hits
// in the right wires. Resolve hits undisambiguated by SpacePoints using
// assignments of neighboring hits.
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
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes
#include "TTree.h"

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
        fhicl::Sequence<size_t> ExcludeTPCs { Name("ExcludeTPCs"), Comment("TPC inndexes where hits are not allowed.") };
        fhicl::Atom<bool> UseNeighbors { Name("UseNeighbors"), Comment("Use neighboring hits to complete hits unresolved with spacepoints.") };
        fhicl::Atom<size_t> NumNeighbors { Name("NumNeighbors"), Comment("Number of neighboring hits to complete hits unresolved with spacepoints.") };
        fhicl::Atom<float> MaxDistance { Name("MaxDistance"), Comment("Distance [cm] used to complete hits unresolved with spacepoints.") };
        fhicl::Atom<std::string> MoveLeftovers { Name("MoveLeftovers"), Comment("Mode of dealing with undisambiguated hits.") };
        fhicl::Atom<bool> MonitoringPlots { Name("MonitoringPlots"), Comment("Create histograms of no. of unresolved hits at eacch stage, per plane.") };
    };
    using Parameters = art::EDProducer::Table<Config>;

    explicit DisambigFromSpacePoints(Parameters const& config);

    void produce(art::Event& evt);

  private:
    int runOnSpacePoints(
            const std::vector< art::Ptr<recob::Hit> > & eventHits,
            const art::FindManyP< recob::SpacePoint > & spFromHit,
            const std::unordered_map< size_t, size_t > & spToTPC,
            std::unordered_map< size_t, geo::WireID > & assignments,
            cryo_tpc_plane_keymap & indHits,
            std::vector<size_t> & unassigned
            );

    int resolveUnassigned(
            detinfo::DetectorPropertiesData const& detProp,
            std::unordered_map< size_t, geo::WireID > & assignments,
            const std::vector< art::Ptr<recob::Hit> > & eventHits,
            cryo_tpc_plane_keymap & indHits,
            std::vector<size_t> & unassigned,
            size_t nNeighbors
            );

    void assignFirstAllowedWire(
            std::unordered_map< size_t, geo::WireID > & assignments,
            const std::vector< art::Ptr<recob::Hit> > & eventHits,
            const std::vector<size_t> & unassigned
            ) const;

    void assignEveryAllowedWire(
            std::unordered_map< size_t, std::vector<geo::WireID> > & assignments,
            const std::vector< art::Ptr<recob::Hit> > & eventHits,
            const std::vector<size_t> & unassigned
            ) const;


    geo::GeometryCore const* fGeom;

    int fRun, fEvent;
    int fNHits[3];                // n all hits in each plane
    int fNMissedBySpacePoints[2]; // n hits unresolved by SpacePoints in induction planes
    int fNMissedByNeighbors[2];   // n hits unresolved by using neighboring hits
    TTree *fTree;

    const bool fMonitoringPlots;
    const bool fUseNeighbors;
    const size_t fNumNeighbors;
    const float fMaxDistance;
    const std::string fMoveLeftovers;
    const std::vector< size_t > fExcludeTPCs;
    const art::InputTag fHitModuleLabel;
    const art::InputTag fSpModuleLabel;
  };

  DisambigFromSpacePoints::DisambigFromSpacePoints(DisambigFromSpacePoints::Parameters const& config) :
    EDProducer(config),
    fTree(0),
    fMonitoringPlots(config().MonitoringPlots()),
    fUseNeighbors(config().UseNeighbors()),
    fNumNeighbors(config().NumNeighbors()),
    fMaxDistance(config().MaxDistance()),
    fMoveLeftovers(config().MoveLeftovers()),
    fExcludeTPCs(config().ExcludeTPCs()),
    fHitModuleLabel(config().HitModuleLabel()),
    fSpModuleLabel(config().SpModuleLabel())
  {
    if (fNumNeighbors < 1)
    {
        throw cet::exception("DisambigFromSpacePoints") << "NumNeighbors should be at least 1." << std::endl;
    }

    // let HitCollectionCreator declare that we are going to produce
    // hits and associations with wires and raw digits
    // (with no particular product label)
    recob::HitCollectionCreator::declare_products(producesCollector(), "", true, false);

    // will also copy associations of SpacePoints to original hits
    produces<art::Assns<recob::Hit, recob::SpacePoint>>();

    fGeom = &*(art::ServiceHandle<geo::Geometry>());

    if (fMonitoringPlots)
    {
        art::ServiceHandle<art::TFileService> tfs;
        fTree = tfs->make<TTree>("hitstats", "Unresolved hits statistics");
        fTree->Branch("fRun", &fRun, "fRun/I");
        fTree->Branch("fEvent", &fEvent, "fEvent/I");
        fTree->Branch("fNHits", fNHits, "fNHits[3]/I");
        fTree->Branch("fNMissedBySpacePoints", fNMissedBySpacePoints, "fNMissedBySpacePoints[2]/I");
        fTree->Branch("fNMissedByNeighbors", fNMissedByNeighbors, "fNMissedByNeighbors[2]/I");
    }
  }

  void DisambigFromSpacePoints::produce(art::Event& evt)
  {
    fRun = evt.run();
    fEvent = evt.id().event();

    auto hitsHandle = evt.getValidHandle< std::vector<recob::Hit> >(fHitModuleLabel);
    auto spHandle = evt.getValidHandle< std::vector<recob::SpacePoint> >(fSpModuleLabel);

    // also get the associated wires and raw digits;
    // we assume they have been created by the same module as the hits
    //art::FindOneP<raw::RawDigit> channelHitRawDigits(hitsHandle, evt, fHitModuleLabel);
    art::FindOneP<recob::Wire>   channelHitWires    (hitsHandle, evt, fHitModuleLabel);

    // this object contains the hit collection
    // and its associations to wires and raw digits:
    recob::HitCollectionCreator hcol(evt,
       channelHitWires.isValid(), // doWireAssns
     //channelHitRawDigits.isValid() // doRawDigitAssns
       false // do not save raw digits
      );

    // here is the copy of associations to hits, based on original hit assns
    auto assns = std::make_unique<art::Assns<recob::Hit, recob::SpacePoint>>(); 

    // all hits in the collection
    std::vector< art::Ptr<recob::Hit> > eventHits;
    art::fill_ptr_vector(eventHits, hitsHandle);

    art::FindManyP< recob::SpacePoint > spFromHit(hitsHandle, evt, fSpModuleLabel);
    art::FindManyP< recob::Hit > hitsFromSp(spHandle, evt, fSpModuleLabel);

    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
    // map induction spacepoints to TPC by collection hits
    std::unordered_map< size_t, size_t > spToTPC;
    for (size_t i = 0; i < spHandle->size(); ++i)
    {
        auto hits = hitsFromSp.at(i);
        size_t tpc = geo::WireID::InvalidID;
        for (const auto & h : hits) // find Collection hit, assume one is enough
        {
            if (h->SignalType() == geo::kCollection) { tpc = h->WireID().TPC; break; }
        }
        if (tpc == geo::WireID::InvalidID)
        {
	  //mf::LogWarning("DisambigFromSpacePoints") << "No collection hit for this spacepoint.";
            continue;
        }
        for (const auto & h : hits) // set mapping for Induction hits
        {
            if (h->SignalType() == geo::kInduction) { spToTPC[i] = tpc; }
        }
    }

    cryo_tpc_plane_keymap indHits;                        // induction hits resolved with spacepoints
    std::vector<size_t> unassignedHits;                   // hits to resolve by neighoring assignments
    std::unordered_map< size_t, geo::WireID > hitToWire;                 // final hit-wire assignments
    std::unordered_map< size_t, std::vector<geo::WireID> > hitToNWires;  // final hit-many-wires assignments

    hitToWire.reserve(eventHits.size());

    int n = runOnSpacePoints(eventHits, spFromHit, spToTPC, hitToWire, indHits, unassignedHits);
    mf::LogInfo("DisambigFromSpacePoints") << n << " hits undisambiguated by space points.";

    if (fUseNeighbors)
    {
        n = resolveUnassigned(detProp, hitToWire, eventHits, indHits, unassignedHits, fNumNeighbors);
        mf::LogInfo("DisambigFromSpacePoints") << n << " hits undisambiguated by neighborhood.";
    }

    if (fMoveLeftovers == "repeat")     { assignEveryAllowedWire(hitToNWires, eventHits, unassignedHits);      }
    else if (fMoveLeftovers == "first") { assignFirstAllowedWire(hitToWire, eventHits, unassignedHits);        }
    else                { mf::LogInfo("DisambigFromSpacePoints") << "Remaining undisambiguated hits dropped."; }

    auto const hitPtrMaker = art::PtrMaker<recob::Hit>(evt);

    for (auto const & hw : hitToWire)
    {
        size_t key = hw.first;
        geo::WireID wid = hw.second;

        recob::HitCreator new_hit(*(eventHits[key]), wid);

        hcol.emplace_back(new_hit.move(), channelHitWires.at(key));//, channelHitRawDigits.at(key));

        auto hitPtr = hitPtrMaker(hcol.size() - 1);
        auto sps = spFromHit.at(eventHits[key].key());
        for (auto const & spPtr : sps)
        {
            assns->addSingle(hitPtr, spPtr);
        }
    }

    for (auto const & hws : hitToNWires)
    {
        size_t key = hws.first;
        for (auto const & wid : hws.second)
        {
            recob::HitCreator new_hit(*(eventHits[key]), wid);

            hcol.emplace_back(new_hit.move(), channelHitWires.at(key));//, channelHitRawDigits.at(key));

            auto hitPtr = hitPtrMaker(hcol.size() - 1);
            auto sps = spFromHit.at(eventHits[key].key());
            for (auto const & spPtr : sps)
            {
                assns->addSingle(hitPtr, spPtr);
            }
        }
    }

    if (fMonitoringPlots && fTree) { fTree->Fill(); } // save statistics if MonitoringPlots was set to true

    // put the hit collection and associations into the event
    hcol.put_into(evt);
    evt.put(std::move(assns));
  }

  int DisambigFromSpacePoints::runOnSpacePoints(
    const std::vector< art::Ptr<recob::Hit> > & eventHits,
    const art::FindManyP< recob::SpacePoint > & spFromHit,
    const std::unordered_map< size_t, size_t > & spToTPC,
    std::unordered_map< size_t, geo::WireID > & assignments,
    cryo_tpc_plane_keymap & indHits,
    std::vector<size_t> & unassigned
    )
  {
    fNHits[0] = 0; fNHits[1] = 0; fNHits[2] = 0;
    fNMissedBySpacePoints[0] = 0;
    fNMissedBySpacePoints[1] = 0;

    for (size_t i = 0; i < eventHits.size(); ++i)
    {
        const art::Ptr<recob::Hit> & hit = eventHits[i];
        std::vector<geo::WireID> cwids = fGeom->ChannelToWire(hit->Channel());
        if (cwids.empty()) { mf::LogWarning("DisambigFromSpacePoints") << "No wires for this channel???"; continue; }
        if (hit->SignalType() == geo::kCollection)
        {
            assignments[hit.key()] = cwids.front();
            fNHits[2]++; // count collection hit
        }
        else
        {
            geo::WireID id = hit->WireID();
            size_t cryo = id.Cryostat, plane = id.Plane;
            fNHits[plane]++; // count induction hit

            if (spFromHit.at(hit.key()).size() == 0)
            {
                unassigned.push_back(hit.key());
                fNMissedBySpacePoints[plane]++; //count unresolved hit
            }
            else
            {
                std::unordered_map< size_t, geo::WireID > tpcBestWire;
                std::unordered_map< size_t, size_t > tpcScore;

                for (const auto & sp : spFromHit.at(hit.key()))
                {
                    auto search = spToTPC.find(sp.key());
                    if (search == spToTPC.end()) { continue; }
                    size_t spTpc = search->second;

                    const float max_dw = 1.; // max dist to wire [wire pitch]
                    for (size_t w = 0; w < cwids.size(); ++w)
                    {
                        if (cwids[w].TPC != spTpc) { continue; } // not that side of APA

                        float sp_wire = fGeom->WireCoordinate(sp->position(), id);
                        float dw = std::fabs(sp_wire - cwids[w].Wire);
                        if (dw < max_dw)
                        {
                            tpcBestWire[spTpc] = cwids[w];
                            tpcScore[spTpc]++;
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
                    assignments[hit.key()] = bestId;
                }
                else
                {
		  //mf::LogWarning("DisambigFromSpacePoints") << "Did not find matching wire (plane:" << plane << ").";
                    unassigned.push_back(hit.key());
                    fNMissedBySpacePoints[plane]++; //count unresolved hit
                }
            }
        }
    }

    return fNMissedBySpacePoints[0] + fNMissedBySpacePoints[1];
  }

  int DisambigFromSpacePoints::resolveUnassigned(
    detinfo::DetectorPropertiesData const& detProp,
    std::unordered_map< size_t, geo::WireID > & assignments,
    const std::vector< art::Ptr<recob::Hit> > & eventHits,
    cryo_tpc_plane_keymap & allIndHits,
    std::vector<size_t> & unassigned,
    size_t nNeighbors
    )
  {
    fNMissedByNeighbors[0] = 0;
    fNMissedByNeighbors[1] = 0;

    std::unordered_map< size_t, geo::WireID > result;

    for (const size_t key : unassigned)
    {
        const auto & hit = eventHits[key];
        geo::WireID id = hit->WireID();
        size_t cryo = id.Cryostat, plane = id.Plane;
        float hitDrift = hit->PeakTime();

        std::vector<geo::WireID> cwids = fGeom->ChannelToWire(hit->Channel());
        if (cwids.empty()) { mf::LogWarning("DisambigFromSpacePoints") << "No wires for this channel???"; continue; }

        const float dwMax = fMaxDistance / fGeom->TPC().Plane(plane).WirePitch(); // max distance in wires to look for neighbors
        const float ddMax = dwMax * fGeom->TPC().Plane(plane).WirePitch() / std::fabs(detProp.GetXTicksCoefficient(0, 0));

        float bestScore = 0;
        geo::WireID bestId;
        for (size_t w = 0; w < cwids.size(); ++w)
        {
            const size_t tpc = cwids[w].TPC;
            const size_t hitWire = cwids[w].Wire;

            bool allowed = true;
            for (auto t : fExcludeTPCs) { if (t == tpc) { allowed = false; break; } }
            if (!allowed) { continue; }

            geo::PlaneID const planeID(cryo, tpc, plane);
            const float wirePitch = fGeom->Plane(planeID).WirePitch();
            const float driftPitch = std::fabs(detProp.GetXTicksCoefficient(tpc, cryo));

            float maxDValue = fMaxDistance*fMaxDistance;
            std::vector<float> distBuff(nNeighbors, maxDValue); // distance to n closest hits
            const auto & keys = allIndHits[cryo][tpc][plane];
            for (const size_t keyInd : keys)
            {
                const auto & hitInd = eventHits[keyInd];

                auto search = assignments.find(keyInd);        // find resolved wire id
                if (search == assignments.end())
                {
                    mf::LogWarning("DisambigFromSpacePoints") << "Did not find resolved wire id.";
                    continue;
                }

                //float dWire = std::abs(hitWire - search->second.Wire);
                float dWire = std::abs(float(hitWire) - float(search->second.Wire));
                float dDrift = std::fabs(hitDrift - hitInd->PeakTime());

                if ((dWire > dwMax) || (dDrift > ddMax)) { continue; }

                dWire *= wirePitch;
                dDrift *= driftPitch;
                float dist2 = dWire * dWire + dDrift * dDrift;

                float maxd2 = 0;
                size_t maxIdx = 0;
                for (size_t i = 0; i < distBuff.size(); ++i )
                {
                    if (distBuff[i] > maxd2) { maxd2 = distBuff[i]; maxIdx = i; }
                }
                if (dist2 < maxd2) { distBuff[maxIdx] = dist2; }
            }
            float score = 0;
            size_t nhits = 0;
            for (size_t i = 0; i < distBuff.size(); ++i )
            {
                if (distBuff[i] < maxDValue) { score += 1./(1. + distBuff[i]); ++nhits; }
            }
            if (nhits > 0)
            {
                if (score > bestScore)
                {
                    bestScore = score;
                    bestId = cwids[w];
                }
            }
        }
        if (bestId.isValid) { result[key] = bestId; }
        else { fNMissedByNeighbors[plane]++; }
    }

    // remove from list of unassigned hits
    for (const auto & entry : result)
    {
        unassigned.erase(std::remove(unassigned.begin(), unassigned.end(), entry.first), unassigned.end());
    }

    //assignments.merge(result); // use this with C++ 17
    assignments.insert(result.begin(), result.end());
    return fNMissedByNeighbors[0] + fNMissedByNeighbors[1];
  }

  void DisambigFromSpacePoints::assignFirstAllowedWire(
    std::unordered_map< size_t, geo::WireID > & assignments,
    const std::vector< art::Ptr<recob::Hit> > & eventHits,
    const std::vector<size_t> & unassigned
    ) const
  {
    for (const size_t key : unassigned)
    {
        const auto & hit = eventHits[key];
        std::vector<geo::WireID> cwids = fGeom->ChannelToWire(hit->Channel());
        if (cwids.empty()) { mf::LogWarning("DisambigFromSpacePoints") << "No wires for this channel???"; continue; }

        geo::WireID bestId;
        for (size_t w = 0; w < cwids.size(); ++w)
        {
            const size_t tpc = cwids[w].TPC;

            bool allowed = true;
            for (auto t : fExcludeTPCs) { if (t == tpc) { allowed = false; break; } }
            if (allowed) { bestId = cwids[w]; break; }
        }
        if (bestId.isValid) { assignments[key] = bestId; }
        else { mf::LogWarning("DisambigFromSpacePoints") << "None of wires is allowed for this hit???"; }
    }
  }

  void DisambigFromSpacePoints::assignEveryAllowedWire(
    std::unordered_map< size_t, std::vector<geo::WireID> > & assignments,
    const std::vector< art::Ptr<recob::Hit> > & eventHits,
    const std::vector<size_t> & unassigned
    ) const
  {
    for (const size_t key : unassigned)
    {
        const auto & hit = eventHits[key];
        std::vector<geo::WireID> cwids = fGeom->ChannelToWire(hit->Channel());
        if (cwids.empty()) { continue; } // add warning here

        for (size_t w = 0; w < cwids.size(); ++w)
        {
            const size_t tpc = cwids[w].TPC;

            bool allowed = true;
            for (auto t : fExcludeTPCs) { if (t == tpc) { allowed = false; break; } }
            if (allowed) { assignments[key].push_back(cwids[w]); }
        }
    }
  }

  DEFINE_ART_MODULE(DisambigFromSpacePoints)

} // end of dune namespace
#endif 
