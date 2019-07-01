#include <vector>
#include <iostream>
#include <ctime>

#include "dune/CVN/func/GCNFeatureUtils.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"

cvn::GCNFeatureUtils::GCNFeatureUtils() {

}

cvn::GCNFeatureUtils::~GCNFeatureUtils(){

}

// Get the number of neighbours within rangeCut cm of this space point
const unsigned int cvn::GCNFeatureUtils::GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, 
                                                                 const float rangeCut, const std::string &spLabel) const{

  return GetAllNeighbours(evt,rangeCut,spLabel).at(sp.ID());

}

// Get the number of neighbours within rangeCut cm of this space point
const unsigned int cvn::GCNFeatureUtils::GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, 
                                                                 const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{

  return GetAllNeighbours(evt,rangeCut,sps).at(sp.ID());

}

// Get a map of the number of neighbours for each space point ID. Using this function is much less wasteful
// than repeated calls to the above function
const std::map<int,unsigned int> cvn::GCNFeatureUtils::GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::string &spLabel) const{

  // Get the space points from the event and make the map
  art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
  std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
  if(evt.getByLabel(spLabel,spacePointHandle)){
    art::fill_ptr_vector(allSpacePoints, spacePointHandle);
  }

  return GetAllNeighbours(evt,rangeCut,allSpacePoints);

}

const std::map<int,unsigned int> cvn::GCNFeatureUtils::GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{

  std::map<int,unsigned int> neighbourMap;

  for(const art::Ptr<recob::SpacePoint> sp0 : sps){
    // We want an entry even if it ends up being zero
    neighbourMap[sp0->ID()] = 0;

    for(const art::Ptr<recob::SpacePoint> sp1 : sps){

      if(sp0->ID() == sp1->ID()) continue;

      // For some reason we have to use arrays
      const double *p0 = sp0->XYZ();
      const double *p1 = sp1->XYZ();

      // Get the distance between the points
      const float dx = p1[0] - p0[0];
      const float dy = p1[1] - p0[1];
      const float dz = p1[2] - p0[2];
      const float dist = sqrt(dx*dx + dy*dy + dz*dz);

      if(dist < rangeCut){
        ++neighbourMap[sp0->ID()];
      }
    }
  }
  return neighbourMap;
} // function GetAllNeighbours

// Use the association between space points and hits to return a charge
const std::map<unsigned int, float> cvn::GCNFeatureUtils::GetSpacePointChargeMap(
  art::Event const &evt, const std::string &spLabel) const {

  std::map<unsigned int, float> ret;

  // Get the hits associated to the space points
  auto allSP = evt.getValidHandle<std::vector<recob::SpacePoint>>(spLabel);
  const art::FindManyP<recob::Hit> sp2Hit(allSP, evt, spLabel);

  for (size_t itSP = 0; itSP < allSP->size(); ++itSP) {
    const recob::SpacePoint& sp = allSP->at(itSP);
    float charge = 0.0;
    const std::vector<art::Ptr<recob::Hit>> spHits = sp2Hit.at(itSP);
    for (auto const hit : spHits) {
      charge += hit->Integral();
    }
    ret[sp.ID()] = charge;
  }

  return ret;

} // function GetSpacePointChargeMap

const std::map<unsigned int, unsigned int> cvn::GCNFeatureUtils::GetTrueG4ID(
  art::Event const& evt, const std::string &spLabel) const {

  art::ServiceHandle<cheat::BackTrackerService> bt;
  std::map<unsigned int, unsigned int> ret;

  // Get the hits associated to the space points
  auto allSP = evt.getValidHandle<std::vector<recob::SpacePoint>>(spLabel);
  const art::FindManyP<recob::Hit> sp2Hit(allSP, evt, spLabel);

  for (size_t itSP = 0; itSP < allSP->size(); ++itSP) {
    const recob::SpacePoint& sp = allSP->at(itSP);
    // Use the backtracker to find the G4 IDs associated with these hits
    std::map<unsigned int, float> trueParticles;
    for (auto hit : sp2Hit.at(itSP)) {
      auto ides = bt->HitToTrackIDEs(hit);
      for (auto ide : ides) {
        unsigned int id = abs(ide.trackID);
        if (trueParticles.count(id)) trueParticles[id] += ide.energy;
        else trueParticles[id] = ide.energy;
      }
    }
    ret[sp.ID()] = std::max_element(trueParticles.begin(), trueParticles.end(), [](const std::pair<unsigned int, unsigned int> &lhs,
      const std::pair<unsigned int, unsigned int> &rhs) { return lhs.second < rhs.second; })->first;
  }
  return ret;
} // function GetTrueG4ID
