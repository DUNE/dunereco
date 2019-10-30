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
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNGraphNode.h"
#include "dune/CVN/func/PixelMap.h"

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

// Sometimes we might want to know the number of neighbours within various radii
const std::vector<std::map<int,unsigned int>> cvn::GCNFeatureUtils::GetNeighbourForRadii(art::Event const &evt,
                                          const std::vector<float> rangeCuts, const std::string &spLabel) const{

  std::vector<std::map<int,unsigned int>> result;
  // Initialise a map for each range cut value
  for(unsigned int m = 0; m < rangeCuts.size(); ++m){
    result.push_back(std::map<int,unsigned int>());
  }

  // Now loop over all of the space points and simultaneously fill our maps
  // Get the space points from the event and make the map
  art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
  std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
  if(evt.getByLabel(spLabel,spacePointHandle)){
    art::fill_ptr_vector(allSpacePoints, spacePointHandle);
  }
  for(const art::Ptr<recob::SpacePoint> sp0 : allSpacePoints){
    // We want an entry even if it ends up being zero
    for(auto &m : result){
      m[sp0->ID()] = 0;
    }

    for(const art::Ptr<recob::SpacePoint> sp1 : allSpacePoints){

      if(sp0->ID() == sp1->ID()) continue;

      // For some reason we have to use arrays
      const double *p0 = sp0->XYZ();
      const double *p1 = sp1->XYZ();

      // Get the distance between the points
      const float dx = p1[0] - p0[0];
      const float dy = p1[1] - p0[1];
      const float dz = p1[2] - p0[2];
      const float dist = sqrt(dx*dx + dy*dy + dz*dz);

      // Fill the maps if we satify the criteria
      for(unsigned int r = 0; r < rangeCuts.size(); ++r){
        if(dist < rangeCuts[r]){
          result[r][sp0->ID()] = result[r][sp0->ID()] + 1;
        }
      }

    }
  }
  return result;

}

const std::map<int,int> cvn::GCNFeatureUtils::GetNearestNeighbours(art::Event const &evt, const std::string &spLabel) const{

  std::map<int,int> closestID;
  std::map<int,float> closestDistance;
  // Now loop over all of the space points and simultaneously fill our maps
  // Get the space points from the event and make the map
  art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
  std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
  if(evt.getByLabel(spLabel,spacePointHandle)){
    art::fill_ptr_vector(allSpacePoints, spacePointHandle);
  }

  for(const art::Ptr<recob::SpacePoint> sp0 : allSpacePoints){
    // We want an entry even if it ends up being zero
    closestID[sp0->ID()] = 0;
    closestDistance[sp0->ID()] = 99999;

    for(const art::Ptr<recob::SpacePoint> sp1 : allSpacePoints){

      if(sp0->ID() == sp1->ID()) continue;

      // For some reason we have to use arrays
      const double *p0 = sp0->XYZ();
      const double *p1 = sp1->XYZ();

      // Get the distance between the points
      const float dx = p1[0] - p0[0];
      const float dy = p1[1] - p0[1];
      const float dz = p1[2] - p0[2];
      const float dist = sqrt(dx*dx + dy*dy + dz*dz);

      if(dist < closestDistance[sp0->ID()]){
        closestDistance[sp0->ID()] = dist;
        closestID[sp0->ID()] = sp1->ID();
      }

    }
  }
  return closestID;

}

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

// Convert a pixel map into three 2D GCNGraph objects
std::vector<cvn::GCNGraph> cvn::GCNFeatureUtils::ExtractGraphsFromPixelMap(const cvn::PixelMap &pm, const float chargeThreshold) const{

  // Each pixel map has three vectors of length (nWires*nTDCs)
  // Each value is the hit charge, and we will make GCNGraph for each view
  std::vector<std::vector<float>> allViews;
  allViews.push_back(pm.fPEX); 
  allViews.push_back(pm.fPEY);
  allViews.push_back(pm.fPEZ);

  const unsigned int nWires = pm.fNWire;
  const unsigned int nTDCs = pm.fNTdc;

  std::vector<cvn::GCNGraph> outputGraphs;

  for(unsigned int v = 0; v < allViews.size(); ++v){

    cvn::GCNGraph newGraph;

    for(unsigned int w = 0; w < nWires; ++w){

      for(unsigned int t = 0; t < nTDCs; ++t){

        const unsigned int index = w*nTDCs + t;
        const float charge = allViews[v][index];

        // If the charge is very small then ignore this pixel
        if(charge < chargeThreshold) continue;        

        // Put the positons into a vector
        std::vector<float> pos = {static_cast<float>(w),static_cast<float>(t)};
        // and the charge
        std::vector<float> features = {charge};

        cvn::GCNGraphNode newNode(pos,features);
        newGraph.AddNode(newNode);
      } // loop over TDCs
    } // loop over wires
    outputGraphs.push_back(newGraph);
  } // loop over views

  return outputGraphs;
}

const std::map<unsigned int,unsigned int> cvn::GCNFeatureUtils::Get2DGraphNeighbourMap(const cvn::GCNGraph &g, const unsigned int npixel) const{

  std::map<unsigned int,unsigned int> neighbourMap;

  // Loop over the nodes and get the wire and tdc
  for(unsigned int n1 = 0; n1 < g.GetNumberOfNodes(); ++n1){

    neighbourMap[n1] = 0;
    const cvn::GCNGraphNode &node1 = g.GetNode(n1);
    const unsigned int w1 = static_cast<unsigned int>(node1.GetPosition()[0]);
    const unsigned int t1 = static_cast<unsigned int>(node1.GetPosition()[1]);

    // Loop over the nodes again to compare w,t coordinates
    for(unsigned int n2 = 0; n2 < g.GetNumberOfNodes(); ++n2){

      if(n1 == n2) continue;

      const cvn::GCNGraphNode &node2 = g.GetNode(n2);
      const unsigned int w2 = static_cast<unsigned int>(node2.GetPosition()[0]);
      const unsigned int t2 = static_cast<unsigned int>(node2.GetPosition()[1]);

      // Check if the box around the node for neighbours
      // In this example npixel = 2.
      //
      // |---|---|---|---|---|---|---|            
      // |   |   |   |   |   |   |   |
      // |---|---|---|---|---|---|---|            
      // |   | x | x | x | x | x |   |
      // |---|---|---|---|---|---|---|            
      // |   | x | x | x | x | x |   |
      // |---|---|---|---|---|---|---|            
      // |   | x | x |n1 | x | x |   | t
      // |---|---|---|---|---|---|---|            
      // |   | x | x | x | x | x |   |
      // |---|---|---|---|---|---|---|            
      // |   | x | x | x | x | x |   |
      // |---|---|---|---|---|---|---|            
      // |   |   |   |   |   |   |   |
      // |---|---|---|---|---|---|---|            
      //               w

      if(w2 <= w1 + npixel && w2 >= w1 - npixel){
        if(t2 <= t1 + npixel && t2 >= t1 - npixel){
          neighbourMap[n1]++;
        }
      }
     
    }
  }

  return neighbourMap;

}

const std::vector<bool> cvn::GCNFeatureUtils::GetNodeGroundTruth(art::Event const &evt,
  const std::string &spLabel, float distCut) const{

  // Fetch cheat services
  art::ServiceHandle<cheat::BackTrackerService> bt;
  art::ServiceHandle<cheat::ParticleInventoryService> pi;

  // Get the hits associated to the space points
  auto allSP = evt.getValidHandle<std::vector<recob::SpacePoint>>(spLabel);
  const art::FindManyP<recob::Hit> sp2Hit(allSP, evt, spLabel);

  std::vector<bool> ret(allSP->size(), false);

  // Loop over each spacepoint
  for (size_t itSP = 0; itSP < allSP->size(); ++itSP) {
    bool done = false;
    const recob::SpacePoint& sp = allSP->at(itSP);

    // Loop over each hit associated with this spacepoint. At this point we
    // ask three questions to figure out if this should be considered a "true"
    // spacepoint:
    //   1) Do all the hits associated with this spacepoint come from a true
    //      particle (ie. are not noise)?
    //   2) If yes to 1), do all the hits associated with this spacepoint
    //      derive the majority of their energy from the same true particle?
    //   3) Does the reconstructed position of the spacepoint fall within
    //      some FHICL-configurable distance to the true particle's
    //      trajectory?
    // If yes to all three of these, the spacepoint is true. Otherwise, it's
    // false.

    int trueParticleID = -1;
    for (auto hit : sp2Hit.at(itSP)) {
      auto ides = bt->HitToTrackIDEs(hit);
      if (ides.empty()) {
        std::cout << "No true particles for this hit!" << std::endl;
        done = true;
        break;
      }
      int trueID = -1;
      float trueEnergy = -1;
      // Loop over IDEs and find the particle that contributed the most energy
      for (auto ide : ides) {
        if (ide.energy > trueEnergy) {
          trueID = abs(ide.trackID);
          trueEnergy = ide.energy;
        }
      } // for ide

      // Check whether this hit's true particle matches other hits from this spacepoint
      if (trueParticleID == -1) trueParticleID = trueID;
      else if (trueParticleID != trueID) {
        // If true IDs don't match, quit
        std::cout << "True particles do not match!" << std::endl;
        done = true;
        break;
      }

    } // for hit

    if (done) continue;

    // Now use the particle inventory service to get the true particle by id
    simb::MCParticle p = pi->TrackIdToParticle(trueParticleID);
    simb::MCTrajectory traj = p.Trajectory();
    // Now loop over the trajectory and find the smallest distance!
    TVector3 pos(sp.XYZ());
    for (size_t it = 0; it < traj.size(); ++it) {
      if (it > 0)
        std::cout << "Space between trajectory points is " << (traj.Position(it).Vect() - traj.Position(it-1).Vect()).Mag() << std::endl;
      TVector3 diff = traj.Position(it).Vect() - pos;
      float dist = diff.Mag();
      if (dist < distCut) {
        std::cout << "Found a trajectory point with distance " << dist << "!" << std::endl;
        ret[itSP] = true;
        break;
      }
    } // for trajectory point it

    if (!ret[itSP]) std::cout << "No trajectory point within range." << std::endl;

  } // for itSP

  return ret;

} // function GCNFeatureUtils::GetNodeGroundTruth



