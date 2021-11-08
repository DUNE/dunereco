#include <vector>
#include <iostream>
#include <ctime>

#include "dune/CVN/func/GCNFeatureUtils.h"
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

#include "TVector3.h"

using std::string;
using std::vector;
using std::pair;
using std::map;
using std::set;

using art::Ptr;
using art::Event;
using art::ServiceHandle;

using cheat::BackTrackerService;
using cheat::ParticleInventoryService;

using simb::MCParticle;

using recob::Hit;
using recob::SpacePoint;

namespace cvn
{

  GCNFeatureUtils::GCNFeatureUtils() {}

  GCNFeatureUtils::~GCNFeatureUtils() {}

  int GCNFeatureUtils::GetTrackIDFromHit(detinfo::DetectorClocksData const& clockData,
                                         const recob::Hit& hit) const {

    // Get the backtracker, then find the TrackIDE with the largest energy
    // deposit and return it
    art::ServiceHandle<cheat::BackTrackerService> bt;
    auto ides = bt->HitToTrackIDEs(clockData, hit);
    if (ides.empty()) return std::numeric_limits<int>::max(); // Return invalid number if no true tracks
    int id = 0;
    float energy = -1;
    for (auto ide : ides) {
      if (ide.energy > energy) {
        energy = ide.energy;
        id = ide.trackID;
      }
    }
    return id;

  } // function GCNFeatureUtils::GetTrackIDFromHit

  pair<unsigned int, float> GCNFeatureUtils::GetClosestApproach(const SpacePoint sp,
    const MCParticle p) const {

    // Spacepoint 3D position
    geoalgo::Point_t spPos(sp.XYZ()[0], sp.XYZ()[1], sp.XYZ()[2]);

    // Loop over trajectory segments and find point of closest approach
    unsigned int id = std::numeric_limits<unsigned int>::max();
    float dist = std::numeric_limits<float>::max();
    simb::MCTrajectory traj = p.Trajectory();
    for (size_t it = 1; it < traj.size(); ++it) {
      geoalgo::Point_t p1(TVector3(traj.Position(it-1).Vect()));
      geoalgo::Point_t p2(TVector3(traj.Position(it).Vect()));
      geoalgo::LineSegment_t ls(p1, p2);
      float distTmp = fGeoAlgo.SqDist(spPos, ls);
      if (distTmp < dist) {
        dist = distTmp;
        id = it-1;
      }
    }
    return std::make_pair(id, dist);

  } // function GCNFeatureUtils::GetClosestApproach

  // Get the number of neighbours within rangeCut cm of this space point
  unsigned int GCNFeatureUtils::GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt,
                                                                   const float rangeCut, const std::string &spLabel) const{

    return GetAllNeighbours(evt,rangeCut,spLabel).at(sp.ID());

  }

  // Get the number of neighbours within rangeCut cm of this space point
  unsigned int GCNFeatureUtils::GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt,
                                                                   const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{

    return GetAllNeighbours(evt,rangeCut,sps).at(sp.ID());

  }

  // Get a map of the number of neighbours for each space point ID. Using this function is much less wasteful
  // than repeated calls to the above function
  std::map<int,unsigned int> GCNFeatureUtils::GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::string &spLabel) const{

    // Get the space points from the event and make the map
    vector<Ptr<SpacePoint>> spacePoints;
    auto spacePointHandle = evt.getHandle<vector<SpacePoint>>(spLabel);
    if (spacePointHandle) {
      art::fill_ptr_vector(spacePoints, spacePointHandle);
    }
    return GetAllNeighbours(evt, rangeCut, spacePoints);
  }

  std::map<int,unsigned int> GCNFeatureUtils::GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{

    map<int,unsigned int> neighbourMap;

    for(const Ptr<SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      neighbourMap[sp0->ID()] = 0;

      for(const Ptr<SpacePoint> sp1 : sps){

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
  std::vector<std::map<int,unsigned int>> GCNFeatureUtils::GetNeighboursForRadii(art::Event const &evt, const std::vector<float>& rangeCuts, const std::string &spLabel) const{
    // Get the space points from the event and make the map
    vector<Ptr<SpacePoint>> allSpacePoints;
    auto spacePointHandle = evt.getHandle<vector<SpacePoint>>(spLabel);
    if (spacePointHandle) {
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetNeighboursForRadii(evt,rangeCuts,allSpacePoints);
  }

  std::vector<std::map<int,unsigned int>> GCNFeatureUtils::GetNeighboursForRadii(art::Event const &evt,
    const std::vector<float>& rangeCuts, const std::map<unsigned int,art::Ptr<recob::SpacePoint>> &sps) const{
    std::vector<art::Ptr<recob::SpacePoint>> vec;
    for(auto m : sps){
      vec.push_back(m.second);
    }
    return GetNeighboursForRadii(evt, rangeCuts, vec);
  }

  std::vector<std::map<int,unsigned int>> GCNFeatureUtils::GetNeighboursForRadii(art::Event const &evt,
    const std::vector<float>& rangeCuts, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{
    std::vector<std::map<int,unsigned int>> result;
    // Initialise a map for each range cut value
    for(unsigned int m = 0; m < rangeCuts.size(); ++m){
      result.push_back(map<int,unsigned int>());
    }
    for(const Ptr<SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      for(auto &m : result){
        m[sp0->ID()] = 0;
      }
      for(const Ptr<SpacePoint> sp1 : sps){
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

  std::map<int,int> GCNFeatureUtils::GetNearestNeighbours(art::Event const &evt,
    const std::string &spLabel) const{
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    auto spacePointHandle = evt.getHandle<std::vector<recob::SpacePoint>>(spLabel);
    if (spacePointHandle) {
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetNearestNeighbours(evt,allSpacePoints);
  }

  std::map<int,int> GCNFeatureUtils::GetNearestNeighbours(art::Event const &evt,
    const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{
    std::map<int,int> closestID;
    std::map<int,float> closestDistance;
    // Now loop over all of the space points and simultaneously fill our maps
    // Get the space points from the event and make the map
    for(const Ptr<SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      closestID[sp0->ID()] = 0;
      closestDistance[sp0->ID()] = 99999;
      for(const Ptr<SpacePoint> sp1 : sps){
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

  // Get the two nearest neighbours to use for calcuation of angles between them and the node in question
  std::map<int,std::pair<int,int>> GCNFeatureUtils::GetTwoNearestNeighbours(art::Event const &evt,
    const std::string &spLabel) const{
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    auto spacePointHandle = evt.getHandle<std::vector<recob::SpacePoint>>(spLabel);
    if (spacePointHandle) {
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetTwoNearestNeighbours(evt,allSpacePoints);
  }

  std::map<int,std::pair<int,int>> GCNFeatureUtils::GetTwoNearestNeighbours(art::Event const &evt,
    const std::map<unsigned int, art::Ptr<recob::SpacePoint>> &sps) const{

    vector<Ptr<SpacePoint>> vec;
    for(auto m : sps){
      vec.push_back(m.second);
    }
    return GetTwoNearestNeighbours(evt,vec);
  }

  std::map<int,std::pair<int,int>> GCNFeatureUtils::GetTwoNearestNeighbours(art::Event const &evt,
    const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{
    std::map<int,int> closestID;
    std::map<int,int> secondID;
    // Now loop over all of the space points and simultaneously fill our maps
    // Get the space points from the event and make the map
    for(const Ptr<SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      int thisSP = sp0->ID();
      int closest = -1;
      int second = -1;
      float closestDist = 99999;
      float secondDist = 99999;
      for(const Ptr<SpacePoint> sp1 : sps){
        if(thisSP == sp1->ID()) continue;
        // For some reason we have to use arrays
        const double *p0 = sp0->XYZ();
        const double *p1 = sp1->XYZ();
        // Get the distance between the points
        const float dx = p1[0] - p0[0];
        const float dy = p1[1] - p0[1];
        const float dz = p1[2] - p0[2];
        const float dist = sqrt(dx*dx + dy*dy + dz*dz);
        if(dist < closestDist){
          secondDist = closestDist;
          closestDist = dist;
          second = closest;
          closest = sp1->ID();
        }
        else if(dist < secondDist){
          secondDist = dist;
          second = sp1->ID();
        }
      }
      closestID.insert(std::make_pair(thisSP,closest));
      secondID.insert(std::make_pair(thisSP,second));
    }
    map<int,pair<int,int>> finalMap;
    for(unsigned int m = 0; m < closestID.size(); ++m){
      finalMap[m] = std::make_pair(closestID[m],secondID[m]);
    }
    return finalMap;
  }

  // Get the angle and the dot product between the vector from the base node to its neighbours
  void GCNFeatureUtils::GetAngleAndDotProduct(const SpacePoint &baseNode,
    const SpacePoint &n1, const SpacePoint &n2, float &dotProduct, float &angle) const{
    TVector3 basePos(baseNode.XYZ());
    TVector3 neighbour1Pos(n1.XYZ());
    TVector3 neighbour2Pos(n2.XYZ());
    TVector3 baseToNeighbour1 = neighbour1Pos - basePos;
    TVector3 baseToNeighbour2 = neighbour2Pos - basePos;
    dotProduct = baseToNeighbour1.Dot(baseToNeighbour2);
    angle = baseToNeighbour1.Angle(baseToNeighbour2);
    return;
  }

  // Use the association between space points and hits to return a charge
  std::map<unsigned int, float> GCNFeatureUtils::GetSpacePointChargeMap(
    std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const {

    map<unsigned int, float> ret;

    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      float charge = 0.0;
      for (Ptr<Hit> hit : sp2Hit[spIdx]) {
        charge += hit->Integral();
      }
      ret[spacePoints[spIdx]->ID()] = charge;
    }

    return ret;

  } // function GetSpacePointChargeMap

  std::map<unsigned int, float> GCNFeatureUtils::GetSpacePointChargeMap(
    art::Event const &evt, const std::string &spLabel) const {

    vector<Ptr<SpacePoint>> spacePoints;
    auto spacePointHandle = evt.getHandle<vector<SpacePoint>>(spLabel);
    if (!spacePointHandle) {
      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << spLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<Hit> fmp(spacePointHandle, evt, spLabel);
    vector<vector<Ptr<Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint

    return GetSpacePointChargeMap(spacePoints, sp2Hit);

  } // function GetSpacePointChargeMap

  std::map<unsigned int, float> GCNFeatureUtils::GetSpacePointMeanHitRMSMap(
    art::Event const &evt, const std::string &spLabel) const {

    map<unsigned int, float> ret;

    // Get the hits associated to the space points
    auto allSP = evt.getValidHandle<vector<SpacePoint>>(spLabel);
    const art::FindManyP<Hit> sp2Hit(allSP, evt, spLabel);

    for (size_t itSP = 0; itSP < allSP->size(); ++itSP) {
      const SpacePoint& sp = allSP->at(itSP);
      float chargeRMS = 0.0;
      const vector<Ptr<Hit>> spHits = sp2Hit.at(itSP);
      unsigned int nHits = spHits.size();
      for (auto const hit : spHits) {
        chargeRMS += hit->RMS();
      }
      ret[sp.ID()] = chargeRMS / static_cast<float>(nHits);
    }

    return ret;

  } // function GetSpacePointMeanHitRMSMap

  std::map<unsigned int, int> GCNFeatureUtils::GetTrueG4ID(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const {

    ServiceHandle<BackTrackerService> bt;
    map<unsigned int, int> ret;

    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      // Use the backtracker to find the G4 IDs associated with these hits
      std::map<unsigned int, float> trueParticles;
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
        auto ides = bt->HitToTrackIDEs(clockData, hit);
        for (auto ide : ides) {
          int id = ide.trackID;
          if (trueParticles.count(id)) trueParticles[id] += ide.energy;
          else trueParticles[id] = ide.energy;
        }
      }
      ret[spacePoints[spIdx]->ID()] = std::max_element(trueParticles.begin(), trueParticles.end(), [](const pair<unsigned int, float> &lhs,
        const pair<unsigned int, float> &rhs) { return lhs.second < rhs.second; })->first;
    }
    return ret;

  } // function GetTrueG4ID

  std::map<unsigned int, int> GCNFeatureUtils::GetTrueG4ID(
    detinfo::DetectorClocksData const& clockData,
    art::Event const &evt, const std::string &spLabel) const {

    vector<Ptr<SpacePoint>> spacePoints;
    auto spacePointHandle = evt.getHandle<vector<SpacePoint>>(spLabel);
    if (!spacePointHandle) {

      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << spLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<Hit> fmp(spacePointHandle, evt, spLabel);
    vector<vector<Ptr<Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint
    return GetTrueG4ID(clockData, spacePoints, sp2Hit);

  } // function GetTrueG4ID

  std::map<unsigned int, int> GCNFeatureUtils::GetTrueG4IDFromHits(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const {

    ServiceHandle<BackTrackerService> bt;
    map<unsigned int, int> ret;

    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      // Use the backtracker to find the G4 IDs associated with these hits
      std::map<unsigned int, unsigned int> trueParticleHits;
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
        auto ides = bt->HitToTrackIDEs(clockData, hit);
        for (auto ide : ides) {
          int id = ide.trackID;
          if (trueParticleHits.count(id)) trueParticleHits[id] += 1;
          else trueParticleHits[id] = 1;
        }
      }
      ret[spacePoints[spIdx]->ID()] = std::max_element(trueParticleHits.begin(), trueParticleHits.end(), [](const pair<unsigned int, unsigned> &lhs,
        const pair<unsigned int, unsigned int> &rhs) { return lhs.second < rhs.second; })->first;
    }
    return ret;

  } // function GetTrueG4IDFromHits

  std::map<unsigned int, int> GCNFeatureUtils::GetTrueG4IDFromHits(
    detinfo::DetectorClocksData const& clockData,
    art::Event const &evt, const std::string &spLabel) const {

    vector<Ptr<SpacePoint>> spacePoints;
    auto spacePointHandle = evt.getHandle<vector<SpacePoint>>(spLabel);
    if (!spacePointHandle) {

      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << spLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<Hit> fmp(spacePointHandle, evt, spLabel);
    vector<vector<Ptr<Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint
    return GetTrueG4IDFromHits(clockData, spacePoints, sp2Hit);

  } // function GetTrueG4IDFromHits

  std::map<unsigned int, int> GCNFeatureUtils::GetTruePDG(
    detinfo::DetectorClocksData const& clockData,
    art::Event const& evt, const std::string &spLabel, bool useAbsoluteTrackID, bool useHits) const {

    std::map<unsigned int, int> idMap;
    if(useHits) idMap  = GetTrueG4IDFromHits(clockData, evt,spLabel);
    else idMap = GetTrueG4ID(clockData, evt,spLabel);

    map<unsigned int,int> pdgMap;

    ServiceHandle<ParticleInventoryService> pi;

    // Now we need to get the true pdg code for each GEANT track ID in the map
    for(const pair<unsigned int, int> m : idMap){
      int pdg = 0;
      if(m.second == 0) std::cout << "Getting particle with ID " << m.second << " for space point " << m.first << std::endl;
      else{
        int trackID = m.second;
        if(useAbsoluteTrackID || trackID >= 0){
          trackID = abs(trackID);
          pdg = pi->TrackIdToParticle_P(trackID)->PdgCode();
        }
        else pdg = 11; // Dummy value to flag EM activity
      }
      pdgMap.insert(std::make_pair(m.first,pdg));
    }

    return pdgMap;
  } // function GetTrueG4ID

  std::map<unsigned int, std::vector<float>> GCNFeatureUtils::Get2DFeatures(
    std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const {

    // For each spacepoint, we want a size-nine float vector containing wire,
    // time and charge for each of the associated hits. These features assume
    // each spacepoint has a maximum of one hit associated from each plane,
    // and will throw an exception if it finds a spacepoint derived from
    // more than one hit on the same plane. Be warned!

    map<unsigned int, vector<float>> ret;
    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      vector<float> feat(9, 0.);
      // Loop over hits
      for (Ptr<Hit> hit : sp2Hit[spIdx]) {
        int offset = 3 * hit->WireID().Plane;
        // Throw an error if this plane's info has already been filled
        if (feat[offset+2] != 0) {
          std::ostringstream err;
          err << "2D features for plane " << hit->WireID().Plane << " have already been "
          << "filled.";
          throw std::runtime_error(err.str());
        }
        feat[offset] = hit->WireID().Wire;
        feat[offset+1] = hit->PeakTime();
        feat[offset+2] = hit->Integral();
      } // for hit
      ret[spacePoints[spIdx]->ID()] = feat;
    } // for spacepoint
    return ret;

  } // function GCNFeatureUtils::Get2DFeatures

  // Convert a pixel map into three 2D GCNGraph objects
  vector<GCNGraph> GCNFeatureUtils::ExtractGraphsFromPixelMap(const PixelMap &pm, const float chargeThreshold) const{

    // Each pixel map has three vectors of length (nWires*nTDCs)
    // Each value is the hit charge, and we will make GCNGraph for each view
    vector<vector<float>> allViews;
    allViews.push_back(pm.fPEX);
    allViews.push_back(pm.fPEY);
    allViews.push_back(pm.fPEZ);

    const unsigned int nWires = pm.fNWire;
    const unsigned int nTDCs = pm.fNTdc;

    vector<GCNGraph> outputGraphs;

    for(unsigned int v = 0; v < allViews.size(); ++v){

      GCNGraph newGraph;

      for(unsigned int w = 0; w < nWires; ++w){

        for(unsigned int t = 0; t < nTDCs; ++t){

          const unsigned int index = w*nTDCs + t;
          const float charge = allViews[v][index];

          // If the charge is very small then ignore this pixel
          if(charge < chargeThreshold) continue;

          // Put the positons into a vector
          vector<float> pos = {static_cast<float>(w),static_cast<float>(t)};
          // and the charge
          vector<float> features = {charge};

          GCNGraphNode newNode(pos,features);
          newGraph.AddNode(newNode);
        } // loop over TDCs
      } // loop over wires
      outputGraphs.push_back(newGraph);
    } // loop over views

    return outputGraphs;
  }

  std::map<unsigned int,unsigned int> GCNFeatureUtils::Get2DGraphNeighbourMap(const GCNGraph &g, const unsigned int npixel) const{

    map<unsigned int,unsigned int> neighbourMap;

    // Loop over the nodes and get the wire and tdc
    for(unsigned int n1 = 0; n1 < g.GetNumberOfNodes(); ++n1){

      neighbourMap[n1] = 0;
      const GCNGraphNode &node1 = g.GetNode(n1);
      const unsigned int w1 = static_cast<unsigned int>(node1.GetPosition()[0]);
      const unsigned int t1 = static_cast<unsigned int>(node1.GetPosition()[1]);

      // Loop over the nodes again to compare w,t coordinates
      for(unsigned int n2 = 0; n2 < g.GetNumberOfNodes(); ++n2){

        if(n1 == n2) continue;

        const GCNGraphNode &node2 = g.GetNode(n2);
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

  std::vector<float> GCNFeatureUtils::GetNodeGroundTruth(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit, float distCut,
    std::vector<std::vector<float>>* dirTruth) const{

    // Fetch cheat services
    ServiceHandle<BackTrackerService> bt;
    ServiceHandle<ParticleInventoryService> pi;

    vector<float> ret(spacePoints.size(), -1);
    if (dirTruth) {
      dirTruth->clear();
      dirTruth->resize(spacePoints.size());
    }
    vector<pair<size_t, float>> dist;
    vector<int> trueIDs(spacePoints.size(), -1);
    vector<int> closestTrajPoint(spacePoints.size(), -1);
    // vector<float> dist(spacePoints.size(), -1);

    // Debug counters!
    size_t nMismatched(0), nNoHit(0), nGood(0);//, nOutOfRange(0), nHitUsed(0), nGood(0);
    // Loop over each spacepoint
    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      // Check index and ID agree
      if (spIdx != (size_t)spacePoints[spIdx]->ID()) {
        throw art::Exception(art::errors::LogicError) << "Spacepoint index "
          << spIdx << " mismatched with spacepoint ID "
          << spacePoints[spIdx]->ID();
      }
      bool done = false;
      Ptr<SpacePoint> sp = spacePoints[spIdx];

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
      //   4) Is the spacepoint's hit already used by a different spacepoint?
      // If yes to all of these, the spacepoint is true. Otherwise, it's false.

      int trueParticleID = std::numeric_limits<int>::max();
      bool firstHit = true;
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
        int trueID = GetTrackIDFromHit(clockData, *hit);
        if (trueID == std::numeric_limits<int>::max()) {
          ++nNoHit;
          done = true;
          break;
        }
        // Check whether this hit's true particle matches other hits from this spacepoint
        if (firstHit) {
          trueParticleID = trueID;
          firstHit = false;
        }
        else if (trueParticleID != trueID) {
          // If true IDs don't match, quit
          ++nMismatched;
          done = true;
          break;
        }
      } // for hit

      if (done) continue; // if this isn't a valid spacepoint then stop here
      if (trueParticleID == std::numeric_limits<int>::max()) {
        ++nNoHit;
        continue;
      }

      // Get the true particle, and find the closest MC trajectory point to spacepoint
      MCParticle p = pi->TrackIdToParticle(abs(trueParticleID));
      pair<unsigned int, float> closest = GetClosestApproach(*sp, p);
      dist.push_back(std::make_pair(spIdx, closest.second));
      trueIDs[spIdx] = abs(trueParticleID);
      closestTrajPoint[spIdx] = closest.first;
      ret[spIdx] = closest.second;
      ++nGood;
      //if (closest.second < distCut) {
        // ++nGood;
        //ret[spIdx] = -1;
      //}
      //else {
      //  ++nOutOfRange;
      //}

    } // for spIdx

    // Sort spacepoints in ascending order of distance to true trajectory
    //std::sort(std::begin(dist), std::end(dist),
    //  [](const auto& a, const auto&b) { return a.second < b.second; });
    //set<size_t> usedHits;
    // Loop over all truth-matched spacepoints
    //for (auto d : dist) {
    //  // Check whether any of the hits associated with this spacepoint have already been used
    //  for (Ptr<Hit> hit : sp2Hit[d.first]) {
    //    if (usedHits.count(hit.key())) {
    //      ret[d.first] = false; // Remove this spacepoint from the list of good ones
    //      ++nHitUsed;
    //      break;
    //    }
    //  }
    //  if (!ret[d.first]) continue; // If this spacepoint is bad, move on to the next one
    //  // None of the hits have been used, so we'll let this one through but
    //  // add its hits to the list so no other spacepoints can use them
    //  for (Ptr<Hit> hit : sp2Hit[d.first]) {
    //    usedHits.insert(hit.key());
    //  }
    //  ++nGood;

    //} // for spacepoint

    // Loop over spacepoints to set directionality
    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      // If this was a true node & we're storing direction, store it!
      if (ret[spIdx] && dirTruth) {
        MCParticle p = pi->TrackIdToParticle(trueIDs[spIdx]);
        TVector3 dir = p.Trajectory().Momentum(closestTrajPoint[spIdx]).Vect().Unit();
        for (size_t it = 0; it < 3; ++it) {
          dirTruth->at(spIdx).push_back(dir[it]);
        }
      }

    } // for spIdx

    std::cout << "There were " << nGood << " valid spacepoints, " << nNoHit//<< nHitUsed
//      << " spacepoints with a hit already used, " << nNoHit
      << " spacepoints with no associated hits or noise hits, and " << nMismatched
      << " spacepoints produced from mismatched hits." << std::endl;// and " << nOutOfRange
//      << " spacepoints reconstructed too far away from the true particle trajectory."
//      << std::endl;

    return ret;

  } // function GCNFeatureUtils::GetNodeGroundTruth

  std::map<unsigned int, unsigned int> GCNFeatureUtils::GetParticleFlowMap(
    const std::set<unsigned int>& particles) const {

    ServiceHandle<ParticleInventoryService> pi;
    map<unsigned int, unsigned int> ret;
    for (int p : particles) {
      if (p == 0) continue; // No parent for the event primary
      ret[p] = pi->TrackIdToParticle(abs(p)).Mother();
    }
    return ret;

  } // function GCNFeatureUtils::GetParticleFlowMap

  vector<ptruth> GCNFeatureUtils::GetParticleTree(
    const cvn::GCNGraph* g) {

    ServiceHandle<ParticleInventoryService> pi;

    set<int> trackIDs;
    for (unsigned int i = 0; i < g->GetNumberOfNodes(); ++i) {
      trackIDs.insert(g->GetNode(i).GetGroundTruth()[0]);
    }

    set<int> allIDs = trackIDs; // Copy original so we can safely modify it

    vector<ptruth> ret;

    // Add invisible particles to hierarchy
    for (int id : trackIDs) {
      const MCParticle* p = pi->TrackIdToParticle_P(abs(id));
      while (p->Mother() != 0) {
        allIDs.insert(abs(p->Mother()));
        p = pi->TrackIdToParticle_P(abs(p->Mother()));
      }
    }

    for (int id : allIDs) {
      const MCParticle* p = pi->TrackIdToParticle_P(abs(id));

      ret.push_back(std::make_tuple(abs(id), p->PdgCode(), p->Mother(),
        p->P(), p->Vx(), p->Vy(), p->Vz(), p->EndX(), p->EndY(), p->EndZ(),
        p->Process(), p->EndProcess()));
    }

    return ret;

  } // function GCNFeatureUtils::GetParticleTree

} // namespace cvn
