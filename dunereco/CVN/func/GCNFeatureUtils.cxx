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

namespace cvn
{

  GCNFeatureUtils::GCNFeatureUtils() {}

  GCNFeatureUtils::~GCNFeatureUtils() {}

  int GCNFeatureUtils::GetTrackIDFromHit(const recob::Hit hit) const {

    // Get the backtracker, then find the TrackIDE with the largest energy
    // deposit and return it
    art::ServiceHandle<cheat::BackTrackerService> bt;
    auto ides = bt->HitToTrackIDEs(hit);
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

  size_t GCNFeatureUtils::GetClosestTrajectoryPoint(const recob::SpacePoint sp,
    const simb::MCParticle p) const {

    // Spacepoint 3D position
    float dist = std::numeric_limits<float>::max();
    TVector3 spPos = TVector3(sp.XYZ());
    size_t ret = std::numeric_limits<size_t>::max();

    // Use particle inventory service to get the trajectory of MC particle
    simb::MCTrajectory traj = p.Trajectory();
    for (size_t it = 0; it < traj.size(); ++it) {
      TVector3 trajPos = traj.Position(it).Vect();
      float distTmp = (spPos - trajPos).Mag();
      if (distTmp < dist) {
        dist = distTmp;
        ret = it;
      }
    } // for trajectory point it

    // Check we actually have a valid ID here
    if (ret == std::numeric_limits<size_t>::max()) {
      throw std::runtime_error("No valid trajectory point found!");
    }

    return ret;

  } // function GCNFeatureUtils::GetClosestTrajectoryPoint

  // Get the number of neighbours within rangeCut cm of this space point
  const unsigned int GCNFeatureUtils::GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, 
                                                                   const float rangeCut, const std::string &spLabel) const{

    return GetAllNeighbours(evt,rangeCut,spLabel).at(sp.ID());

  }

  // Get the number of neighbours within rangeCut cm of this space point
  const unsigned int GCNFeatureUtils::GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, 
                                                                   const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{

    return GetAllNeighbours(evt,rangeCut,sps).at(sp.ID());

  }

  // Get a map of the number of neighbours for each space point ID. Using this function is much less wasteful
  // than repeated calls to the above function
  const std::map<int,unsigned int> GCNFeatureUtils::GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::string &spLabel) const{

    // Get the space points from the event and make the map
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints;
    if(evt.getByLabel(spLabel,spacePointHandle)){
      art::fill_ptr_vector(spacePoints, spacePointHandle);
    }
    return GetAllNeighbours(evt, rangeCut, spacePoints);
  }

  const std::map<int,unsigned int> GCNFeatureUtils::GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{

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
  const std::vector<std::map<int,unsigned int>> GCNFeatureUtils::GetNeighboursForRadii(art::Event const &evt, const std::vector<float> rangeCuts, const std::string &spLabel) const{
    // Get the space points from the event and make the map
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    if(evt.getByLabel(spLabel,spacePointHandle)){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetNeighboursForRadii(evt,rangeCuts,allSpacePoints);
  }

  const std::vector<std::map<int,unsigned int>> GCNFeatureUtils::GetNeighboursForRadii(art::Event const &evt,
    const std::vector<float> rangeCuts, const std::map<unsigned int,art::Ptr<recob::SpacePoint>> &sps) const{
    std::vector<art::Ptr<recob::SpacePoint>> vec;
    for(auto m : sps){
      vec.push_back(m.second);
    }
    return GetNeighboursForRadii(evt,rangeCuts,vec);
  }

  const std::vector<std::map<int,unsigned int>> GCNFeatureUtils::GetNeighboursForRadii(art::Event const &evt,
    const std::vector<float> rangeCuts, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{
    std::vector<std::map<int,unsigned int>> result;
    // Initialise a map for each range cut value
    for(unsigned int m = 0; m < rangeCuts.size(); ++m){
      result.push_back(std::map<int,unsigned int>());
    }
    for(const art::Ptr<recob::SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      for(auto &m : result){
        m[sp0->ID()] = 0;
      }
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

  const std::map<int,int> GCNFeatureUtils::GetNearestNeighbours(art::Event const &evt,
    const std::string &spLabel) const{
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    if(evt.getByLabel(spLabel,spacePointHandle)){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetNearestNeighbours(evt,allSpacePoints);
  }

  const std::map<int,int> GCNFeatureUtils::GetNearestNeighbours(art::Event const &evt,
    const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{
    std::map<int,int> closestID;
    std::map<int,float> closestDistance;
    // Now loop over all of the space points and simultaneously fill our maps
    // Get the space points from the event and make the map
    for(const art::Ptr<recob::SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      closestID[sp0->ID()] = 0;
      closestDistance[sp0->ID()] = 99999;
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
        if(dist < closestDistance[sp0->ID()]){
          closestDistance[sp0->ID()] = dist;
          closestID[sp0->ID()] = sp1->ID();
        }
      }
    }
    return closestID;
  }

  // Get the two nearest neighbours to use for calcuation of angles between them and the node in question
  const std::map<int,std::pair<int,int>> GCNFeatureUtils::GetTwoNearestNeighbours(art::Event const &evt,
    const std::string &spLabel) const{
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    if(evt.getByLabel(spLabel,spacePointHandle)){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    return GetTwoNearestNeighbours(evt,allSpacePoints);
  }

  const std::map<int,std::pair<int,int>> GCNFeatureUtils::GetTwoNearestNeighbours(art::Event const &evt,
    const std::map<unsigned int, art::Ptr<recob::SpacePoint>> &sps) const{

    std::vector<art::Ptr<recob::SpacePoint>> vec;
    for(auto m : sps){
      vec.push_back(m.second);
    }
    return GetTwoNearestNeighbours(evt,vec);
  }

  const std::map<int,std::pair<int,int>> GCNFeatureUtils::GetTwoNearestNeighbours(art::Event const &evt,
    const std::vector<art::Ptr<recob::SpacePoint>> &sps) const{
    std::map<int,int> closestID;
    std::map<int,int> secondID;
    // Now loop over all of the space points and simultaneously fill our maps
    // Get the space points from the event and make the map
    for(const art::Ptr<recob::SpacePoint> sp0 : sps){
      // We want an entry even if it ends up being zero
      int thisSP = sp0->ID();
      int closest = -1;
      int second = -1;
      float closestDist = 99999;
      float secondDist = 99999;
      for(const art::Ptr<recob::SpacePoint> sp1 : sps){
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
    std::map<int,std::pair<int,int>> finalMap;
    for(unsigned int m = 0; m < closestID.size(); ++m){
      finalMap[m] = std::make_pair(closestID[m],secondID[m]);
    }
    return finalMap;
  }

  // Get the angle and the dot product between the vector from the base node to its neighbours
  void GCNFeatureUtils::GetAngleAndDotProduct(const recob::SpacePoint &baseNode,
    const recob::SpacePoint &n1, const recob::SpacePoint &n2, float &dotProduct, float &angle) const{
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
  const std::map<unsigned int, float> GCNFeatureUtils::GetSpacePointChargeMap(
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit) const {

    std::map<unsigned int, float> ret;

    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      float charge = 0.0;
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
        charge += hit->Integral();
      }
      ret[spacePoints[spIdx]->ID()] = charge;
    }

    return ret;

  } // function GetSpacePointChargeMap

  const std::map<unsigned int, float> GCNFeatureUtils::GetSpacePointChargeMap(
    art::Event const &evt, const std::string &spLabel) const {

    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints;
    if (!evt.getByLabel(spLabel, spacePointHandle)) {
      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << spLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<recob::Hit> fmp(spacePointHandle, evt, spLabel);
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint

    return GetSpacePointChargeMap(spacePoints, sp2Hit);

  } // function GetSpacePointChargeMap

  const std::map<unsigned int, int> GCNFeatureUtils::GetTrueG4ID(
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit) const {

    art::ServiceHandle<cheat::BackTrackerService> bt;
    std::map<unsigned int, int> ret;

    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      // Use the backtracker to find the G4 IDs associated with these hits
      std::map<unsigned int, float> trueParticles;
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
        auto ides = bt->HitToTrackIDEs(hit);
        for (auto ide : ides) {
          int id = ide.trackID;
          if (trueParticles.count(id)) trueParticles[id] += ide.energy;
          else trueParticles[id] = ide.energy;
        }
      }
      ret[spacePoints[spIdx]->ID()] = std::max_element(trueParticles.begin(), trueParticles.end(), [](const std::pair<unsigned int, float> &lhs,
        const std::pair<unsigned int, float> &rhs) { return lhs.second < rhs.second; })->first;
    }
    return ret;

  } // function GetTrueG4ID

  const std::map<unsigned int, int> GCNFeatureUtils::GetTrueG4ID(
    art::Event const &evt, const std::string &spLabel) const {

    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints;
    if (!evt.getByLabel(spLabel, spacePointHandle)) {

      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << spLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<recob::Hit> fmp(spacePointHandle, evt, spLabel);
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint
    return GetTrueG4ID(spacePoints, sp2Hit);

  } // function GetTrueG4ID

  const std::map<unsigned int, int> cvn::GCNFeatureUtils::GetTruePDG(
    art::Event const& evt, const std::string &spLabel) const {
  
    const std::map<unsigned int, int> idMap = GetTrueG4ID(evt,spLabel);
    std::map<unsigned int,int> pdgMap;  
    
    art::ServiceHandle<cheat::ParticleInventoryService> pi;
  
    // Now we need to get the true pdg code for each GEANT track ID in the map
    for(const std::pair<unsigned int, int> m : idMap){
      int pdg = 0;
      if(m.second == 0) std::cout << "Getting particle with ID " << m.second << " for space point " << m.first << std::endl;
      else pdg = pi->TrackIdToParticle_P(abs(m.second))->PdgCode();
      pdgMap.insert(std::make_pair(m.first,pdg));
    }
  
    return pdgMap;
  } // function GetTrueG4ID

  const std::map<unsigned int, std::vector<float>> GCNFeatureUtils::Get2DFeatures(
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit) const {

    // For each spacepoint, we want a size-nine float vector containing wire,
    // time and charge for each of the associated hits. These features assume
    // each spacepoint has a maximum of one hit associated from each plane,
    // and will throw an exception if it finds a spacepoint derived from
    // more than one hit on the same plane. Be warned!

    std::map<unsigned int, std::vector<float>> ret;
    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      std::vector<float> feat(9, 0.);
      // Loop over hits
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
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
  std::vector<GCNGraph> GCNFeatureUtils::ExtractGraphsFromPixelMap(const PixelMap &pm, const float chargeThreshold) const{

    // Each pixel map has three vectors of length (nWires*nTDCs)
    // Each value is the hit charge, and we will make GCNGraph for each view
    std::vector<std::vector<float>> allViews;
    allViews.push_back(pm.fPEX); 
    allViews.push_back(pm.fPEY);
    allViews.push_back(pm.fPEZ);

    const unsigned int nWires = pm.fNWire;
    const unsigned int nTDCs = pm.fNTdc;

    std::vector<GCNGraph> outputGraphs;

    for(unsigned int v = 0; v < allViews.size(); ++v){

      GCNGraph newGraph;

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

          GCNGraphNode newNode(pos,features);
          newGraph.AddNode(newNode);
        } // loop over TDCs
      } // loop over wires
      outputGraphs.push_back(newGraph);
    } // loop over views

    return outputGraphs;
  }

  const std::map<unsigned int,unsigned int> GCNFeatureUtils::Get2DGraphNeighbourMap(const GCNGraph &g, const unsigned int npixel) const{

    std::map<unsigned int,unsigned int> neighbourMap;

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

  const std::vector<bool> GCNFeatureUtils::GetNodeGroundTruth(
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints,
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit, float distCut,
    std::vector<std::vector<float>>* dirTruth) const{

    // Fetch cheat services
    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi;

    std::vector<bool> ret(spacePoints.size(), false);
    if (dirTruth) {
      dirTruth->clear();
      dirTruth->resize(spacePoints.size());
    }

    // Debug counters!
    size_t nIDEMismatched(0), nMismatched(0), nNoHit(0), nOutOfRange(0), nGood(0);
    // Loop over each spacepoint
    for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
      bool done = false;
      art::Ptr<recob::SpacePoint> sp = spacePoints[spIdx];

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

      // Extra piece of code to try and match IDEs!
      bool firstHit = true;
      TVector3 spPos(sp->XYZ());
      TVector3 pos(0., 0., 0.);
      int trueID = std::numeric_limits<int>::max();
      for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
        float energy = -1;
        int trueParticleID = std::numeric_limits<int>::max();
        TVector3 hitPos(0., 0., 0.);
        for (const sim::IDE* ide : bt->HitToSimIDEs_Ps(hit)) {
          if (ide->energy > energy) {
            energy = ide->energy;
            trueParticleID = ide->trackID;
            hitPos.SetX(ide->x);
            hitPos.SetY(ide->y);
            hitPos.SetZ(ide->z);
          } // if higher-energy IDE
        } // for IDE
        if (firstHit) {
          firstHit = false;
          pos = hitPos;
          trueID = trueParticleID;
        }
        else {
          if (trueID != trueParticleID) {
            ++nMismatched;
            done = true;
            break;
          }
          float dist = (pos-hitPos).Mag();
          // Check IDE distance is small
          if (dist > 0.1) {
            ++nIDEMismatched;
            done = true;
            break;
          }
          // else {
          //   std::cout << "Here we have a distance of " << dist << "between IDEs!" << std::endl;
          // }
        }
      } // for hit

      if (trueID == std::numeric_limits<int>::max()) {
        ++nNoHit;
        done = true;
      }

      if (done) continue; // if it isn't a valid spacepoint then stop here

      // int trueParticleID = std::numeric_limits<int>::max();
      // firstHit = true;
      // for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
      //   int trueID = GetTrackIDFromHit(*hit);
      //   if (trueID == std::numeric_limits<int>::max()) {
      //     ++nNoHit;
      //     done = true;
      //     break;
      //   }
      //   // Check whether this hit's true particle matches other hits from this spacepoint
      //   if (firstHit) {
      //     trueParticleID = trueID;
      //     firstHit = false;
      //   }
      //   else if (trueParticleID != trueID) {
      //     // If true IDs don't match, quit
      //     ++nMismatched;
      //     done = true;
      //     break;
      //   }
      // } // for hit

      // if (done) continue; // if this isn't a valid spacepoint then stop here

      // Get the true particle, and find the closest MC trajectory point to spacepoint
      simb::MCParticle p = pi->TrackIdToParticle(abs(trueID));
      size_t closestPointIdx = GetClosestTrajectoryPoint(*sp, p);
      TVector3 closestPoint = p.Trajectory().Position(closestPointIdx).Vect();
      float dist = (TVector3(sp->XYZ())-closestPoint).Mag();
      if (dist < distCut) {
        ++nGood;
        ret[spIdx] = true;
      }
      else {
        ++nOutOfRange;
      }

      // If this was a true node & we're storing direction, store it!
      if (ret[spIdx] && dirTruth) {
        TVector3 dir = p.Trajectory().Momentum(closestPointIdx).Vect().Unit();
        for (size_t it = 0; it < 3; ++it) {
          dirTruth->at(spIdx).push_back(dir[it]);
        }
      }

    } // for spIdx

    std::cout << "There were " << nGood << " valid spacepoints, " << nIDEMismatched
      << " spacepoints with IDE mismatch, " << nNoHit
      << " spacepoints with no associated hits or noise hits, " << nMismatched
      << " spacepoints produced from mismatched hits and " << nOutOfRange
      << " spacepoints reconstructed too far away from the true particle trajectory."
      << std::endl;

    return ret;

  } // function GCNFeatureUtils::GetNodeGroundTruth

  std::map<unsigned int, unsigned int> GCNFeatureUtils::GetParticleFlowMap(
    const std::set<unsigned int> particles) const {

    art::ServiceHandle<cheat::ParticleInventoryService> pi;
    std::map<unsigned int, unsigned int> ret;
    for (int p : particles) {
      if (p == 0) continue; // No parent for the event primary
      ret[p] = pi->TrackIdToParticle(abs(p)).Mother();
    }
    return ret;

  } // function GCNFeatureUtils::GetParticleFlowMap

} // namespace cvn

