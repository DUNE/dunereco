#include <vector>
#include <iostream>

#include "dune/CVN/func/GCNFeatureUtils.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

cvn::GCNFeatureUtils::GCNFeatureUtils(){

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

}

// Use the association between space points and hits to return a charge
const float cvn::GCNFeatureUtils::GetSpacePointCharge(const recob::SpacePoint &sp, art::Event const &evt, const std::string &spLabel) const{

  float charge = 0.0;

  // Get the hits associated to the space points
  auto allSP = evt.getValidHandle<std::vector<recob::SpacePoint>>(spLabel);
  const art::FindManyP<recob::Hit> findHits(allSP,evt,spLabel);
  const std::vector<art::Ptr<recob::Hit>> spHits = findHits.at(sp.ID());

  for(auto const hit : spHits){
    charge += hit->Integral();
  }

  return charge;

}



