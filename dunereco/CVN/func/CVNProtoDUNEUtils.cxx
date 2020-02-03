#include "dune/CVN/func/CVNProtoDUNEUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

cvn::CVNProtoDUNEUtils::CVNProtoDUNEUtils(){

}

cvn::CVNProtoDUNEUtils::~CVNProtoDUNEUtils(){

}

// Get the hits from a given reco slice
const std::vector<const recob::Hit*> cvn::CVNProtoDUNEUtils::GetRecoSliceHits(const recob::Slice &slice, art::Event const &evt, const std::string sliceModule) const{

  return GetRecoSliceHits(slice.ID(),evt,sliceModule);  

}

// Get the reco hits but using the slice id instead
const std::vector<const recob::Hit*> cvn::CVNProtoDUNEUtils::GetRecoSliceHits(const unsigned int sliceID, art::Event const &evt, const std::string sliceModule) const{

  auto recoSlices = evt.getValidHandle<std::vector<recob::Slice> >(sliceModule);
  art::FindManyP<recob::Hit> findHits(recoSlices,evt,sliceModule);
  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(sliceID);

  std::vector<const recob::Hit*> sliceHits;

  for(const art::Ptr<recob::Hit> hit : inputHits){

    sliceHits.push_back(hit.get());

  }

  return sliceHits;

}

// Get a map of a slice number and all hits in the slice
const std::map<unsigned int, std::vector<const recob::Hit*>> cvn::CVNProtoDUNEUtils::GetRecoSliceHitMap(art::Event const &evt, const std::string sliceModule) const{

  auto recoSlices = evt.getValidHandle<std::vector<recob::Slice> >(sliceModule);
  std::map<unsigned int, std::vector<const recob::Hit*>> hitMap;

  for(auto const slice : *recoSlices){

    const std::vector<const recob::Hit*> constvec = GetRecoSliceHits(slice.ID(),evt,sliceModule);
    for(auto const h : constvec){
      hitMap[slice.ID()].push_back(h);
    }

  }

  return hitMap;

}

// Try to get the slice tagged as beam
unsigned short cvn::CVNProtoDUNEUtils::GetBeamSlice(art::Event const &evt, const std::string particleLabel) const{

  const std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap = GetPFParticleSliceMap(evt,particleLabel);

  for(auto slice : sliceMap){
    for(auto particle : slice.second){
      if(IsBeamParticle(*particle,evt,particleLabel)){
        return slice.first;
      }
    }
  }

  return 9999;

}

// Return a map of primary particles grouped by their reconstructed slice. Useful for finding slices with multiple particles
const std::map<unsigned int,std::vector<const recob::PFParticle*>> cvn::CVNProtoDUNEUtils::GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const{

  return SliceMapHelper(evt,particleLabel,true);

}

// Return a map of all particles grouped by their reconstructed slice. 
const std::map<unsigned int,std::vector<const recob::PFParticle*>> cvn::CVNProtoDUNEUtils::GetAllPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const{

  return SliceMapHelper(evt,particleLabel,false);

}

// Helper to get slice maps and avoid duplicate code
const std::map<unsigned int,std::vector<const recob::PFParticle*>> cvn::CVNProtoDUNEUtils::SliceMapHelper(art::Event const &evt, const std::string particleLabel, bool primaryOnly) const{

  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);

  std::map<unsigned int, std::vector<const recob::PFParticle*>> sliceMap;

  for(unsigned int p = 0; p < pfParticles->size(); ++p){
    const recob::PFParticle* particle = &(pfParticles->at(p));

    //  Only the primary particles have the slice association
    if(primaryOnly && !particle->IsPrimary()) continue;

    unsigned int thisSlice = GetPFParticleSliceIndex(*particle,evt,particleLabel);

    if(thisSlice != 9999){
      sliceMap[thisSlice].push_back(particle);
    }
  }

  return sliceMap;

}

// Get the space points associated to the PFParticle
const std::vector<const recob::SpacePoint*> cvn::CVNProtoDUNEUtils::GetPFParticleSpacePoints(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Get the particles and their associations
  auto particles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindManyP<recob::SpacePoint> findSpacePoints(particles,evt,particleLabel);
  const std::vector<art::Ptr<recob::SpacePoint>> pfpSpacePoints = findSpacePoints.at(particle.Self());

  // We don't want the art::Ptr so we need to get rid of it
  std::vector<const recob::SpacePoint*> sp;
  for(auto pointer : pfpSpacePoints){
    sp.push_back(pointer.get());
  }

  return sp;
}

// Get the reconstructed slice associated with a particle
const recob::Slice* cvn::CVNProtoDUNEUtils::GetPFParticleSlice(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Perhaps we should use the associations to do this? 
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  const art::FindOneP<recob::Slice> findSlice(pfParticles,evt,particleLabel);

  const recob::Slice* slice = findSlice.at(particle.Self()).get();

  return slice;
}

// Get the reconstructed slice associated with a particle
unsigned short cvn::CVNProtoDUNEUtils::GetPFParticleSliceIndex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{

  // Try to use slices if we can
  try{
    const recob::Slice* slice = GetPFParticleSlice(particle,evt,particleLabel);
    return slice->ID();
  }
  // Otherwise fall back on metadata
  catch(...){
    std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);
    std::string search = "SliceIndex";
    if(mdMap.find(search) != mdMap.end()){
      return static_cast<unsigned short>(mdMap.at(search));
    }
    else{
//    std::cerr << "Object has no slice index... returning 9999" << std::endl;
      return 9999;
    }
  }

}

const std::map<std::string,float> cvn::CVNProtoDUNEUtils::GetPFParticleMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const {
  // Get the particles
  auto pfParticles = evt.getValidHandle<std::vector<recob::PFParticle>>(particleLabel);
  // And their meta data
  const art::FindManyP<larpandoraobj::PFParticleMetadata> findMetaData(pfParticles,evt,particleLabel);

  const larpandoraobj::PFParticleMetadata metaData = *((findMetaData.at(particle.Self())).at(0));

  return metaData.GetPropertiesMap();
}

// Use the pandora metadata to tell us if this is a beam particle or not
bool cvn::CVNProtoDUNEUtils::IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const{
  std::map<std::string,float> mdMap = GetPFParticleMetaData(particle,evt,particleLabel);
  if(mdMap.find("IsTestBeam") != mdMap.end()){
    return true;
  }
  else{
    return false;
  }
}
