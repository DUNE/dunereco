#ifndef CVN_PROTODUNE_UTILS_H
#define CVN_PROTODUNE_UTILS_H

///////////////////////////////////////////////////////////////
// CVNProtoDUNEUtils
//  - Class to help producing images and graphs for specific
//    particles and slices for ProtoDUNE events
//  - Had to copy some of these from protoduneana to avoid
//    having a circular dependency 
//
// Leigh Whitehead - leigh.howard.whitehead@cern.ch
///////////////////////////////////////////////////////////////

#include <map>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "art/Framework/Principal/Event.h"

namespace cvn {

  class CVNProtoDUNEUtils {

  public:

    CVNProtoDUNEUtils();
    ~CVNProtoDUNEUtils();

    /*** ---------- Slice functions ---------- ***/

    // Get the hits from the slice (using the slice object)
    const std::vector<const recob::Hit*> GetRecoSliceHits(const recob::Slice &slice, art::Event const &evt, const std::string sliceModule) const;
    const std::vector<const recob::Hit*> GetRecoSliceHits(unsigned int sliceID, art::Event const &evt, const std::string sliceModule) const;

    // A map of all hits in each slice
    const std::map<unsigned int, std::vector<const recob::Hit*>> GetRecoSliceHitMap(art::Event const &evt, const std::string sliceModule) const;
    /*** ------------------------------------- ***/
  
    /*** -------- Particle functions --------- ***/

    /// Try to get the slice tagged as beam. Returns 9999 if no beam slice was found
    unsigned short GetBeamSlice(art::Event const &evt, const std::string particleLabel) const; 

    /// Get a map of slice index to the primary PFParticles within it
    const std::map<unsigned int,std::vector<const recob::PFParticle*>> GetPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const;

    /// Get a map of slice index to all of the PFParticles within it
    const std::map<unsigned int,std::vector<const recob::PFParticle*>> GetAllPFParticleSliceMap(art::Event const &evt, const std::string particleLabel) const;

    /// Get the reconstructed slice number associated with a particle
    unsigned short GetPFParticleSliceIndex(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the reconstructed slice associated with a particle
    const recob::Slice* GetPFParticleSlice(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Get the SpacePoints associated to the PFParticle
    const std::vector<const recob::SpacePoint*> GetPFParticleSpacePoints(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /// Use the pandora metadata to tell us if this is a beam particle or not
    bool IsBeamParticle(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;

    /*** ------------------------------------- ***/
  private:

    /// Helper to get the slice map and avoid code repetition
    const std::map<unsigned int,std::vector<const recob::PFParticle*>> SliceMapHelper(art::Event const &evt, const std::string particleLabel, bool primaryOnly) const;

    /// Get the metadata associated to a PFParticle from pandora
    const std::map<std::string,float> GetPFParticleMetaData(const recob::PFParticle &particle, art::Event const &evt, const std::string particleLabel) const;
  };

}

#endif
