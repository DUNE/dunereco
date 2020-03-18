/**
 *
 * @file dune/AnaUtils/DUNEAnaPFParticleUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#ifndef DUNE_ANA_PFPARTICLE_UTILS_H
#define DUNE_ANA_PFPARTICLE_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "dune/AnaUtils/DUNEAnaUtilsBase.h"

#include <string>
#include <vector>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"

namespace dune_ana
{
/**
 *
 * @brief DUNEAnaPFParticleUtils class
 *
*/
class DUNEAnaPFParticleUtils:DUNEAnaUtilsBase
{
public:
    /**
    * @brief Get the T0(s) associated with the particle.
    *
    * @param pParticle particle for which we want the T0
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return vector of art::Ptrs to the T0(s)
    */ 
    static std::vector<art::Ptr<anab::T0>> GetT0(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the Cosmic Tag(s) associated with the particle.
    *
    * @param pParticle particle for which we want the cosmic tag
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return vector of art::Ptrs to the cosmic tag(s) 
    */
    static std::vector<art::Ptr<anab::CosmicTag>> GetCosmicTag(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the child particles (one step down in the hierarchy) of this particle.
    *
    * @param pParticle particle for which we want the child particles
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return vector of art::Ptrs to the child particles 
    */
    static std::vector<art::Ptr<recob::PFParticle>> GetChildParticles(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the hits associated to this particle.
    *
    * @param pParticle particle for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return vector of art::Ptrs to the hits 
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the hits associated to this particle in a given view.
    *
    * @param pParticle particle for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * @param view is the view for which we want the hits 
    *
    * @return vector of art::Ptrs to the hits 
    */
    static std::vector<art::Ptr<recob::Hit>> GetViewHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label, const unsigned short &view);

    /**
    * @brief Get the spacepoints associated to this particle.
    *
    * @param pParticle particle for which we want the spacepoints
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return vector of art::Ptrs to the spacepoints 
    */
    static std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label);

    /**
    * @brief Get the track associated to this particle. Should only be called if IsTrack method succeeds
    *
    * @param pParticle particle for which we want the track
    * @param evt is the underlying art event
    * @param particleLabel is the label for the PFParticle producer
    * @param trackLabel is the label for the Track producer
    *
    * @return art::Ptr to the track 
    */
    static art::Ptr<recob::Track> GetTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel, const std::string &trackLabel);

    /**
    * @brief Get the shower associated to this particle. Should only be called if IsShower method succeeds
    *
    * @param pParticle particle for which we want the shower
    * @param evt is the underlying art event
    * @param particleLabel is the label for the PFParticle producer
    * @param trackLabel is the label for the Shower producer
    *
    * @return art::Ptr to the shower
    */
    static art::Ptr<recob::Shower> GetShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel, const std::string &showerLabel);
    
    /**
    * @brief Get the vertex associated to this particle.
    *
    * @param pParticle particle for which we want the vertex
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return art::Ptr to the vertex
    */
    static art::Ptr<recob::Vertex> GetVertex(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel);

    /**
    * @brief Get the slice associated to this particle.
    *
    * @param pParticle particle for which we want the slice
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return art::Ptr to the slice
    */
    static art::Ptr<recob::Slice> GetSlice(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel);

    /**
    * @brief Get the metadata associated to this particle.
    *
    * @param pParticle particle for which we want the slice
    * @param evt is the underlying art event
    * @param label is the label for the PFParticle producer
    * 
    * @return art::Ptr to the metadata
    */
    static art::Ptr<larpandoraobj::PFParticleMetadata> GetMetadata(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label);

    /**
    * @brief Check if this particle has an associated track
    *
    * @param pParticle particle for which we want the track-like confirmation
    * @param evt is the underlying art event
    * @param particleLabel is the label for the PFParticle producer
    * @param trackLabel is the label for the track producer
    *
    * @return bool stating if the object is track-like
    */
    static bool IsTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel, const std::string &trackLabel);

    /**
    * @brief Check if this particle has an associated shower
    *
    * @param pParticle particle for which we want the shower-like confirmation
    * @param evt is the underlying art event
    * @param particleLabel is the label for the PFParticle producer
    * @param showerLabel is the label for the shower producer
    *
    * @return bool stating if the object is shower-like
    */
    static bool IsShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel, const std::string &showerLabel);

    /**
    * @brief Check if this particle is a clear cosmic ray
    *
    * @param pParticle particle for which we want the clear cosmic confirmation
    * @param evt is the underlying art event
    * @param particleLabel is the label for the PFParticle producer
    *
    * @return bool stating if the object is a clear cosmic ray
    */
    static bool IsClearCosmic(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &pParticleLabel);

    /**
    * @brief Check if this particle is a neutrino
    *
    * @param pParticle particle for which we want the neutrino confirmation
    *
    * @return bool stating if the object is a neutrino
    */
    static bool IsNeutrino(const art::Ptr<recob::PFParticle> &pParticle);
};

} // namespace dune_ana

#endif // DUNE_ANA_PFPARTICLE_UTILS_H

