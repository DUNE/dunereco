/**
 *
 * @file dunereco/AnaUtils/DUNEAnaEventUtils.h
 *
 * @brief Utility containing helpful functions for end users to access products from events
*/

#ifndef DUNE_ANA_EVENT_UTILS_H
#define DUNE_ANA_EVENT_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "dunereco/AnaUtils/DUNEAnaUtilsBase.h"
#include "dunereco/CVN/func/Result.h"
#include "dunereco/TrackPID/products/CTPResult.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <string>
#include <vector>

namespace dune_ana
{
/**
 *
 * @brief DUNEAnaEventUtils class
 *
*/
class DUNEAnaEventUtils:DUNEAnaUtilsBase
{
public:
    /**
    * @brief Get the particles from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return vector of art::Ptrs to particles
    */
    static std::vector<art::Ptr<recob::PFParticle>> GetPFParticles(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the tracks from the event. This function shouldn't be used as the basis of an analysis.
    *
    * @param evt is the underlying art event
    * @param label is the label for the track producer
    *
    * @return vector of art::Ptrs to tracks
    */
    static std::vector<art::Ptr<recob::Track>> GetTracks(const art::Event &evt, const std::string &label);
    
    /**
    * @brief Get the showers from the event. This function shouldn't be used as the basis of an analysis.
    *
    * @param evt is the underlying art event
    * @param label is the label for the shower producer
    *
    * @return vector of art::Ptrs to showers
    */
    static std::vector<art::Ptr<recob::Shower>> GetShowers(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the vertices from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the vertex producer
    *
    * @return vector of art::Ptrs to vertices
    */
    static std::vector<art::Ptr<recob::Vertex>> GetVertices(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the hits from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the hit producer
    *
    * @return vector of art::Ptrs to spacepoints
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the wires from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the wire producer
    *
    * @return vector of art::Ptrs to wires
    */
    static std::vector<art::Ptr<recob::Wire>> GetWires(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the spacepoints from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the spacepoint producer
    *
    * @return vector of art::Ptrs to spacepoints
    */
    static std::vector<art::Ptr<recob::SpacePoint>> GetSpacePoints(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the slices from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the slice producer
    *
    * @return vector of art::Ptrs to slices
    */
    static std::vector<art::Ptr<recob::Slice>> GetSlices(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the MC truths from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the MC truth producer
    *
    * @return vector of art::Ptrs to MC truths
    */
    static std::vector<art::Ptr<simb::MCTruth>> GetMCTruths(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the MC particles from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the MC particle producer
    *
    * @return vector of art::Ptrs to MC particles
    */
    static std::vector<art::Ptr<simb::MCParticle>> GetMCParticles(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the clear cosmic ray primaries from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return vector of art::Ptrs to cosmic particles
    */
    static std::vector<art::Ptr<recob::PFParticle>> GetClearCosmics(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the neutrino from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return atr::Ptr to the neutrino
    */
    static art::Ptr<recob::PFParticle> GetNeutrino(const art::Event &evt, const std::string &label);

    /**
    * @brief Check to see if the event has a reconstructed neutrino
    *
    * @param evt is the underlying art event
    * @param label is the label for the particle producer
    *
    * @return true if the event has a reconstructed neutrino
    */
    static bool HasNeutrino(const art::Event &evt, const std::string &label);    

    /**
    * @brief Get the CVN result objects from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the cvn evaluator
    *
    * @return vector of art::Ptrs to CVN result objects
    */
    static std::vector<art::Ptr<cvn::Result>> GetCVNResults(const art::Event &evt, const std::string &label);

    /**
    * @brief Get the SquID track PID result objects from the event
    *
    * @param evt is the underlying art event
    * @param label is the label for the SquID evaluator
    *
    * @return vector of art::Ptrs to SquID result objects
    */
    static std::vector<art::Ptr<ctp::CTPResult>> GetSquIDResults(const art::Event &evt, const std::string &label);


};

} // namespace dune_ana

#endif // DUNE_ANA_EVENT_UTILS_H

