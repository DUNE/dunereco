/**
*
* @file dune/AnaUtils/DUNEAnaEventUtils.cxx
*
* @brief Utility containing helpful functions for end users to access products from events
*/

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace dune_ana
{

std::vector<art::Ptr<recob::PFParticle>> DUNEAnaEventUtils::GetPFParticles(const art::Event &evt, const std::string &label)
{
    return DUNEAnaEventUtils::GetProductVector<recob::PFParticle>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Track>> DUNEAnaEventUtils::GetTracks(const art::Event &evt, const std::string &label)
{
    mf::LogWarning("DUNEAna") << " Please note: accessing PFParticle tracks through this method is not the recommended workflow.\n"
                                 << " Please use DUNEAnaEventUtils::GetPFParticles and access the tracks with DUNEAnaPFParticleUtils::GetTrack."
                                 << std::endl;

    return DUNEAnaEventUtils::GetProductVector<recob::Track>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Shower>> DUNEAnaEventUtils::GetShowers(const art::Event &evt, const std::string &label)
{
    mf::LogWarning("DUNEAna") << " Please note: accessing PFParticle showers through this method is not the recommended workflow.\n"
                                 << " Please use DUNEAnaEventUtils::GetPFParticles and access the tracks with DUNEAnaPFParticleUtils::GetShower."
                                 << std::endl;

    return DUNEAnaEventUtils::GetProductVector<recob::Shower>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Vertex>> DUNEAnaEventUtils::GetVertices(const art::Event &evt, const std::string &label)
{
    return DUNEAnaEventUtils::GetProductVector<recob::Vertex>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> DUNEAnaEventUtils::GetSpacePoints(const art::Event &evt, const std::string &label)
{
    return DUNEAnaEventUtils::GetProductVector<recob::SpacePoint>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Slice>> DUNEAnaEventUtils::GetSlices(const art::Event &evt, const std::string &label)
{
    return DUNEAnaEventUtils::GetProductVector<recob::Slice>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<simb::MCParticle>> DUNEAnaEventUtils::GetMCParticles(const art::Event &evt, const std::string &label)
{
    return DUNEAnaEventUtils::GetProductVector<simb::MCParticle>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::PFParticle>> DUNEAnaEventUtils::GetClearCosmics(const art::Event &evt, const std::string &label)
{
    std::vector<art::Ptr<recob::PFParticle>> theseParticles = DUNEAnaEventUtils::GetProductVector<recob::PFParticle>(evt,label);

    std::vector<art::Ptr<recob::PFParticle>> theseCosmics;

    // We only want primary cosmic rays
    for (art::Ptr<recob::PFParticle> pParticle : theseParticles)
    {
        if (!pParticle->IsPrimary())
        {
            continue;
        } 

        if (DUNEAnaPFParticleUtils::IsClearCosmic(pParticle, evt, label))
        {
            theseCosmics.push_back(pParticle);
        }
    }

    return theseCosmics;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> DUNEAnaEventUtils::GetNeutrino(const art::Event &evt, const std::string &label)
{
    if (!HasNeutrino(evt,label))
    {
        throw cet::exception("DUNEAna") << "DUNEAnaEventUtils::GetNeutrino --- No neutrino found";
    }
    
    art::Ptr<recob::PFParticle> pNeutrino;
    std::vector<art::Ptr<recob::PFParticle>> particles = DUNEAnaEventUtils::GetPFParticles(evt,label);
    for (art::Ptr<recob::PFParticle> pParticle : particles)
    {
        if (DUNEAnaPFParticleUtils::IsNeutrino(pParticle))
        {
            pNeutrino = pParticle;
            break; 
        }
    }
    return pNeutrino;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaEventUtils::HasNeutrino(const art::Event &evt, const std::string &label)
{
    bool hasNeutrino = false;
    std::vector<art::Ptr<recob::PFParticle>> particles = DUNEAnaEventUtils::GetPFParticles(evt,label);
    for (art::Ptr<recob::PFParticle> pParticle : particles)
    {
        if(!pParticle->IsPrimary())
        {
            continue;
        }
        if (DUNEAnaPFParticleUtils::IsNeutrino(pParticle))
        {
            hasNeutrino = true;
            break;
        }
    }
    return hasNeutrino;
}

} // namespace dune_ana

