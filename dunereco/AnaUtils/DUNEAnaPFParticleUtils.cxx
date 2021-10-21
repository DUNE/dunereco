/**
*
* @file dune/AnaUtils/DUNEAnaPFParticleUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"

namespace dune_ana
{

std::vector<art::Ptr<anab::T0>> DUNEAnaPFParticleUtils::GetT0(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProductVector<anab::T0>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<anab::CosmicTag>> DUNEAnaPFParticleUtils::GetCosmicTag(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProductVector<anab::CosmicTag>(pParticle,evt,label,label); 
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::PFParticle>> DUNEAnaPFParticleUtils::GetChildParticles(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    auto theseParticles = evt.getHandle<std::vector<recob::PFParticle>>(label);
    bool success = theseParticles.isValid();
    
    if (!success)
    {   
        mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<recob::PFParticle>>();
    }
    
    std::vector<art::Ptr<recob::PFParticle>> children;

    for (unsigned int iPart = 0; iPart < theseParticles->size(); ++iPart)
    {     
        if (theseParticles->at(iPart).Parent() == pParticle.key())
        {
            art::Ptr<recob::PFParticle> pChild(theseParticles,iPart);
            children.push_back(pChild);
        }
    }

    return children;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> DUNEAnaPFParticleUtils::GetHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{    
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<art::Ptr<recob::Cluster>> theseClusters = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<art::Ptr<recob::Hit>> theseHits;
    for (const art::Ptr<recob::Cluster> pCluster : theseClusters)
    {
      const std::vector<art::Ptr<recob::Hit>> tempHits = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
      theseHits.insert(theseHits.end(),tempHits.begin(),tempHits.end());
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> DUNEAnaPFParticleUtils::GetViewHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label, const unsigned short &view)
{
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<art::Ptr<recob::Cluster>> theseClusters = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<art::Ptr<recob::Hit>> theseHits;
    for (const art::Ptr<recob::Cluster> pCluster : theseClusters)
    { 
        const std::vector<art::Ptr<recob::Hit>> tempHits = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
        for(const art::Ptr<recob::Hit> pHit : tempHits)
        {
            if (pHit->View() == view)
            {
                theseHits.push_back(pHit);
            }
        }
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> DUNEAnaPFParticleUtils::GetSpacePoints(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProductVector<recob::SpacePoint>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Track> DUNEAnaPFParticleUtils::GetTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<recob::Track>(pParticle,evt,particleLabel,trackLabel);
}    

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Shower> DUNEAnaPFParticleUtils::GetShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Vertex> DUNEAnaPFParticleUtils::GetVertex(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<recob::Vertex>(pParticle,evt,particleLabel,particleLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Slice> DUNEAnaPFParticleUtils::GetSlice(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> pMetadata = DUNEAnaPFParticleUtils::GetMetadata(pParticle,evt,label);
    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    unsigned int sliceIndex;
    std::map<std::string,float>::iterator mapItr = metaMap.find("SliceIndex");

    if (mapItr == metaMap.end())   
    {
        throw cet::exception("DUNEAna") << " DUNEAnaPFParticleUtils::GetSlice --- No associated slice found";
    }
    else
    {
        sliceIndex = mapItr->second;
    }

    return DUNEAnaPFParticleUtils::GetProductVector<recob::Slice>(evt,label).at(sliceIndex);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<larpandoraobj::PFParticleMetadata> DUNEAnaPFParticleUtils::GetMetadata(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<larpandoraobj::PFParticleMetadata>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    // This function needs to fail if GetTrack would fail
    const std::vector<art::Ptr<recob::Track>> theseTracks = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Track>(pParticle,evt,particleLabel,trackLabel);

    return !theseTracks.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    const std::vector<art::Ptr<recob::Shower>> theseShowers = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
   
    return !theseShowers.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsClearCosmic(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> pMetadata = DUNEAnaPFParticleUtils::GetMetadata(pParticle,evt,particleLabel);

    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    return metaMap.find("IsClearCosmic") != metaMap.end();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsNeutrino(const art::Ptr<recob::PFParticle> &pParticle)
{
    const int pdg = pParticle->PdgCode();
    return ((std::abs(pdg) == 12) || (std::abs(pdg) == 14) || (std::abs(pdg) ==16));
}

} // namespace dune_ana


