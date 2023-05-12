/**
*
* @file dunereco/AnaUtils/DUNEAnaPFParticleUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

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

std::vector<anab::T0> DUNEAnaPFParticleUtils::GetT0(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProductVector<anab::T0>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<anab::CosmicTag> DUNEAnaPFParticleUtils::GetCosmicTag(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProductVector<anab::CosmicTag>(pParticle,evt,label,label); 
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<recob::PFParticle> DUNEAnaPFParticleUtils::GetChildParticles(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{
    auto theseParticles = evt.getHandle<std::vector<recob::PFParticle>>(label);
    bool success = theseParticles.isValid();
    
    if (!success)
    {   
        mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<recob::PFParticle>();
    }
    
    std::vector<recob::PFParticle> children;

    for (unsigned int iPart = 0; iPart < theseParticles->size(); ++iPart)
    {     
        if (theseParticles->at(iPart).Parent() == pParticle.key())
        {
            recob::PFParticle pChild(theseParticles,iPart);
            children.push_back(pChild);
        }
    }

    return children;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<recob::Hit> DUNEAnaPFParticleUtils::GetHits(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{    
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<recob::Cluster> theseClusters = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<recob::Hit> theseHits;
    for (const recob::Cluster pCluster : theseClusters)
    {
      const std::vector<recob::Hit> tempHits = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
      theseHits.insert(theseHits.end(),tempHits.begin(),tempHits.end());
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<recob::Hit> DUNEAnaPFParticleUtils::GetViewHits(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label, const unsigned short &view)
{
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<recob::Cluster> theseClusters = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<recob::Hit> theseHits;
    for (const recob::Cluster pCluster : theseClusters)
    { 
        const std::vector<recob::Hit> tempHits = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
        for(const recob::Hit pHit : tempHits)
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

std::vector<recob::SpacePoint> DUNEAnaPFParticleUtils::GetSpacePoints(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProductVector<recob::SpacePoint>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

recob::Track DUNEAnaPFParticleUtils::GetTrack(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<recob::Track>(pParticle,evt,particleLabel,trackLabel);
}    

//-----------------------------------------------------------------------------------------------------------------------------------------

recob::Shower DUNEAnaPFParticleUtils::GetShower(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

recob::Vertex DUNEAnaPFParticleUtils::GetVertex(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<recob::Vertex>(pParticle,evt,particleLabel,particleLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

recob::Slice DUNEAnaPFParticleUtils::GetSlice(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{
    const larpandoraobj::PFParticleMetadata pMetadata = DUNEAnaPFParticleUtils::GetMetadata(pParticle,evt,label);
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

larpandoraobj::PFParticleMetadata DUNEAnaPFParticleUtils::GetMetadata(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &label)
{
    return DUNEAnaPFParticleUtils::GetAssocProduct<larpandoraobj::PFParticleMetadata>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsTrack(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    // This function needs to fail if GetTrack would fail
    const std::vector<recob::Track> theseTracks = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Track>(pParticle,evt,particleLabel,trackLabel);

    return !theseTracks.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsShower(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    const std::vector<recob::Shower> theseShowers = DUNEAnaPFParticleUtils::GetAssocProductVector<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
   
    return !theseShowers.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsClearCosmic(const recob::PFParticle &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    const larpandoraobj::PFParticleMetadata pMetadata = DUNEAnaPFParticleUtils::GetMetadata(pParticle,evt,particleLabel);

    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    return metaMap.find("IsClearCosmic") != metaMap.end();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool DUNEAnaPFParticleUtils::IsNeutrino(const recob::PFParticle &pParticle)
{
    const int pdg = pParticle->PdgCode();
    return ((std::abs(pdg) == 12) || (std::abs(pdg) == 14) || (std::abs(pdg) ==16));
}

} // namespace dune_ana


