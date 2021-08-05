/**
*
* @file dune/AnaUtils/DUNEAnaTrackUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Tracks
*/

#include "dune/AnaUtils/DUNEAnaTrackUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "dune/TrackPID/products/CTPResult.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace dune_ana
{

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> DUNEAnaTrackUtils::GetHits(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &label)
{    
    return DUNEAnaTrackUtils::GetAssocProductVector<recob::Hit>(pTrack,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> DUNEAnaTrackUtils::GetSpacePoints(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &label)
{
    return DUNEAnaTrackUtils::GetAssocProductVector<recob::SpacePoint>(pTrack,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> DUNEAnaTrackUtils::GetPFParticle(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &label)
{
    return DUNEAnaTrackUtils::GetAssocProduct<recob::PFParticle>(pTrack,evt,label,label);
}    

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<anab::Calorimetry> DUNEAnaTrackUtils::GetCalorimetry(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &trackLabel, const std::string &caloLabel)
{
    std::vector<art::Ptr<anab::Calorimetry>> calos = DUNEAnaTrackUtils::GetAssocProductVector<anab::Calorimetry>(pTrack,evt,trackLabel,caloLabel);

    for (const art::Ptr<anab::Calorimetry> &caloObject : calos)
    {
        if (caloObject->PlaneID().Plane == geo::kW) return caloObject;
    }
    
    // If we don't have the collection plane then give up
    mf::LogError("DUNEAna") << "DUNEAnaTrackUtils::GetCalorimetry - no collection view calorimetry found. Returning first alternate view." << std::endl;

    return calos.at(0);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<ctp::CTPResult> DUNEAnaTrackUtils::GetSquIDResult(const art::Ptr<recob::Track> &pTrack, const art::Event &evt, const std::string &trackLabel, const std::string &squidLabel)
{
    return DUNEAnaTrackUtils::GetAssocProduct<ctp::CTPResult>(pTrack,evt,trackLabel,squidLabel);
}

} // namespace dune_ana

