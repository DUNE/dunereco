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
    return DUNEAnaTrackUtils::GetAssocProduct<anab::Calorimetry>(pTrack,evt,trackLabel,caloLabel);
}

} // namespace dune_ana

