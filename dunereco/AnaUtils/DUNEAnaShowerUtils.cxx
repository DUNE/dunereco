/**
*
* @file dune/AnaUtils/DUNEAnaShowerUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Showers
*/

#include "dune/AnaUtils/DUNEAnaShowerUtils.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace dune_ana
{

std::vector<art::Ptr<recob::Hit>> DUNEAnaShowerUtils::GetHits(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label)
{    
    return DUNEAnaShowerUtils::GetAssocProductVector<recob::Hit>(pShower,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> DUNEAnaShowerUtils::GetSpacePoints(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label)
{
    return DUNEAnaShowerUtils::GetAssocProductVector<recob::SpacePoint>(pShower,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> DUNEAnaShowerUtils::GetPFParticle(const art::Ptr<recob::Shower> &pShower, const art::Event &evt, const std::string &label)
{
    return DUNEAnaShowerUtils::GetAssocProduct<recob::PFParticle>(pShower,evt,label,label);
}    

} // namespace dune_ana


