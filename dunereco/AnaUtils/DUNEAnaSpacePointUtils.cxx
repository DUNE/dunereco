/**
*
* @file dune/AnaUtils/DUNEAnaSpacePointUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about SpacePoints
*/

#include "dune/AnaUtils/DUNEAnaSpacePointUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace dune_ana
{

std::vector<art::Ptr<recob::Hit>> DUNEAnaSpacePointUtils::GetHits(const art::Ptr<recob::SpacePoint> &pSpacepoint, const art::Event &evt, const std::string &label)
{    
    return DUNEAnaSpacePointUtils::GetAssocProductVector<recob::Hit>(pSpacepoint,evt,label,label);
}

} // namespace dune_ana


