/**
*
* @file dune/AnaUtils/DUNEAnaClusterUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Clusters
*/

#include "dune/AnaUtils/DUNEAnaClusterUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

namespace dune_ana
{

std::vector<art::Ptr<recob::Hit>> DUNEAnaClusterUtils::GetHits(const art::Ptr<recob::Cluster> &pCluster, const art::Event &evt, const std::string &label)
{    
    return DUNEAnaClusterUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
}


} // namespace dune_ana


