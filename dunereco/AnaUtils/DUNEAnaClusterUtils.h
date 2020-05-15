/**
 *
 * @file dune/AnaUtils/DUNEAnaClusterUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Clusters
*/

#ifndef DUNE_ANA_CLUSTER_UTILS_H
#define DUNE_ANA_CLUSTER_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "dune/AnaUtils/DUNEAnaUtilsBase.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"

#include <string>
#include <vector>

namespace dune_ana
{
/**
 *
 * @brief DUNEAnaClusterUtils class
 *
*/
class DUNEAnaClusterUtils:DUNEAnaUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the cluster.
    *
    * @param cluster is the cluster for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the cluster producer
    * 
    * @return vector of art::Ptrs to the hits 
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Cluster> &pCluster, const art::Event &evt, const std::string &label);
};

} // namespace dune_ana

#endif // DUNE_ANA_CLUSTER_UTILS_H

