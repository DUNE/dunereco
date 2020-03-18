/**
 *
 * @file dune/AnaUtils/DUNEAnaSpacePointUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about SpacePoints
*/

#ifndef DUNE_ANA_SPACEPOINT_UTILS_H
#define DUNE_ANA_SPACEPOINT_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "dune/AnaUtils/DUNEAnaUtilsBase.h"

#include <string>
#include <vector>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"


namespace dune_ana
{
/**
 *
 * @brief DUNEAnaSpacePointUtils class
 *
*/
class DUNEAnaSpacePointUtils:DUNEAnaUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the spacepoint.
    *
    * @param spacepoint is the spacepoint for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the spacepoint producer
    * 
    * @return vector of art::Ptrs to the hits 
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::SpacePoint> &pSpacepoint, const art::Event &evt, const std::string &label);
};

} // namespace dune_ana


#endif // DUNE_ANA_SPACEPOINT_UTILS_H

