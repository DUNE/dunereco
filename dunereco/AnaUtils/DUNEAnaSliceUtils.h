/**
 *
 * @file dune/AnaUtils/DUNEAnaSliceUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Slices
*/

#ifndef DUNE_ANA_SLICE_UTILS_H
#define DUNE_ANA_SLICE_UTILS_H

#include "art/Framework/Principal/Event.h"

#include "dune/AnaUtils/DUNEAnaUtilsBase.h"

#include <string>
#include <vector>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

namespace dune_ana
{
/**
 *
 * @brief DUNEAnaSliceUtils class
 *
*/
class DUNEAnaSliceUtils:DUNEAnaUtilsBase
{
public:
    /**
    * @brief Get the hits associated with the slice.
    *
    * @param slice is the slice for which we want the hits
    * @param evt is the underlying art event
    * @param label is the label for the slice producer
    * 
    * @return vector of art::Ptrs to the hits
    */
    static std::vector<art::Ptr<recob::Hit>> GetHits(const art::Ptr<recob::Slice> &pSlice, const art::Event &evt, const std::string &label);
};

} // namespace dune_ana

#endif // DUNE_ANA_SLICE_UTILS_H

