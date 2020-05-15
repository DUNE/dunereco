/**
*
* @file dune/AnaUtils/DUNEAnaSliceUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about Slices
*/

#include "dune/AnaUtils/DUNEAnaSliceUtils.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

namespace dune_ana
{

std::vector<art::Ptr<recob::Hit>> DUNEAnaSliceUtils::GetHits(const art::Ptr<recob::Slice> &pSlice, const art::Event &evt, const std::string &label)
{    
    return DUNEAnaSliceUtils::GetAssocProductVector<recob::Hit>(pSlice,evt,label,label);
}

} // namespace dune_ana


