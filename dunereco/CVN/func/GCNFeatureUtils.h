////////////////////////////////////////////////////////////////////////
/// \file    GCNFeatureUtils.h
/// \brief   Utilities for calculating feature values for the GCN
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef GCN_FEATURE_UTILS_H
#define GCN_FEATURE_UTILS_H

#include <vector>
#include <string>
#include <map>

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/SpacePoint.h"

namespace cvn
{

  /// Class containing some utility functions for all things CVN
  class GCNFeatureUtils
  {
  public:
    GCNFeatureUtils();
    ~GCNFeatureUtils();

    /// Get the number of neighbours within rangeCut cm of this space point
    const unsigned int GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, const float rangeCut, const std::string &spLabel) const;
    /// Get a map of the number of neighbours for each space point ID. Using this function is much less wasteful
    /// than repeated calls to the above function
    const std::map<int,unsigned int> GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::string &spLabel) const;
    /// Use the association between space points and hits to return a charge
    const float GetSpacePointCharge(const recob::SpacePoint &sp, art::Event const &evt, const std::string &spLabel) const;

  private:

  };

}

#endif  // GCN_FEATURE_UTILS_H
