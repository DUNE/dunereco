////////////////////////////////////////////////////////////////////////////
//
// Definition of cluster object for low energy neutrino studies in DUNE
//  
//
// \author sergio.manthey@ciemat.es
//
////////////////////////////////////////////////////////////////////////////

#include "LowECluster.h"

namespace solar {

  LowECluster::LowECluster()
    : averagePosition(std::vector<float>(3, -1e6))
    , nHits(0)
    , mainChannel(-1)
    , totalCharge(-1)
    , averageTime(-1e6)
    , averagePurity(-1)
    , totalCompleteness(-1)
  {
    clusterVector.clear();
    clusterVector.resize(3, recob::Cluster());
  }

  LowECluster::LowECluster(
    const std::vector<float>& position,
    int   nhits,
    int   channel,
    float charge,
    float time,
    float purity,
    float completeness,
    const std::vector<recob::Cluster>& clusters)
    : averagePosition(position)
    , nHits(nhits)
    , mainChannel(channel)
    , totalCharge(charge)
    , averageTime(time)
    , averagePurity(purity)
    , totalCompleteness(completeness)
  {
    clusterVector.resize(3, recob::Cluster());
    std::copy(clusters.begin(), clusters.end(), clusterVector.begin());
  }

  LowECluster::LowECluster(const LowECluster& toCopy)
  {
    averagePosition = toCopy.averagePosition;
    nHits = toCopy.nHits;
    mainChannel = toCopy.mainChannel;
    totalCharge = toCopy.totalCharge;
    averageTime = toCopy.averageTime;
    averagePurity = toCopy.averagePurity;
    totalCompleteness = toCopy.totalCompleteness;
    clusterVector = toCopy.clusterVector;
  }

  LowECluster& LowECluster::operator=(LowECluster const& toCopy)
  {
    using std::swap;
    auto tmp = toCopy;
    swap(tmp, *this);
    return *this;
  }

  void LowECluster::initialize(
    const std::vector<float>& position,
    int   nhits,
    int   channel,
    float charge,
    float time,
    float purity,
    float completeness,
    const std::vector<recob::Cluster>& clusters)
  {
    averagePosition = position;
    nHits = nhits;
    mainChannel = channel;
    totalCharge = charge;
    averageTime = time;
    averagePurity = purity;
    totalCompleteness = completeness;
    clusterVector = clusters;

    return;
  }
} // namespace
