////////////////////////////////////////////////////////////////////////////
//
// \brief Definition of cluster object for Solar neutrino studies in DUNE
//
// \author sergio.manthey@ciemat.es
//
////////////////////////////////////////////////////////////////////////////

#include "LowECluster.h"

namespace solar {

  LowECluster::LowECluster()
    : fPosition(std::vector<float>(3, -1e6))
    , fTotalCharge(-1)
    , fAvePeakTime(-1e6)
    , fPurity(-1)
    , fCompleteness(-1)
  {
    fClusterVector.clear();
    fClusterVector.resize(3, recob::Cluster());
  }

  LowECluster::LowECluster(
    const std::vector<float>& position,
    float totalCharge,
    float avePeakTime,
    float purity,
    float completeness,
    const std::vector<recob::Cluster>& clusterVec)
    : fPosition(position)
    , fTotalCharge(totalCharge)
    , fAvePeakTime(avePeakTime)
    , fPurity(purity)
    , fCompleteness(completeness)
  {
    fClusterVector.resize(3, recob::Cluster());
    std::copy(clusterVec.begin(), clusterVec.end(), fClusterVector.begin());
  }

  LowECluster::LowECluster(const LowECluster& toCopy)
  {
    fPosition = toCopy.fPosition;
    fTotalCharge = toCopy.fTotalCharge;
    fAvePeakTime = toCopy.fAvePeakTime;
    fPurity = toCopy.fPurity;
    fCompleteness = toCopy.fCompleteness;
    fClusterVector = toCopy.fClusterVector;
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
    float totalCharge,
    float avePeakTime,
    float purity,
    float completeness,
    const std::vector<recob::Cluster>& clusterVec)
  {
    fPosition = position;
    fTotalCharge = totalCharge;
    fAvePeakTime = avePeakTime;
    fPurity = purity;
    fCompleteness = completeness;
    fClusterVector = clusterVec;

    return;
  }
} // namespace
