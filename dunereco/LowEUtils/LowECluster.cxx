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
    : fID(std::numeric_limits<size_t>::max())
    , fPosition(std::vector<float>(3, -1e6))
    , fTotalCharge(0.)
    , fAvePeakTime(-1.)
    , fDeltaPeakTime(0.)
    , fSigmaPeakTime(0.)
    , fChargeAsymmetry(0.)
  {
    fHitDelTSigVec.clear();
    fWireIDVector.clear();
    fHitVector.clear();
    fHitDelTSigVec.resize(3, 0.);
    fWireIDVector.resize(3, geo::WireID());
    fHitVector.resize(3, std::vector<recob::Hit>());
  }

  LowECluster::LowECluster(size_t id,
                          const std::vector<float>& position,
                          float totalCharge,
                          float avePeakTime,
                          float deltaPeakTime,
                          float sigmaPeakTime,
                          float chargeAsymmetry,
                          const std::vector<std::vector<recob::Hit>>& hitVec,
                          const std::vector<float>& hitDelTSigVec,
                          const std::vector<geo::WireID>& wireIDs)
    : fID(id)
    , fPosition(position)
    , fTotalCharge(totalCharge)
    , fAvePeakTime(avePeakTime)
    , fDeltaPeakTime(deltaPeakTime)
    , fSigmaPeakTime(sigmaPeakTime)
    , fChargeAsymmetry(chargeAsymmetry)
    , fHitDelTSigVec(hitDelTSigVec)
    , fWireIDVector(wireIDs)
  {
    fHitVector.resize(3, std::vector<recob::Hit>());
    std::copy(hitVec.begin(), hitVec.end(), fHitVector.begin());
  }

  LowECluster::LowECluster(const LowECluster& toCopy)
  {
    fID = toCopy.fID;
    fPosition = toCopy.fPosition;
    fTotalCharge = toCopy.fTotalCharge;
    fAvePeakTime = toCopy.fAvePeakTime;
    fDeltaPeakTime = toCopy.fDeltaPeakTime;
    fSigmaPeakTime = toCopy.fSigmaPeakTime;
    fChargeAsymmetry = toCopy.fChargeAsymmetry;
    fHitVector = toCopy.fHitVector;
    fHitDelTSigVec = toCopy.fHitDelTSigVec;
    fWireIDVector = toCopy.fWireIDVector;
  }

  LowECluster& LowECluster::operator=(LowECluster const& toCopy)
  {
    using std::swap;
    auto tmp = toCopy;
    swap(tmp, *this);
    return *this;
  }

  void LowECluster::initialize(size_t id,
                                const std::vector<float>& position,
                                float totalCharge,
                                float avePeakTime,
                                float deltaPeakTime,
                                float sigmaPeakTime,
                                float chargeAsymmetry,
                                const std::vector<std::vector<recob::Hit>>& hitVec,
                                const std::vector<float>& hitDelTSigVec,
                                const std::vector<geo::WireID>& wireIDs)
  {
    fID = id;
    fPosition = position;
    fTotalCharge = totalCharge;
    fAvePeakTime = avePeakTime;
    fDeltaPeakTime = deltaPeakTime;
    fSigmaPeakTime = sigmaPeakTime;
    fChargeAsymmetry = chargeAsymmetry;
    fHitVector = hitVec;
    fHitDelTSigVec = hitDelTSigVec;
    fWireIDVector = wireIDs;

    return;
  }

  void LowECluster::setWireID(const geo::WireID& wid) const
  {
    fWireIDVector[wid.Plane] = wid;
  }

} // namespace
