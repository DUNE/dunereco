/**
 *
 *  @brief Definition of reco objects for use in the Solar neutrino studies for DUNE
 *
 *  @author sergio.manthey@ciemat.es
 *
 */

#ifndef LOWECLUSTER_H
#define LOWECLUSTER_H

#include <vector>

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace solar {
  // Now define an object with the recob::Hit information that will comprise the 3D cluster
  class LowECluster {
  public:
    LowECluster(); // Default constructor

    LowECluster(size_t id,
                 const std::vector<float>& position,
                 float totalCharge,
                 float avePeakTime,
                 float deltaPeakTime,
                 float sigmaPeakTime,
                 float chargeAsymmetry,
                 const std::vector<std::vector<recob::Hit>>& hitVec,
                 const std::vector<float>& hitDelTSigVec,
                 const std::vector<geo::WireID>& wireIDVec);

    LowECluster(const LowECluster&);
    LowECluster& operator=(LowECluster const&);

    void initialize(size_t id,
                 const std::vector<float>& position,
                 float totalCharge,
                 float avePeakTime,
                 float deltaPeakTime,
                 float sigmaPeakTime,
                 float chargeAsymmetry,
                 const std::vector<std::vector<recob::Hit>>& hitVec,
                 const std::vector<float>& hitDelTSigVec,
                 const std::vector<geo::WireID>& wireIDVec);

    size_t getID() const { return fID; }
    const std::vector<float> getPosition() const { return fPosition; }
    float getX() const { return fPosition[0]; }
    float getY() const { return fPosition[1]; }
    float getZ() const { return fPosition[2]; }
    float getTotalCharge() const { return fTotalCharge; }
    float getAvePeakTime() const { return fAvePeakTime; }
    float getDeltaPeakTime() const { return fDeltaPeakTime; }
    float getSigmaPeakTime() const { return fSigmaPeakTime; }
    float getChargeAsymmetry() const { return fChargeAsymmetry; }
    const std::vector<std::vector<recob::Hit>>& getHits() const { return fHitVector; }
    const std::vector<float> getHitDelTSigVec() const { return fHitDelTSigVec; }
    const std::vector<geo::WireID>& getWireIDs() const { return fWireIDVector; }

    std::vector<std::vector<recob::Hit>>& getHits() { return fHitVector; }

    void setID(const size_t& id) const { fID = id; }
    void setWireID(const geo::WireID& wid) const;

    void setPosition(const std::vector<float>& pos) const { fPosition = pos; }

    bool operator<(const solar::LowECluster& other) const
    {
      if (fPosition[2] != other.fPosition[2])
        return fPosition[2] < other.fPosition[2];
      else
        return fPosition[0] < other.fPosition[0];
    }

    bool operator==(const solar::LowECluster& other) const { return fID == other.fID; }

    friend std::ostream& operator<<(std::ostream& o, const LowECluster& c);
    //friend bool          operator <  (const LowECluster & a, const LowECluster & b);

  private:
    mutable size_t fID;                              ///< "id" of this hit (useful for indexing)
    mutable std::vector<float> fPosition;            ///< position of this hit combination in world coordinates
    float fTotalCharge;                              ///< Sum of charges of all associated recob::Hits
    float fAvePeakTime;                              ///< Average peak time of all associated recob::Hits
    float fDeltaPeakTime;                            ///< Largest delta peak time of associated recob::Hits
    float fSigmaPeakTime;                            ///< Quad sum of peak time sigmas
    float fChargeAsymmetry;                          ///< Assymetry of average of two closest to third charge
    std::vector<std::vector<recob::Hit>> fHitVector; ///< Hits comprising this 3D hit
    mutable std::vector<float> fHitDelTSigVec;       ///< Delta t of hit to matching pair / sig
    mutable std::vector<geo::WireID> fWireIDVector;  ///< Wire ID's for the planes making up hit
  };
} // namespace

#endif //LOWECLUSTER_H
