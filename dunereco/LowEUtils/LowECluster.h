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

    LowECluster(
      const std::vector<float>& position,
      float totalCharge,
      float avePeakTime,
      float purity,
      float completeness,
      const std::vector<recob::Cluster>& clusterVec);

    LowECluster(const LowECluster&);
    LowECluster& operator=(LowECluster const&);

    void initialize(
      const std::vector<float>& position,
      float totalCharge,
      float avePeakTime,
      float purity,
      float completeness,
      const std::vector<recob::Cluster>& clusterVec);

    const std::vector<float> getPosition() const { return fPosition; }
    float getX() const { return fPosition[0]; }
    float getY() const { return fPosition[1]; }
    float getZ() const { return fPosition[2]; }
    float getTotalCharge() const { return fTotalCharge; }
    float getAvePeakTime() const { return fAvePeakTime; }
    float getPurity() const { return fPurity; }
    float getCompleteness() const { return fCompleteness; }
    
    const std::vector<recob::Cluster>& getClsuters() const { return fClusterVector; }
    std::vector<recob::Cluster>& getClusters() { return fClusterVector; }

    void setPosition(const std::vector<float>& pos) const { fPosition = pos; }

    bool operator<(const solar::LowECluster& other) const
    {
      if (fPosition[2] != other.fPosition[2])
        return fPosition[2] < other.fPosition[2];
      else
        return fPosition[0] < other.fPosition[0];
    }

    friend std::ostream& operator<<(std::ostream& o, const LowECluster& c);
    //friend bool          operator <  (const LowECluster & a, const LowECluster & b);

  private:
    mutable std::vector<float> fPosition;       ///< position of this hit combination in world coordinates
    float fTotalCharge;                         ///< Sum of charges of all associated recob::Hits
    float fAvePeakTime;                         ///< Average peak time of all associated recob::Hits
    float fPurity;                              ///< Purity of the cluster wrt input signal label
    float fCompleteness;                        ///< Completeness of the cluster wrt main signal particle
    std::vector<recob::Cluster> fClusterVector; ///< Clusters comprising this 3D cluster
  };
} // namespace

#endif //LOWECLUSTER_H
