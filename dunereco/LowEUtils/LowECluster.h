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
      int   id,
      int   nhits,
      int   channel,
      float charge,
      float time,
      float purity,
      float completeness,
      const std::vector<recob::Cluster>& clusters);

    LowECluster(const LowECluster&);
    LowECluster& operator=(LowECluster const&);

    // Destructor
    ~LowECluster() = default; // Default destructor

    void initialize(
      const std::vector<float>& position,
      int   id,
      int   nhits,
      int   channel,
      float charge,
      float time,
      float purity,
      float completeness,
      const std::vector<recob::Cluster>& clusters);

    // Getters
    const std::vector<float> getPosition() const { return averagePosition; }
    float getX() const { return averagePosition[0]; }
    float getY() const { return averagePosition[1]; }
    float getZ() const { return averagePosition[2]; }
    int   getMainID() const { return mainID; }
    int   getNHits() const { return nHits; }
    int   getMainChannel() const { return mainChannel; }
    float getTotalCharge() const { return totalCharge; }
    float getAverageTime() const { return averageTime; }
    float getPurity() const { return averagePurity; }
    float getCompleteness() const { return totalCompleteness; }
    
    const std::vector<recob::Cluster>& getClusters() const { return clusterVector; }
    std::vector<recob::Cluster>& getClusters() { return clusterVector; }

    // Setters
    void setPosition(const std::vector<float>& pos) const { averagePosition = pos; }
    void setMainID(int id) { mainID = id; }
    void setNHits(int nhits) { nHits = nhits; }
    void setMainChannel(int channel) { mainChannel = channel; }
    void setTotalCharge(float charge) { totalCharge = charge; }
    void setAverageTime(float time) { averageTime = time; }
    void setPurity(float purity) { averagePurity = purity; }
    void setCompleteness(float completeness) { totalCompleteness = completeness; }
    void setClusters(const std::vector<recob::Cluster>& clusters) { clusterVector = clusters; }
    void addCluster(const recob::Cluster& cluster) { clusterVector.push_back(cluster); }

    bool operator<(const solar::LowECluster& other) const
    {
      if (averagePosition[2] != other.averagePosition[2])
        return averagePosition[2] < other.averagePosition[2];
      else
        return averagePosition[0] < other.averagePosition[0];
    }

    friend std::ostream& operator<<(std::ostream& o, const LowECluster& c);
    //friend bool          operator <  (const LowECluster & a, const LowECluster & b);

  private:
    mutable std::vector<float> averagePosition; ///< position of this hit combination in world coordinates
    int mainID;                                 ///< Main track ID contributing to this cluster
    int nHits;                                  ///< Number of hits in this cluster
    int mainChannel;                            ///< Main channel of the cluster
    float totalCharge;                          ///< Sum of charges of all associated recob::Hits
    float averageTime;                          ///< Average peak time of all associated recob::Hits
    float averagePurity;                        ///< Purity of the cluster wrt input signal label
    float totalCompleteness;                    ///< Completeness of the cluster wrt main signal particle
    std::vector<recob::Cluster> clusterVector;  ///< Clusters comprising this 3D cluster
  };
} // namespace

#endif //LOWECLUSTER_H
