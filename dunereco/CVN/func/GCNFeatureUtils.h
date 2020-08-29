////////////////////////////////////////////////////////////////////////
/// \file    GCNFeatureUtils.h
/// \brief   Utilities for calculating feature values for the GCN
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef GCN_FEATURE_UTILS_H
#define GCN_FEATURE_UTILS_H

#include <vector>
#include <tuple>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcorealg/GeoAlgo/GeoAlgo.h"
namespace detinfo {
  class DetectorClocksData;
}

#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/PixelMap.h"

namespace cvn
{
  typedef std::tuple<int, int, int, float, float, float, float, float, float,
    float, std::string, std::string> ptruth;

  /// Class containing some utility functions for all things CVN
  class GCNFeatureUtils
  {
  public:
    GCNFeatureUtils();
    ~GCNFeatureUtils();

    /// Get primary true G4 ID for hit
    int GetTrackIDFromHit(detinfo::DetectorClocksData const& clockData, const recob::Hit&) const;
    /// Get closest trajectory point for true particle given reconstructed spacepoint
    std::pair<unsigned int, float> GetClosestApproach(const recob::SpacePoint sp, const simb::MCParticle p) const;

    /// Get the number of neighbours within rangeCut cm of this space point
    unsigned int GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, const float rangeCut, const std::string &spLabel) const;
    unsigned int GetSpacePointNeighbours(const recob::SpacePoint &sp, art::Event const &evt, const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const;
    /// Get a map of the number of neighbours for each space point ID. Using this function is much less wasteful
    /// than repeated calls to the above function
    std::map<int,unsigned int> GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::string &spLabel) const;
    std::map<int,unsigned int> GetAllNeighbours(art::Event const &evt, const float rangeCut, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const;

    /// Gets the number of nearest neigbours for each space point for a vector of cut values. Much more efficient that using the above functions multiple times.
    std::vector<std::map<int,unsigned int>> GetNeighboursForRadii(art::Event const &evt, const std::vector<float>& rangeCuts, const std::string &spLabel) const;
    std::vector<std::map<int,unsigned int>> GetNeighboursForRadii(art::Event const &evt, const std::vector<float>& rangeCuts, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const;
    std::vector<std::map<int,unsigned int>> GetNeighboursForRadii(art::Event const &evt, const std::vector<float>& rangeCuts, const std::map<unsigned int,art::Ptr<recob::SpacePoint>> &sps) const;

    /// Get the nearest neighbour map for the spacepoints. Returns a map of <nose_spacepoint_id,nearest_neighbour_spacepoint_id>
    std::map<int,int> GetNearestNeighbours(art::Event const &evt, const std::string &spLabel) const;
    std::map<int,int> GetNearestNeighbours(art::Event const &evt, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const;

    /// Get the two nearest neighbours to use for calcuation of angles between them and the node in question
    std::map<int,std::pair<int,int>> GetTwoNearestNeighbours(art::Event const &evt, const std::string &spLabel) const;
    std::map<int,std::pair<int,int>> GetTwoNearestNeighbours(art::Event const &evt, const std::vector<art::Ptr<recob::SpacePoint>> &sps) const;
    std::map<int,std::pair<int,int>> GetTwoNearestNeighbours(art::Event const &evt, const std::map<unsigned int,art::Ptr<recob::SpacePoint>> &sps) const;

    /// Get the angle and the dot product between the vector from the base node to its neighbours
    void GetAngleAndDotProduct(const recob::SpacePoint &baseNode, const recob::SpacePoint &n1, const recob::SpacePoint &n2, float &dotProduct, float &angle) const;

    /// Use the association between space points and hits to return a charge
    std::map<unsigned int, float> GetSpacePointChargeMap(std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
                                                         std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const;
    std::map<unsigned int, float> GetSpacePointChargeMap(art::Event const &evt, const std::string &spLabel) const;

    // Mean RMS of the hits making up a spacepoint
    std::map<unsigned int, float> GetSpacePointMeanHitRMSMap(art::Event const &evt, const std::string &spLabel) const;

    /// Get the true G4 ID for each spacepoint using energy matching
    std::map<unsigned int, int> GetTrueG4ID(detinfo::DetectorClocksData const& clockData,
                                            std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
                                            std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const;
    std::map<unsigned int, int> GetTrueG4ID(detinfo::DetectorClocksData const& clockData,
                                            art::Event const& evt, const std::string &spLabel) const;
    /// Get the true G4 ID for each spacepoint using energy matching
    std::map<unsigned int, int> GetTrueG4IDFromHits(
      detinfo::DetectorClocksData const& clockData,
      std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
      std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const;
    std::map<unsigned int, int> GetTrueG4IDFromHits(
      detinfo::DetectorClocksData const& clockData,
      art::Event const& evt, const std::string &spLabel) const;

    /// Get the true pdg code for each spacepoint
    std::map<unsigned int, int> GetTruePDG(
      detinfo::DetectorClocksData const& clockData,
      art::Event const& evt, const std::string &spLabel, bool useAbsoluteTrackID, bool useHits) const;
    /// Get 2D hit features for a given spacepoint
    std::map<unsigned int, std::vector<float>> Get2DFeatures(
      std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
      std::vector<std::vector<art::Ptr<recob::Hit>>> const& sp2Hit) const;

    /// Convert a pixel map into three 2D GCNGraph objects
    std::vector<cvn::GCNGraph> ExtractGraphsFromPixelMap(const cvn::PixelMap &pm, const float chargeThreshold) const;
    /// Get the neighbours map <graph node, neighbours> for the three 2D graph in 2 box (npixel+1) around the pixel
    std::map<unsigned int,unsigned int> Get2DGraphNeighbourMap(const cvn::GCNGraph &g, const unsigned int npixel) const;

    /// Get ground truth for spacepoint deghosting graph network
    std::vector<float> GetNodeGroundTruth(detinfo::DetectorClocksData const& clockData,
                                          std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints,
                                          std::vector<std::vector<art::Ptr<recob::Hit>>> const& spToHit,
                                          float distCut,
                                          std::vector<std::vector<float>>* dirTruth=nullptr) const;
    /// Get hierarchy map from set of particles
    std::map<unsigned int, unsigned int> GetParticleFlowMap(const std::set<unsigned int>& particles) const;

    static std::vector<ptruth> GetParticleTree(const cvn::GCNGraph* g);

  private:

    geoalgo::GeoAlgo fGeoAlgo;

  };

}

#endif  // GCN_FEATURE_UTILS_H
