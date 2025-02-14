////////////////////////////////////////////////////////////////////////////////////
// Class:       SolarCluster                                                      //
// Module Type: producer                                                          //
// File:        SolarCluster_module.cc                                            //
//                                                                                //
// Hit clustering, based on plane-matching information.                           //
// To be used for 2D cluster reconstrution SolarNuAna_module & LowEAna_module.    //
////////////////////////////////////////////////////////////////////////////////////

#ifndef SolarCluster_H
#define SolarCluster_H 1

#include "dunereco/LowEUtils/LowEUtils.h"

namespace solar
{

  class SolarCluster : public art::EDProducer
  {
  public:
    void ProduceCluster(
      const std::vector<LowEUtils::RawSolarCluster> &RawSolarCluster,
      std::vector<reco::ClusterHit3D> &SolarClusters,
      detinfo::DetectorClocksData const &ts) const;

    // Standard constructor and destructor for an ART module.
    explicit SolarCluster(const fhicl::ParameterSet &);
    virtual ~SolarCluster();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const &pset);

    // The producer routine, called once per event.
    void produce(art::Event &);

  private:
    // The parameters we'll read from the .fcl file.
    art::ServiceHandle<geo::Geometry> geo;
    std::string fHitLabel; // Input tag for Hit collection
    std::string fClusterLabel; // Input tag for LowECluster collection
    std::string fGeometry;
    double fDetectorSizeX;

    float fClusterAlgoTime;
    int fClusterAlgoAdjChannel;
    std::string fClusterChargeVariable;

    int fClusterMatchNHit;
    float fClusterMatchCharge;
    float fClusterInd0MatchTime;
    float fClusterInd1MatchTime;
    float fClusterMatchTime;
    int fClusterPreselectionNHits;
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    // std::unique_ptr<solar::SolarAuxUtils> solaraux;
    std::unique_ptr<solar::LowEUtils> lowe;
  };
}

namespace solar
{
  DEFINE_ART_MODULE(SolarCluster)
}

#endif

namespace solar
{
  //--------------------------------------------------------------------------
  // Constructor
  SolarCluster::SolarCluster(const fhicl::ParameterSet &p)
  : EDProducer{p},
    fHitLabel(p.get<std::string>("HitLabel")),
    fClusterLabel(p.get<std::string>("ClusterLabel")),
    fGeometry(p.get<std::string>("Geometry")),
    fDetectorSizeX(p.get<double>("DetectorSizeX")),
    fClusterAlgoTime(p.get<float>("ClusterAlgoTime")),
    fClusterAlgoAdjChannel(p.get<int>("ClusterAlgoAdjChannel")),
    fClusterChargeVariable(p.get<std::string>("ClusterChargeVariable")),
    fClusterMatchNHit(p.get<int>("ClusterMatchNHit")),
    fClusterMatchCharge(p.get<float>("ClusterMatchCharge")),
    fClusterInd0MatchTime(p.get<float>("ClusterInd0MatchTime")),
    fClusterInd1MatchTime(p.get<float>("ClusterInd1MatchTime")),
    fClusterMatchTime(p.get<float>("ClusterMatchTime")),
    fClusterPreselectionNHits(p.get<int>("ClusterPreselectionNHits")),
    // solaraux(new solar::SolarAuxUtils(p)),
    lowe(new solar::LowEUtils(p))
  {
    reconfigure(p);
    produces<std::vector<reco::ClusterHit3D>>();
    produces<art::Assns<reco::ClusterHit3D, recob::Cluster>>();
  }

  //--------------------------------------------------------------------------
  void SolarCluster::reconfigure(fhicl::ParameterSet const &p)
  {
  }

  //--------------------------------------------------------------------------
  // Destructor
  SolarCluster::~SolarCluster()
  {
  }
  //--------------------------------------------------------------------------
  void SolarCluster::beginJob()
  {
  }

  //--------------------------------------------------------------------------
  void SolarCluster::endJob()
  {
  }

  //--------------------------------------------------------------------------
  void SolarCluster::produce(art::Event &evt)
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Prepare temporary data
    std::vector<std::vector<art::Ptr<recob::Hit>>> Clusters0, Clusters1, Clusters2;
    std::vector<std::vector<int>>  ClustersIdx = {{}, {}, {}}, MatchedClustersIdx = {{}, {}, {}};
    std::vector<std::vector<std::vector<recob::Hit>>> Clusters = {{}, {}, {}}, MatchedClusters = {{}, {}, {}};
    std::vector<std::vector<int>> ClNHits = {{}, {}, {}};
    std::vector<std::vector<float>> ClT = {{}, {}, {}}, ClCharge = {{}, {}, {}};
    std::vector<LowEUtils::RawSolarCluster> RawSolarClusters;

    // Prepare the output data
    std::unique_ptr<std::vector<reco::ClusterHit3D>>
        SolarClusterPtr(new std::vector<reco::ClusterHit3D>);

    std::unique_ptr<art::Assns<reco::ClusterHit3D, recob::Cluster>>
        assnPtr(new art::Assns<reco::ClusterHit3D, recob::Cluster>);

    // Get input data
    std::vector<art::Ptr<recob::Cluster>> ClusterPtr;
    art::Handle<std::vector<recob::Cluster>> ClusterHandle;
    
    if (evt.getByLabel(fClusterLabel, ClusterHandle))
    {
      art::fill_ptr_vector(ClusterPtr, ClusterHandle);
    }
    int NClusters = ClusterPtr.size();
    art::FindManyP<recob::Hit> HitAssns(ClusterPtr, evt, fClusterLabel);
    
    for (int i = 0; i < NClusters; ++i)
    {
      std::vector<art::Ptr<recob::Hit>> Hits = HitAssns.at(i);

      std::vector<recob::Hit> HitVec;

      int View = 0;
      int SignalType = 0;
      int NHits = Hits.size();

      if (NHits > 0)
      {
        for (auto const &hit : Hits)
        {
          HitVec.push_back(*hit);
        }
        View = HitVec[0].View();
        SignalType = HitVec[0].SignalType();
      }
      else
      {
        std::cout << "No hits in cluster " << i << std::endl;
        continue;
      }
      if (View == 0)
      {
        if (SignalType == 0)
        {
          Clusters0.push_back(Hits);
          Clusters[0].push_back(HitVec);
          ClustersIdx[0].push_back(i);
        }
        else
        {
          Clusters1.push_back(Hits);
          Clusters[1].push_back(HitVec);
          ClustersIdx[1].push_back(i);
        }
      }
      else
      {
        Clusters2.push_back(Hits);
        Clusters[2].push_back(HitVec);
        ClustersIdx[2].push_back(i);
      }
    }

    lowe->MatchClusters(MatchedClustersIdx, MatchedClusters, ClustersIdx, Clusters, ClNHits, ClT, ClCharge, false);
    std::vector<std::vector<float>> ClY = ClT;
    std::vector<std::vector<float>> ClZ = ClT;
    // Reset all vector entries in ClY and ClZ to -1e6
    for (size_t i = 0; i < ClY.size(); i++)
    {
      for (size_t j = 0; j < ClY[i].size(); j++)
      {
        ClY[i][j] = -1e6;
        ClZ[i][j] = -1e6;
      }
    }
    lowe->ComputeCluster3D(RawSolarClusters, MatchedClusters, clockData);
    ProduceCluster(RawSolarClusters, *SolarClusterPtr, clockData);

    // Make the associations which we noted we need
    for (size_t i = 0; i != MatchedClustersIdx[2].size(); ++i)
    {
      art::PtrVector<recob::Cluster> ClusterPtrVector;
      ClusterPtrVector.push_back(art::Ptr<recob::Cluster>(ClusterHandle, MatchedClustersIdx[0][i]));
      ClusterPtrVector.push_back(art::Ptr<recob::Cluster>(ClusterHandle, MatchedClustersIdx[1][i]));
      ClusterPtrVector.push_back(art::Ptr<recob::Cluster>(ClusterHandle, MatchedClustersIdx[2][i]));


      // Create the association between the SolarCluster (reco::ClusterHit3D) and the LowECluster (recob::Cluster)
      util::CreateAssn(*this, evt, *SolarClusterPtr, ClusterPtrVector,
                       *(assnPtr.get()), i);
    }
    
    // Store results into the event
    evt.put(std::move(SolarClusterPtr));
    evt.put(std::move(assnPtr));
  }

  //--------------------------------------------------------------------------
  void SolarCluster::ProduceCluster(const std::vector<LowEUtils::RawSolarCluster> &Clusters,
                                   std::vector<reco::ClusterHit3D> &SolarClusters,
                                   detinfo::DetectorClocksData const &ts) const
  {
    // Loop over the flashes with TheFlash being the flash
    for (int i = 0; i < int(Clusters.size()); i++)
    {
      LowEUtils::RawSolarCluster Cluster = Clusters[i];
      recob::Cluster::ID_t ID = Cluster.ID;
      unsigned int StatusBits = Cluster.StatusBits;
      Eigen::Vector3f Position = Cluster.Position;
      float TotalCharge = Cluster.TotalCharge;
      float AveragePeakTime = Cluster.AveragePeakTime;
      float DeltaPeakTime = Cluster.DeltaPeakTime;
      float SigmaPeakTime = Cluster.SigmaPeakTime;
      float HitChiSquare = Cluster.HitChiSquare;
      float OverlapFraction = Cluster.OverlapFraction;
      float ChargeAsymmetry = Cluster.ChargeAsymmetry;
      float DOCAToAxis = Cluster.DOCAToAxis;
      float ArcLenToPOCA = Cluster.ArcLenToPOCA;
      reco::ClusterHit2DVec HitVec = Cluster.HitVec;
      std::vector<float> HitDeltaTSigmaVec = Cluster.HitDeltaTSigmaVec;
      std::vector<geo::WireID> WireIDVec = Cluster.WireIDVec;
      // size_t id, unsigned int statusBits, const Eigen::Vector3f &position, float totalCharge, float avePeakTime, float deltaPeakTime, float sigmaPeakTime, float hitChiSquare, float overlapFraction, float chargeAsymmetry, float docaToAxis, float arclenToPoca, const ClusterHit2DVec &hitVec, const std::vector< float > &hitDelTSigVec, const std::vector< geo::WireID > &wireIDVec
      SolarClusters.emplace_back(ID, StatusBits, Position, TotalCharge, AveragePeakTime, DeltaPeakTime, SigmaPeakTime, HitChiSquare, OverlapFraction, ChargeAsymmetry, DOCAToAxis, ArcLenToPOCA, HitVec, HitDeltaTSigmaVec, WireIDVec);
    }
  }
} // namespace solar
