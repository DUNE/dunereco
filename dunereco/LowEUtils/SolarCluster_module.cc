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
#include "dunereco/LowEUtils/LowECluster.h"

using namespace producer;

namespace solar
{

  class SolarCluster : public art::EDProducer
  {
  public:
    void FillClusters(
      std::vector<LowEUtils::RawSolarCluster> &Clusters,
      std::vector<std::vector<int>> MatchedClustersIdx,
      std::vector<std::vector<int>> ClNHits,
      std::vector<std::vector<int>> ClChannel,
      std::vector<std::vector<float>> ClT,
      std::vector<std::vector<float>> ClY,
      std::vector<std::vector<float>> ClZ,
      std::vector<std::vector<float>> ClCharge,
      std::vector<std::vector<float>> ClPurity,
      std::vector<std::vector<float>> ClCompleteness);

    void ProduceCluster(
      const std::vector<LowEUtils::RawSolarCluster> &RawSolarCluster,
      const std::vector<art::Ptr<recob::Cluster>> &ClusterPtr,
      std::vector<solar::LowECluster> &SolarClusters,
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
    std::string fSignalLabel; // Input tag for MCTruth label
    std::string fGEANTLabel; // Input tag for GEANT label
    std::string fHitLabel; // Input tag for Hit collection
    std::string fClusterLabel; // Input tag for LowECluster collection
    std::string fGeometry;

    float fClusterAlgoTime;
    int fClusterAlgoAdjChannel;
    std::string fClusterChargeVariable;

    int   fClusterMatchNHit;
    float fClusterMatchCharge;
    float fClusterInd0MatchTime;
    float fClusterInd1MatchTime;
    float fClusterMatchTime;
    int   fClusterPreselectionNHits;
    bool  fDebug;
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unique_ptr<producer::ProducerUtils> producer;
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
    fSignalLabel(p.get<std::string>("SignalLabel")),
    fGEANTLabel(p.get<std::string>("GEANT4Label")),
    fHitLabel(p.get<std::string>("HitLabel")),
    fClusterLabel(p.get<std::string>("ClusterLabel")),
    fGeometry(p.get<std::string>("Geometry")),
    fClusterAlgoTime(p.get<float>("ClusterAlgoTime")),
    fClusterAlgoAdjChannel(p.get<int>("ClusterAlgoAdjChannel")),
    fClusterChargeVariable(p.get<std::string>("ClusterChargeVariable")),
    fClusterMatchNHit(p.get<int>("ClusterMatchNHit")),
    fClusterMatchCharge(p.get<float>("ClusterMatchCharge")),
    fClusterInd0MatchTime(p.get<float>("ClusterInd0MatchTime")),
    fClusterInd1MatchTime(p.get<float>("ClusterInd1MatchTime")),
    fClusterMatchTime(p.get<float>("ClusterMatchTime")),
    fClusterPreselectionNHits(p.get<int>("ClusterPreselectionNHits")),
    fDebug(p.get<bool>("Debug")),
    producer(new producer::ProducerUtils(p)),
    lowe(new solar::LowEUtils(p))
  {
    reconfigure(p);
    produces<std::vector<solar::LowECluster>>();
    produces<art::Assns<solar::LowECluster, recob::Cluster>>();
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
    std::string SolarTruthInfo;
    std::string SolarClusterInfo;
    std::string SolarClusterDebug;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Prepare temporary data
    std::vector<std::set<int>> trackids = {};
    std::map<int, simb::MCParticle> ThisGeneratorParts;
    std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {};
    std::vector<std::vector<art::Ptr<recob::Hit>>> Clusters0, Clusters1, Clusters2;
    std::vector<std::vector<int>>  ClustersIdx = {{}, {}, {}}, MatchedClustersIdx = {{}, {}, {}};
    std::vector<std::vector<float>> ClustersPurity = {{}, {}, {}};
    std::vector<std::vector<std::vector<recob::Hit>>> Clusters = {{}, {}, {}}, MatchedClusters = {{}, {}, {}};
    std::vector<std::vector<int>> MatchedClNHits = {{}, {}, {}}, MatchedClChannel = {{}, {}, {}};
    std::vector<std::vector<float>> MatchedClT = {{}, {}, {}}, MatchedClY = {{}, {}, {}}, MatchedClZ = {{}, {}, {}}, MatchedClDir = {{}, {}, {}}; 
    std::vector<std::vector<float>> MatchedClCharge = {{}, {}, {}}, MatchedClPurity = {{}, {}, {}}, MatchedClCompleteness = {{}, {}, {}};
    std::vector<LowEUtils::RawSolarCluster> RawSolarClusters;

    // Prepare the output data
    std::unique_ptr<std::vector<solar::LowECluster>>
        SolarClusterPtr(new std::vector<solar::LowECluster>);

    std::unique_ptr<art::Assns<solar::LowECluster, recob::Cluster>>
        assnPtr(new art::Assns<solar::LowECluster, recob::Cluster>);

    // Get input data
    std::vector<art::Ptr<recob::Cluster>> ClusterPtr;
    art::Handle<std::vector<recob::Cluster>> ClusterHandle;
    
    if (evt.getByLabel(fClusterLabel, ClusterHandle))
    {
      art::fill_ptr_vector(ClusterPtr, ClusterHandle);
      // producer->PrintInColor("Found valid handle for per-plane clusters!", ProducerUtils::GetColor("green"));
      SolarClusterInfo = "Found valid handle for per-plane clusters!";
    }
    else
    {
      // producer->PrintInColor("Valid handle not found!", ProducerUtils::GetColor("red"), "Debug");
      SolarClusterDebug = "Valid handle not found!";
    }


    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Create maps for ID tracking -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    // --- Fill MC Truth IDs to tracking vectors. Get a list of all of my particles in one chunk. ---
    const sim::ParticleList &PartList = pi_serv->ParticleList();
    SolarTruthInfo = "\nThere are a total of " + ProducerUtils::str(PartList.size()) + " particles in the event";

    // Find all simb::MCTruth objects in the event
    // art::ProductID const prodId = evt.getProductID<std::vector<simb::MCTruth>>();
    // art::EDProductGetter const* prodGetter = evt.productGetter(prodId);
    // art::Ptr<std::vector<simb::MCTruth>> MCTruthPtr{ prodId, 0U, prodGetter };

    art::Handle<std::vector<simb::MCTruth>> ThisHandle;
    evt.getByLabel(fSignalLabel, ThisHandle);

    ThisGeneratorParts.clear();
    GeneratorParticles.push_back(ThisGeneratorParts); // For each label insert empty list

    auto SignalHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(fSignalLabel); // Get generator handle
    art::FindManyP<simb::MCParticle> Assn(SignalHandle, evt, fGEANTLabel);            // Assign labels to MCPArticles
    producer->FillMyMaps(GeneratorParticles[0], Assn, SignalHandle);                  // Fill empty list with previously assigned particles
    SolarTruthInfo = SolarTruthInfo + "\n# of particles " + ProducerUtils::str(int(GeneratorParticles[0].size())) + " for " + fSignalLabel;

    const simb::MCNeutrino &nue = SignalHandle->at(0).GetNeutrino();
    SolarTruthInfo = SolarTruthInfo + "\n\tNeutrino Energy: " + ProducerUtils::str(1e3*nue.Nu().E()) + " with vertex (x,y,z) " + ProducerUtils::str(nue.Nu().Vx()) + ", " + ProducerUtils::str(nue.Nu().Vy()) + ", " + ProducerUtils::str(nue.Nu().Vz());

    if (GeneratorParticles[0].size() > 0)
    {
      float eKE = 0;
      std::string ElectronTruthInfo;
      for (std::map<int, simb::MCParticle>::iterator iter = GeneratorParticles[0].begin(); iter != GeneratorParticles[0].end(); iter++)
      {
        std::set<int> ThisGeneratorIDs = {};
        trackids.push_back(ThisGeneratorIDs);
        trackids[0].insert(iter->first);
        int pdg = iter->second.PdgCode();
        if (pdg == 11 && 1e3*iter->second.E()-1e3*iter->second.Mass() > eKE)
        {
          eKE = 1e3*iter->second.E()-1e3*iter->second.Mass();
          ElectronTruthInfo = "\n\tMain e⁻ Energy: " + ProducerUtils::str(eKE) + "MeV with vertex (x,y,z) " + ProducerUtils::str(iter->second.EndX()) + ", " + ProducerUtils::str(iter->second.EndY()) + ", " + ProducerUtils::str(iter->second.EndZ());
        }
      }
      if (eKE > 0)
      {
        SolarTruthInfo = SolarTruthInfo + ElectronTruthInfo;
      }
    }
    else
    {
      std::set<int> ThisGeneratorIDs = {};
      trackids.push_back(ThisGeneratorIDs);
    }

    producer->PrintInColor(SolarTruthInfo, ProducerUtils::GetColor("bright_red"));


    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //-------------------------------------------------------------------- Evaluate MCParticles ---------------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    std::set<int> SignalTrackIDs = {};
    std::vector<int> SignalPDGDepList = {}, SignalIDDepList = {};
    std::vector<float> SignalEDepList = {}, SignalMaxEDepList = {};

    for (size_t i = 0; i < Assn.size(); i++)
    {
      auto SignalParticles = Assn.at(i);
      for (auto SignalParticle = SignalParticles.begin(); SignalParticle != SignalParticles.end(); SignalParticle++)
      {
        std::map<int, float> SignalMaxEDepMap, SignalMaxEDepXMap, SignalMaxEDepYMap, SignalMaxEDepZMap;
        std::vector<const sim::IDE *> ides = bt_serv->TrackIdToSimIDEs_Ps((*SignalParticle)->TrackId());
        float maxEDep = 0;
        for (auto const &ide : ides)
        {
          if (ide->numElectrons < 1 || ide->energy < 1e-6)
          {
            continue;
          } 
          SignalPDGDepList.push_back((*SignalParticle)->PdgCode());
          SignalIDDepList.push_back((*SignalParticle)->TrackId());
          SignalEDepList.push_back(ide->energy);

          if (ide->energy > maxEDep)
          {
            maxEDep = ide->energy;
          }
        } 
        SignalMaxEDepList.push_back(maxEDep);
        SignalTrackIDs.emplace((*SignalParticle)->TrackId());
      }
    }
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------- Evaluate Per-Plane Clusters -----------------------------------------------------------------//
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------//

    int NClusters = ClusterPtr.size();
    // producer->PrintInColor("Analyzing # Clusters " + ProducerUtils::str(NClusters), ProducerUtils::GetColor("green"));
    SolarClusterInfo = SolarClusterInfo + "\nAnalyzing # Clusters " + ProducerUtils::str(NClusters);
    art::FindManyP<recob::Hit> HitAssns(ClusterPtr, evt, fClusterLabel);

    for (int i = 0; i < NClusters; ++i)
    {
      std::vector<recob::Hit> HitVec;
      art::Ptr<recob::Cluster> ThisCluster = ClusterPtr[i];
      std::vector<art::Ptr<recob::Hit>> Hits = HitAssns.at(i);

      // SolarClusterInfo = SolarClusterInfo + "\n\tCluster " + ProducerUtils::str(i) + " has " + ProducerUtils::str(Hits.size()) + " hits"; 
      // SolarClusterInfo = SolarClusterInfo + " with ID " + ProducerUtils::str(int(ThisCluster->ID())) + ", view " + ProducerUtils::str(int(ThisCluster->View())) + " and plane " + ProducerUtils::str(int(ThisCluster->Plane().Plane));

      for (auto const &hit : Hits)
      {
        HitVec.push_back(*hit);
      }

      if (ThisCluster->View() == 0)
      {
        Clusters0.push_back(Hits);
        Clusters[0].push_back(HitVec);
        ClustersIdx[0].push_back(i);
        // std::cout << "Cluster " << i << " is Ind0" << std::endl;
      }
      if (ThisCluster->View() == 1)
      {
        Clusters1.push_back(Hits);
        Clusters[1].push_back(HitVec);
        ClustersIdx[1].push_back(i);
        // std::cout << "Cluster " << i << " is Ind1" << std::endl;
      }
      if (ThisCluster->View() == 2)
      {
        Clusters2.push_back(Hits);
        Clusters[2].push_back(HitVec);
        ClustersIdx[2].push_back(i);
        // std::cout << "Cluster " << i << " is Col." << std::endl;
      }
    }
    SolarClusterInfo = SolarClusterInfo + " (" + ProducerUtils::str(Clusters0.size()) + "," + ProducerUtils::str(Clusters1.size()) + "," + ProducerUtils::str(Clusters2.size()) + ")";
    lowe->MatchClusters(SignalTrackIDs, MatchedClustersIdx, MatchedClusters, ClustersIdx, Clusters, MatchedClNHits, MatchedClChannel, MatchedClT, MatchedClY, MatchedClZ, MatchedClDir, MatchedClCharge, MatchedClPurity, MatchedClCompleteness, clockData, fDebug);
    SolarClusterInfo = SolarClusterInfo + "\nFound " + ProducerUtils::str(int(MatchedClustersIdx[2].size())) + " MatchedClusters (from col. plane loop)!";

    producer->PrintInColor(SolarClusterInfo, ProducerUtils::GetColor("blue"));
    producer->PrintInColor(SolarClusterDebug, ProducerUtils::GetColor("red"), "Debug");
    FillClusters(RawSolarClusters, MatchedClustersIdx, MatchedClNHits, MatchedClChannel, MatchedClT, MatchedClY, MatchedClZ, MatchedClCharge, MatchedClPurity, MatchedClCompleteness);
    ProduceCluster(RawSolarClusters, ClusterPtr, *SolarClusterPtr, clockData);

    // Make the associations which we noted we need
    for (size_t i = 0; i != MatchedClustersIdx[2].size(); ++i)
    {
      art::PtrVector<recob::Cluster> ClusterPtrVector;
      ClusterPtrVector.push_back(art::Ptr<recob::Cluster>(ClusterHandle, MatchedClustersIdx[0][i]));
      ClusterPtrVector.push_back(art::Ptr<recob::Cluster>(ClusterHandle, MatchedClustersIdx[1][i]));
      ClusterPtrVector.push_back(art::Ptr<recob::Cluster>(ClusterHandle, MatchedClustersIdx[2][i]));


      // Create the association between the SolarCluster (solar::LowECluster) and the LowECluster (recob::Cluster)
      util::CreateAssn(*this, evt, *SolarClusterPtr, ClusterPtrVector,
                       *(assnPtr.get()), i);
    }
    
    // Store results into the event
    evt.put(std::move(SolarClusterPtr));
    evt.put(std::move(assnPtr));
  }

  //--------------------------------------------------------------------------
  void SolarCluster::FillClusters(
    std::vector<LowEUtils::RawSolarCluster> &Clusters,
    std::vector<std::vector<int>> MatchedClustersIdx,
    std::vector<std::vector<int>> ClNHits,
    std::vector<std::vector<int>> ClChannel,
    std::vector<std::vector<float>> ClT,
    std::vector<std::vector<float>> ClY,
    std::vector<std::vector<float>> ClZ,
    std::vector<std::vector<float>> ClCharge,
    std::vector<std::vector<float>> ClPurity,
    std::vector<std::vector<float>> ClCompleteness
    )
  {
    std::string FillClustersDebug;
    FillClustersDebug = "Filling " + ProducerUtils::str(int(ClNHits[2].size())) + " SolarClusters";
    for (int i = 0; i < int(ClNHits[2].size()); i++)
    {
      std::vector<float> Position = {0, float(ClY[2][i]), float(ClZ[2][i])};
      std::vector<int> ClusterIdxVec = {};
      for (size_t ii = 0; ii < 3; ii++)
      {
        ClusterIdxVec.push_back(MatchedClustersIdx[ii][i]);
      }
      Clusters.push_back(LowEUtils::RawSolarCluster{Position, int(ClusterIdxVec.size()), ClChannel[2][1], ClCharge[2][i], ClT[2][i], ClPurity[2][i], ClCompleteness[2][i], ClusterIdxVec});
    }
    producer->PrintInColor(FillClustersDebug, ProducerUtils::GetColor("green"), "Debug");
  }

  void SolarCluster::ProduceCluster(const std::vector<LowEUtils::RawSolarCluster> &Clusters,
                                    const std::vector<art::Ptr<recob::Cluster>> &ClusterPtr,
                                    std::vector<solar::LowECluster> &SolarClusters,
                                    detinfo::DetectorClocksData const &ts) const
  {
    std::string ProduceClusterInfo;
    std::string ProduceClusterDebug;
    ProduceClusterDebug = "Producing " + ProducerUtils::str(int(Clusters.size())) + " SolarClusters";

    std::vector<recob::Cluster> ClusterVec(3);
    for (int i = 0; i < int(Clusters.size()); i++)
    {
      // std::cout << "Producing SolarCluster " << i << std::endl;
      LowEUtils::RawSolarCluster Cluster = Clusters[i];
      std::vector<float> Position = Cluster.Position;
      int NHits = Cluster.NHits;
      int MainChannel = Cluster.MainChannel;
      float TotalCharge = Cluster.TotalCharge;
      float AveragePeakTime = Cluster.AveragePeakTime;
      float Purity = Cluster.Purity;
      float Completeness = Cluster.Completeness;
      std::vector<int> ClusterIdxVec = Cluster.ClusterIdxVec;

      ClusterVec.clear();
      for (int j = 0; j < 3; j++)
      {
        if (ClusterIdxVec[j] < 0)
        {
          ClusterVec.push_back(recob::Cluster());
        }
        else
        {
          // Retrieve the cluster from the input pointer vector
          recob::Cluster ThisCluster = *ClusterPtr[ClusterIdxVec[j]];
          ClusterVec.push_back(*ClusterPtr[ClusterIdxVec[j]]);
          ProduceClusterDebug = ProduceClusterDebug + "\nSolarCluster " + ProducerUtils::str(i) + " has " + ProducerUtils::str(j) + " cluster ID " + ProducerUtils::str(int(ThisCluster.ID())) + " and cluster Idx " + ProducerUtils::str(ClusterIdxVec[j]);
        }
      }
      if ( i < 10)
      {
        ProduceClusterInfo = ProduceClusterInfo + "\nSolarCluster " + ProducerUtils::str(i) + " has position (t,y,z) " + ProducerUtils::str(AveragePeakTime) + ", " + ProducerUtils::str(Position[1]) + ", " + ProducerUtils::str(Position[2]);
        ProduceClusterInfo = ProduceClusterInfo + " with total charge " + ProducerUtils::str(TotalCharge) + ", purity " + ProducerUtils::str(Purity) + " and completeness " + ProducerUtils::str(Completeness);
      }
      SolarClusters.emplace_back(Position, NHits, MainChannel, TotalCharge, AveragePeakTime, Purity, Completeness, ClusterVec);
    }
    producer->PrintInColor(ProduceClusterInfo, ProducerUtils::GetColor("green"));
    producer->PrintInColor(ProduceClusterDebug, ProducerUtils::GetColor("green"), "Debug");
  }
} // namespace solar
