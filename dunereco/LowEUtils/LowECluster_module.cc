////////////////////////////////////////////////////////////////////////////////////
// Class:       LowECluster                                                       //
// Module Type: producer                                                          //
// File:        LowECluster_module.cc                                             //
//                                                                                //
// Hit clustering, based on time and proximity information.                       //
// To be used as per-plane clustering for the SolarCluser_module.                 //
////////////////////////////////////////////////////////////////////////////////////

#ifndef LowECluster_H
#define LowECluster_H 1

// #include "duneopdet/SolarNuUtils/SolarAuxUtils.h"
#include "dunereco/LowEUtils/LowEUtils.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// C++ Includes
#include <vector>
#include <string>
#include <memory>
#include <limits>

namespace solar
{

  class LowECluster : public art::EDProducer
  {
  public:
    void ProduceCluster(const std::vector<LowEUtils::LowEClusterInfo> &LowEClusterinfo,
                        std::vector<recob::Cluster> &LowEClusters,
                        detinfo::DetectorClocksData const &ts) const;

    // Standard constructor and destructor for an ART module.
    explicit LowECluster(const fhicl::ParameterSet &);
    virtual ~LowECluster();

    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const &pset);

    // The producer routine, called once per event.
    void produce(art::Event &);

  private:
    // The parameters we'll read from the .fcl file.
    art::ServiceHandle<geo::Geometry> geo;
    std::string fHitLabel; // Input tag for Hit collection
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
    // std::unique_ptr<solar::SolarAuxUtils> solaraux;
    std::unique_ptr<solar::LowEUtils> lowe;
  };

}

namespace solar
{
  DEFINE_ART_MODULE(LowECluster)
}

#endif

namespace solar
{

  //--------------------------------------------------------------------------
  // Constructor
  LowECluster::LowECluster(const fhicl::ParameterSet &p)
      : EDProducer{p},
        fHitLabel(p.get<std::string>("HitLabel")),
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
        // solaraux(new solar::SolarAuxUtils(p)),
        lowe(new solar::LowEUtils(p))
  {
    reconfigure(p);
    produces<std::vector<recob::Cluster>>();
    produces<art::Assns<recob::Cluster, recob::Hit>>();
  }

  //--------------------------------------------------------------------------
  void LowECluster::reconfigure(fhicl::ParameterSet const &p)
  {
  }

  //--------------------------------------------------------------------------
  // Destructor
  LowECluster::~LowECluster()
  {
  }
  //--------------------------------------------------------------------------
  void LowECluster::beginJob()
  {
  }

  //--------------------------------------------------------------------------
  void LowECluster::endJob()
  {
  }

  //--------------------------------------------------------------------------
  void LowECluster::produce(art::Event &evt)
  {
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    // Prepare the input data
    std::vector<std::vector<int>> HitIdx;
    std::vector<art::Ptr<recob::Hit>> Hits;
    std::vector<std::vector<art::Ptr<recob::Hit>>> Clusters;
    std::vector<LowEUtils::LowEClusterInfo> LowEClusters;

    // Prepare the output data
    std::unique_ptr<std::vector<recob::Cluster>>
        ClusterPtr(new std::vector<recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>>
        assnPtr(new art::Assns<recob::Cluster, recob::Hit>);

    // Get Hits from the event
    auto HitHandle = evt.getHandle<std::vector<recob::Hit>>(fHitLabel);
    int NTotHits = HitHandle->size();
    for (int i = 0; i < NTotHits; ++i)
    {
      art::Ptr<recob::Hit> TPCHit(HitHandle, i);
      Hits.push_back(TPCHit);
    }

    // Run the clustering
    lowe->CalcAdjHits(Hits, Clusters, HitIdx);
    lowe->MakeClusterVector(LowEClusters, Clusters, evt);
    ProduceCluster(LowEClusters, *ClusterPtr, clockData);

    // Make the associations which we noted we need
    for (size_t i = 0; i != HitIdx.size(); ++i)
    {
      // std::cout << "Cluster " << i << " has " << HitIdx.at(i).size() << " hits" << std::endl;
      art::PtrVector<recob::Hit> HitPtrVector;
      for (int const &hitIndex : HitIdx.at(i))
      {
        art::Ptr<recob::Hit> HitPtr(HitHandle, hitIndex);
        HitPtrVector.push_back(HitPtr);
      }
      if (i < 10)
      {
        std::cout << "Genrating Cluster " << i << " with " << HitPtrVector.size() << " hits" << std::endl;
      }
      // Create the association between the flash and the Hits
      util::CreateAssn(*this, evt, *ClusterPtr, HitPtrVector,
                       *(assnPtr.get()), i);
      if (i == HitIdx.size() - 1)
      {
        std::cout << "..." << std::endl;
        std::cout << "Generated " << i + 1 << " Cluster" << std::endl;
      }
    }
    // Store results into the event
    evt.put(std::move(ClusterPtr));
    evt.put(std::move(assnPtr));
  }

  //--------------------------------------------------------------------------
  void LowECluster::ProduceCluster(const std::vector<LowEUtils::LowEClusterInfo> &Clusters,
                                   std::vector<recob::Cluster> &LowEClusters,
                                   detinfo::DetectorClocksData const &ts) const
  {
    // Loop over the flashes with TheFlash being the flash
    for (int i = 0; i < int(Clusters.size()); i++)
    {
      LowEUtils::LowEClusterInfo Cluster = Clusters[i];

      // Time to make the Cluster collect the info from the Clusterinfo struct
      float StartWire = Cluster.StartWire;
      float SigmaStartWire = Cluster.SigmaStartWire;
      float StartTick = Cluster.StartTick;
      float SigmaStartTick = Cluster.SigmaStartTick;
      float StartCharge = Cluster.StartCharge;
      float StartAngle = Cluster.StartAngle;
      float StartOpeningAngle = Cluster.StartOpeningAngle;
      float EndWire = Cluster.EndWire;
      float SigmaEndWire = Cluster.SigmaEndWire;
      float EndTick = Cluster.EndTick;
      float SigmaEndTick = Cluster.SigmaEndTick;
      float EndCharge = Cluster.EndCharge;
      float EndAngle = Cluster.EndAngle;
      float EndOpeningAngle = Cluster.EndOpeningAngle;
      float Integral = Cluster.Integral;
      float IntegralStdDev = Cluster.IntegralStdDev;
      float SummedADC = Cluster.SummedADC;
      float SummedADCstdDev = Cluster.SummedADCstdDev;
      int NHit = Cluster.NHit;
      float MultipleHitDensity = Cluster.MultipleHitDensity;
      float Width = Cluster.Width;
      int ID = Cluster.ID;
      geo::View_t View = Cluster.View;
      geo::PlaneID Plane = Cluster.Plane;

      LowEClusters.emplace_back(StartWire, SigmaStartWire, StartTick, SigmaStartTick, StartCharge, StartAngle, StartOpeningAngle,
                                EndWire, SigmaEndWire, EndTick, SigmaEndTick, EndCharge, EndAngle, EndOpeningAngle,
                                Integral, IntegralStdDev, SummedADC, SummedADCstdDev,
                                NHit, MultipleHitDensity, Width, ID, View, Plane);
    }
  }
} // namespace solar
