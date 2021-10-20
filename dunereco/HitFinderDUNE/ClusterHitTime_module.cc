////////////////////////////////////////////////////////////////////////
// Class:       ClusterHitTime
// Plugin Type: analyzer (art v3_05_00)
// File:        ClusterHitTime_module.cc
//
// Generated at Thu Apr  2 08:17:05 2020 by Tingjun Yang using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art_root_io/TFileService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/CoreUtils/ServiceUtil.h"

#include "TH1D.h"

class ClusterHitTime;


class ClusterHitTime : public art::EDAnalyzer {
public:
  explicit ClusterHitTime(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ClusterHitTime(ClusterHitTime const&) = delete;
  ClusterHitTime(ClusterHitTime&&) = delete;
  ClusterHitTime& operator=(ClusterHitTime const&) = delete;
  ClusterHitTime& operator=(ClusterHitTime&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:
  
  TH1D *clutime[3];
  art::InputTag fClusterModuleLabel;
  std::vector<int> fClusterIDs;
};


ClusterHitTime::ClusterHitTime(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fClusterModuleLabel(p.get<art::InputTag>("ClusterModuleLabel","linecluster")),
  fClusterIDs(p.get<std::vector<int>>("ClusterIDs"))
{
}

void ClusterHitTime::analyze(art::Event const& e)
{

  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  auto clusterListHandle = e.getHandle< std::vector<recob::Cluster> >(fClusterModuleLabel);
  if (clusterListHandle)
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::FindManyP<recob::Hit> fmhc(clusterListHandle, e, fClusterModuleLabel);

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);
  for (const auto & cluster : clusterlist){
    bool match = false;
    for (const auto & id : fClusterIDs){
      if (int(cluster.key()) == id){
        match = true;
      }
    }
    //std::cout<<cluster.key()<<std::endl;
    if (!match) continue;
    auto &hits = fmhc.at(cluster.key());
    for (const auto & hit : hits){
      std::cout<<hit->WireID().Plane<<" "<<hit->PeakTime()<<" "<<detProp.GetXTicksOffset(hit->WireID())<<" "<<hit->Integral()<<std::endl;
      clutime[hit->WireID().Plane]->Fill(hit->PeakTime()-detProp.GetXTicksOffset(hit->WireID())+5000, hit->Integral());
    }
  }

}

void ClusterHitTime::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  clutime[0] = tfs->make<TH1D>("clutime_0","Plane 0;Tick;ADC",2000,0,2000);
  clutime[1] = tfs->make<TH1D>("clutime_1","Plane 1;Tick;ADC",2000,0,2000);
  clutime[2] = tfs->make<TH1D>("clutime_2","Plane 2;Tick;ADC",2000,0,2000);

}

DEFINE_ART_MODULE(ClusterHitTime)
