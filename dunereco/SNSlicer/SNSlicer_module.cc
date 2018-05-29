
// Tutorial module example

#ifndef SNSlicer_H
#define SNSlicer_H 1

// Framework includes

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
//#include "larsim/MCCheater/BackTrackerService.h"
//#include "larsim/MCCheater/ParticleInventoryService.h"
//#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "nutools/ParticleNavigation/EmEveIdCalculator.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// ROOT includes

#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"
#include "TPad.h"

// C++ includes

#include <vector>
#include <string>
#include <iostream>

#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Wire.h"


#include "dunetpc/dune/SNSlicer/SNSlice.h"

namespace
{
  struct Pt2D
  {
    Pt2D(){}
    Pt2D(double _x, double _y, double _q, bool _isSig, int _idx) : x(_x), y(_y), q(_q), isSig(_isSig), idx(_idx) {}
    double x, y;
    double q;
    bool isSig;
    int idx; // index into original hits array
  };

  class SNSlicer: public art::EDProducer 
  {
  public:
    SNSlicer(const fhicl::ParameterSet&);

    void produce(art::Event&) override;
   
  protected:
    std::vector<int> regionQuery(const std::vector<Pt2D>& D, const Pt2D& p) const;

    void expandCluster(const std::vector<Pt2D>& D,
                       int ip,
                       std::vector<int> neiPts, int C,
                       std::vector<bool>& visited,
                       std::vector<int>& inClust) const;

    std::vector<std::vector<Pt2D>> DBSCAN(const std::vector<Pt2D>& D) const;

    //    std::string fDetSimProducerLabel;

//    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
//    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
//    art::ServiceHandle<cheat::PhotonBackTrackerService> pbt;

    const geo::GeometryCore& geom;
    const detinfo::DetectorProperties& dp;

    // DBScan params
    double fEps;
    int fMinPts;
  };

}

#endif 

namespace {

  DEFINE_ART_MODULE(SNSlicer)

}

namespace {

  //---------------------------------------------------------------------------
  // Constructor
  SNSlicer::SNSlicer(const fhicl::ParameterSet& pset)
    : EDProducer(),//pset),
      geom(*lar::providerFrom<geo::Geometry>()),
      dp(*lar::providerFrom<detinfo::DetectorPropertiesService>())
  {
    produces<std::vector<sn::SNSlice>>();

    // Read the fcl-file
    //    fDetSimProducerLabel = pset.get< std::string >("DetSimLabel");

    fEps = pset.get<double>("Eps");
    fMinPts = pset.get<int>("MinPts");
  }

  //---------------------------------------------------------------------------
  std::vector<int> SNSlicer::regionQuery(const std::vector<Pt2D>& D, const Pt2D& p) const
  {
    std::vector<int> ret;
    for(unsigned int i = 0; i < D.size(); ++i){
      if((D[i].x-p.x)*(D[i].x-p.x) + (D[i].y-p.y)*(D[i].y-p.y) < fEps*fEps)
        ret.push_back(i);
    }
    return ret;
  }

  //---------------------------------------------------------------------------
  std::vector<int> Join(const std::vector<int> a, const std::vector<int> b)
  {
    std::set<int> s(a.begin(), a.end());
    s.insert(b.begin(), b.end());

    return std::vector<int>(s.begin(), s.end());
  }

  //---------------------------------------------------------------------------
  void SNSlicer::expandCluster(const std::vector<Pt2D>& D,
                               int ip,
                               std::vector<int> neiPts, int C,
                               std::vector<bool>& visited,
                               std::vector<int>& inClust) const
  {
    inClust[ip] = C;

    bool any = true;
    while(any){
      any = false;
      for(unsigned int iNP = 0; iNP < neiPts.size(); ++iNP){ // P'
        if(inClust[neiPts[iNP]] == -1){
          inClust[neiPts[iNP]] = C;
        }

        if(!visited[neiPts[iNP]]){
          visited[neiPts[iNP]] = true;
          std::vector<int> neiPts2 = regionQuery(D, D[neiPts[iNP]]);
          //        double Q = 0;
          //        for(int i: neiPts2) Q += D[i].q;
          const int Q = neiPts2.size();
          if(Q >= fMinPts){
            neiPts = Join(neiPts, neiPts2);
            any = true;
            break; // restart the while loop
          }
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  std::vector<std::vector<Pt2D>> SNSlicer::DBSCAN(const std::vector<Pt2D>& D) const
  {
    std::vector<bool> visited(D.size(), false);
    //    std::vector<bool> noise(D.size(), false);
    std::vector<int> inclust(D.size(), -1);

    int C = 0; // TODO, should this be -1 since we ++ below?
    for(unsigned int iD = 0; iD < D.size(); ++iD){
      if(visited[iD]) continue;
      visited[iD] = true;

      std::vector<int> neiPts = regionQuery(D, D[iD]);
      //      double Q = 0;
      //      for(int i: neiPts) Q += /D[i].q;
      const int Q = neiPts.size();
      if(Q < fMinPts){
        //        noise[iD] = true;
      }
      else {
        ++C;
        expandCluster(D, iD, neiPts, C, visited, inclust);
      }
    }

    std::vector<std::vector<Pt2D>> ret(C+1);
    for(unsigned int i = 0; i < D.size(); ++i){
      if(inclust[i] >= 0) ret[inclust[i]].push_back(D[i]);
    }
    return ret;
  }

  //---------------------------------------------------------------------------
  void SNSlicer::produce(art::Event& evt)
  {
    std::unique_ptr<std::vector<sn::SNSlice>> slicecol(new std::vector<sn::SNSlice>);

    //    sim::EveIdCalculator eve;
    //    eve.Init(&bt->ParticleList());

    art::Handle<std::vector<recob::Hit>> hits;
    evt.getByLabel("gaushit", hits);
    //    std::cout << "nhits: " << hits->size() << std::endl;

    // art::Handle<std::vector<sim::SimChannel>> simchans;
    // evt.getByLabel("largeant", simchans);
    // std::cout << "nsimchans: " << simchans->size() << std::endl;

    //    art::Handle<std::vector<recob::OpHit>> ophits;
    //    evt.getByLabel("ophit", ophits);

    // std::set<raw::ChannelID_t> hitChans;

    // for(unsigned int hitIdx = 0; hitIdx < hits->size(); ++hitIdx){
    //   const recob::Hit& hit = (*hits)[hitIdx];

    //   if(hit.View() != geo::kZ) continue; // vertical ie collection

    //   std::vector<sim::IDE> ides;
    //   bt->HitToSimIDEs(hit, ides);
    //   for(const sim::IDE& ide: ides){
    //     const art::Ptr<simb::MCTruth>& mct = bt->TrackIDToMCTruth(eve.CalculateEveId(ide.trackID));
    //     if(mct->GetNeutrino().Nu().PdgCode() == 12){
    //       hitChans.insert(hit.Channel());
    //       break;
    //     }
    //   }
    // }

    // gTrueDepositE = 0;
    // gTrueVisE = 0;
    // for(const sim::SimChannel& simchan: *simchans){
    //   // Only count truth on the collection plane
    //   if(geom.View(simchan.Channel()) != geo::kZ) continue;

    //   for(const auto& it: simchan.TDCIDEMap()){
    //     for(const sim::IDE& ide: it.second){
    //       const art::Ptr<simb::MCTruth>& mct = bt->TrackIDToMCTruth(eve.CalculateEveId(ide.trackID));

    //       if(mct->GetNeutrino().Nu().PdgCode() == 12){
    //         gTrueDepositE += ide.energy;
    //         if(hitChans.count(simchan.Channel())) gTrueVisE += ide.energy;
    //       }
    //     }
    //   }
    // }



    std::vector<Pt2D> pts;

    // cm / tick
    const double driftVel = dp.DriftVelocity() * dp.SamplingRate() / 1000;

    //    const double driftLen = geom.Cryostat(0).TPC(0).DriftDistance();
    //    const double driftT = driftLen / dp.DriftVelocity();

    double totTrueQ = 0;
    double totHitQ = 0;

    for(unsigned int hitIdx = 0; hitIdx < hits->size(); ++hitIdx){
      const recob::Hit& hit = (*hits)[hitIdx];

      if(hit.View() != geo::kZ) continue; // vertical ie collection

      //      std::vector<sim::TrackIDE> eveIDs = bt->HitToEveID(art::Ptr<recob::Hit>(hits, hitIdx));

      bool isSig = false;
      // std::vector<sim::IDE> ides;
      // bt->HitToSimIDEs(hit, ides);
      // for(const sim::IDE& ide: ides){
      //   const art::Ptr<simb::MCTruth>& mct = bt->TrackIDToMCTruth(eve.CalculateEveId(ide.trackID));

      //   if(mct->GetNeutrino().Nu().PdgCode() == 12){
      //     isSig = true;
      //   }
      // }

      const geo::WireGeo* wiregeo = geom.WirePtr(hit.WireID());
      double xyz[3];
      wiregeo->GetCenter(xyz);

      // Apparently this is in ticks. We want it in time units.
      const double t = hit.PeakTime() * dp.SamplingRate() / 1000;

      const double q = hit.Integral();
      totHitQ += q;

      pts.push_back(Pt2D(xyz[2], t*driftVel, q, isSig, hitIdx));
      if(isSig) totTrueQ += q;
    }

    std::vector<std::vector<Pt2D>> slices = DBSCAN(pts);

    for(const std::vector<Pt2D>& slice: slices){
      if(slice.empty()) continue; // TODO - how does this happen?

      sn::SNSlice prod;

      prod.meanT = 0;
      prod.meanZ = 0;
      prod.totQ = 0;

      for(const Pt2D& p: slice){
        prod.meanT += p.q * p.y / driftVel;
        prod.meanZ += p.q * p.x;
        prod.totQ += p.q;

        prod.hits.push_back(art::Ptr<recob::Hit>(hits, p.idx));
      }

      prod.meanT /= prod.totQ;
      prod.meanZ /= prod.totQ;
      
      slicecol->push_back(prod);
    }
    

    // gEfficiency = 0;
    // gPurity = 0;

    // for(unsigned int i = 0; i < slices.size(); ++i){
    //   double sliceQ = 0;
    //   double pur = 0;
    //   double eff = 0;
    //   for(unsigned int j = 0; j < slices[i].size(); ++j){
    //     sliceQ += slices[i][j].q;
    //     if(slices[i][j].isSig){
    //       pur += slices[i][j].q;
    //       eff += slices[i][j].q;
    //     }
    //   }
    // }

    evt.put(std::move(slicecol));
  }

  } // namespace
