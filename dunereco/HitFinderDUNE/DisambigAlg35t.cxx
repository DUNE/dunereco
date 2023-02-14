////////////////////////////////////////////////////////////////////////
//
// DisambigAlg35t.cxx
//
// trj@fnal.gov
// tjyang@fnal.gov
//
// description
//
// Based on Tom Junk's idea of triplets matching
//
//
//
////////////////////////////////////////////////////////////////////////


//Framework includes:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "DisambigAlg35t.h"
#include "larevt/Filters/ChannelFilter.h"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "TH1D.h"

namespace dune{

  DisambigAlg35t::DisambigAlg35t(fhicl::ParameterSet const& pset)
    : fDBScan(pset.get< fhicl::ParameterSet >("DBScanAlg"))
  {
    fTimeCut = pset.get<double>("TimeCut");
    fDistanceCut = pset.get<double>("DistanceCut");
    fDistanceCutClu = pset.get<double>("DistanceCutClu");
    fTimeWiggle     = pset.get<double>("TimeWiggle");
    fColChanWiggle  = pset.get<int>("ColChannelWiggle");
    fDoCleanUpHits  = pset.get<bool>("DoCleanUpHits", true);
  }


  //----------------------------------------------------------
  //----------------------------------------------------------
  void DisambigAlg35t::RunDisambig(detinfo::DetectorClocksData const& clockData,
                                   detinfo::DetectorPropertiesData const& detProp,
                                   const std::vector< art::Ptr<recob::Hit> > &OrigHits   )
  {
    fDisambigHits.clear();

    art::ServiceHandle<geo::Geometry> geo;

    std::vector<std::vector<art::Ptr<recob::Hit> > > hitsUV(2);
    std::vector<art::Ptr<recob::Hit> > hitsZ;

    //  TH1D *histu = new TH1D("histu","histu",4000,0,4000);
    //  TH1D *histv = new TH1D("histv","histv",4000,0,4000);
    //  TH1D *histz = new TH1D("histz","histz",4000,0,4000);

    for (size_t i = 0; i<OrigHits.size(); ++i){
      switch (OrigHits[i]->View()){
      case geo::kU:
        hitsUV[0].push_back(OrigHits[i]);
        //      if (OrigHits[i]->WireID().TPC==1) 
        //	histu->Fill(OrigHits[i]->PeakTime()
        //		    - detProp.GetXTicksOffset(0,
        //					       OrigHits[i]->WireID().TPC,
        //					       OrigHits[i]->WireID().Cryostat)
        //		    ,OrigHits[i]->Charge());
        break;
      case geo::kV:
        hitsUV[1].push_back(OrigHits[i]);
        //      if (OrigHits[i]->WireID().TPC==1) 
        //	histv->Fill(OrigHits[i]->PeakTime()
        //		    - detProp.GetXTicksOffset(1,
        //					       OrigHits[i]->WireID().TPC,
        //					       OrigHits[i]->WireID().Cryostat)
        //		    ,OrigHits[i]->Charge());
        break;
      case geo::kZ:
        hitsZ.push_back(OrigHits[i]);
        //      if (OrigHits[i]->WireID().TPC==1) 
        //	histz->Fill(OrigHits[i]->PeakTime()
        //		    - detProp.GetXTicksOffset(2,
        //					       OrigHits[i]->WireID().TPC,
        //					       OrigHits[i]->WireID().Cryostat)
        //		    ,OrigHits[i]->Charge());
        break;
      default:
        throw cet::exception("DisambigAlg35t")
          <<": hit view unkonw. \n";
      }
    }
    //  std::cout<<histu->GetMean()<<" "<<histv->GetMean()<<" "<<histz->GetMean()<<std::endl;
    //  delete histu;
    //  delete histv;
    //  delete histz;
    std::vector<std::map<unsigned int, unsigned int> >fHasBeenDisambigedUV(2);
    const unsigned int ntpc = geo->NTPC();
    //hits and wireids for DBScan
    std::vector<std::vector<art::Ptr<recob::Hit> > > allhitsu(ntpc);
    std::vector<std::vector<art::Ptr<recob::Hit> > > allhitsv(ntpc);
    std::vector<std::vector<art::Ptr<recob::Hit> > > allhitsz(ntpc);
    std::vector<std::vector<geo::WireID> > wireidsu(ntpc);
    std::vector<std::vector<geo::WireID> > wireidsv(ntpc);
  
    std::vector<std::vector<unsigned int> > cluidu(ntpc);
    std::vector<std::vector<unsigned int> > cluidv(ntpc);
    std::vector<std::vector<unsigned int> > bestwireidu(ntpc);
    std::vector<std::vector<unsigned int> > bestwireidv(ntpc);

    //Look for triplets of U,V,Z hits that are common in time
    //Try all possible wire segments for U and V hits and 
    //see if U, V, Z wires cross
    for (size_t z = 0; z<hitsZ.size(); ++z){//loop over z hits
      //find triplets of U,V,Z hits that are common in time
      unsigned int apaz(0), cryoz(0);
      fAPAGeo.ChannelToAPA(hitsZ[z]->Channel(), apaz, cryoz);

      double tz = hitsZ[z]->PeakTime()
        - detProp.GetXTicksOffset(hitsZ[z]->WireID().Plane,
                                   hitsZ[z]->WireID().TPC,
                                   hitsZ[z]->WireID().Cryostat);
      //loop over u hits
      bool findmatch = false;
      for (size_t u = 0; u<hitsUV[0].size() && !findmatch; ++u){
        if (fHasBeenDisambigedUV[0].find(u)!=fHasBeenDisambigedUV[0].end()) continue;
        unsigned int apau(0), cryou(0);
        fAPAGeo.ChannelToAPA(hitsUV[0][u]->Channel(), apau, cryou);
        if (apau!=apaz) continue;
        double tu = hitsUV[0][u]->PeakTime()
          - detProp.GetXTicksOffset(0,
                                     hitsZ[z]->WireID().TPC,
                                     hitsZ[z]->WireID().Cryostat);
       
        if (std::abs(tu-tz)<fTimeCut){
          //find a matched u hit, loop over v hits
          for (size_t v = 0; v<hitsUV[1].size(); ++v){
            if (fHasBeenDisambigedUV[1].find(v)!=fHasBeenDisambigedUV[1].end()) continue;
            unsigned int apav(0), cryov(0);
            fAPAGeo.ChannelToAPA(hitsUV[1][v]->Channel(), apav, cryov);
            if (apav!=apaz) continue;
            double tv = hitsUV[1][v]->PeakTime()
              - detProp.GetXTicksOffset(1,
                                         hitsZ[z]->WireID().TPC,
                                         hitsZ[z]->WireID().Cryostat);
    
            if (std::abs(tv-tz)<fTimeCut){
	    
              //find a matched v hit, see if the 3 wire segments cross
              geo::WireID zwire = geo->ChannelToWire(hitsZ[z]->Channel())[0];
              std::vector<geo::WireID>  uwires = geo->ChannelToWire(hitsUV[0][u]->Channel());
              std::vector<geo::WireID>  vwires = geo->ChannelToWire(hitsUV[1][v]->Channel());
	    
              unsigned int totalintersections = 0;
              unsigned int bestu = 0;  //index to wires associated with channel
              unsigned int bestv = 0;

              for (size_t uw = 0; uw<uwires.size(); ++uw){
                for (size_t vw = 0; vw<vwires.size(); ++vw){
                  geo::WireIDIntersection widiuv;
                  geo::WireIDIntersection widiuz;
                  geo::WireIDIntersection widivz;

                  if (uwires[uw].TPC!=vwires[vw].TPC) continue;
                  if (uwires[uw].TPC!=zwire.TPC) continue;
                  if (vwires[vw].TPC!=zwire.TPC) continue;

                  if (!geo->WireIDsIntersect(zwire,uwires[uw],widiuz)) continue;
                  if (!geo->WireIDsIntersect(zwire,vwires[vw],widivz)) continue;
                  if (!geo->WireIDsIntersect(uwires[uw],vwires[vw],widiuv)) continue;

                  double dis1 = sqrt(pow(widiuz.y-widivz.y,2)+pow(widiuz.z-widivz.z,2));
                  double dis2 = sqrt(pow(widiuz.y-widiuv.y,2)+pow(widiuz.z-widiuv.z,2));
                  double dis3 = sqrt(pow(widiuv.y-widivz.y,2)+pow(widiuv.z-widivz.z,2));
                  double maxdis = std::max(dis1,dis2);
                  maxdis = std::max(maxdis,dis3);

                  if (maxdis<fDistanceCut){
                    ++totalintersections;
                    bestu = uw;
                    bestv = vw;
                  }
                }
              }
              if (totalintersections==1){
                allhitsu[uwires[bestu].TPC].push_back(hitsUV[0][u]);
                allhitsv[vwires[bestv].TPC].push_back(hitsUV[1][v]);
                allhitsz[hitsZ[z]->WireID().TPC].push_back(hitsZ[z]);
                wireidsu[uwires[bestu].TPC].push_back(uwires[bestu]);
                wireidsv[vwires[bestv].TPC].push_back(vwires[bestv]);
                cluidu[uwires[bestu].TPC].push_back(u);
                cluidv[vwires[bestv].TPC].push_back(v);
                bestwireidu[uwires[bestu].TPC].push_back(1+bestu);
                bestwireidv[vwires[bestv].TPC].push_back(1+bestv);	      
                findmatch = true;
                break;
              }
            }//find v hit consistent in time
          }//loop over all v hits
        }//find u hit consistent in time
      }//loop over all u hits
    }//loop over all z hits
    //Done finding trivial disambiguated hits

    if (fDoCleanUpHits){
      //running DB scan to identify and remove outlier hits
      // get the ChannelFilter
      filter::ChannelFilter chanFilt;
      double CorrectedHitTime = 0;
      int ChannelNumber = 0;
      for (size_t i = 0; i<ntpc; ++i){//loop over all TPCs
        if (!allhitsu[i].size()) continue;

        //initialize dbscan with all hits in u view
        fDBScan.InitScan(clockData, detProp, allhitsu[i], chanFilt.SetOfBadChannels(), wireidsu[i]);
        //run dbscan
        fDBScan.run_cluster();
        //fpointId_to_clusterId maps a hit index to a cluster index (which can be negative if the hit is not associated with any clusters)
        if (allhitsu[i].size()!=fDBScan.fpointId_to_clusterId.size())
          throw cet::exception("DisambigAlg35t") <<"DBScan hits do not match input hits"<<allhitsu[i].size()<<" "<<fDBScan.fpointId_to_clusterId.size()<<"\n";

        //Find out which clusters are unique in time
        std::vector<unsigned int> hitstoaddu;
        std::vector<unsigned int> dbcluhits(fDBScan.fclusters.size(),0);
        std::vector< bool > boolVector(fDBScan.fclusters.size(), true);
        std::map<int, std::pair <double,double> > ClusterStartEndTime;
        std::map<int, std::pair <int, int> > ClusterStartEndColChan;
        for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
	  // for c2: fDBScan.fpointId_to_clusterId[j]>=0 is always true
          //if (fDBScan.fpointId_to_clusterId[j]>=0&&fDBScan.fpointId_to_clusterId[j]<fDBScan.fclusters.size()) {
          if (fDBScan.fpointId_to_clusterId[j]<fDBScan.fclusters.size()) {
            ++dbcluhits[fDBScan.fpointId_to_clusterId[j]];

            CorrectedHitTime = allhitsu[i][j]->PeakTime() - detProp.GetXTicksOffset(allhitsu[i][j]->WireID().Plane, allhitsu[i][j]->WireID().TPC, allhitsu[i][j]->WireID().Cryostat);
            if (ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first==0 || ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first > CorrectedHitTime) 
              ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first = CorrectedHitTime;
            if (ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second==0 || ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second < CorrectedHitTime) 
              ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second = CorrectedHitTime;

            ChannelNumber = allhitsz[i][j]->Channel();
            if (ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first==0 || ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first > ChannelNumber) 
              ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first = ChannelNumber;
            if (ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second==0 || ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second < ChannelNumber) 
              ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second = ChannelNumber;
            /*
              std::cout << "Looking at Cluster " << fDBScan.fpointId_to_clusterId[j] << " in TPC " << (int)i
              << ", It has start time " << ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first 
              << ", and end time " << ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second 
              << ", It has start ColChan " << ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first 
              << ", and end ColChan " << ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second 
              << std::endl;
            //*/
          }
        }
        for ( unsigned int ClusNum = 0; ClusNum < fDBScan.fclusters.size(); ++ClusNum ) {
          for ( unsigned int ClusIt = 0; ClusIt < fDBScan.fclusters.size(); ++ClusIt ) {
            if ( dbcluhits[ClusIt] > dbcluhits[ClusNum] ) { 
              // Smaller cluster, so check not in the same time range...
              if ( ClusterStartEndTime[ClusNum].first > ( ClusterStartEndTime[ClusIt].first + fTimeWiggle ) 
                   && ClusterStartEndTime[ClusNum].first < ( ClusterStartEndTime[ClusIt].second + fTimeWiggle ) ) {
                if ( ClusterStartEndTime[ClusNum].second > ( ClusterStartEndTime[ClusIt].first + fTimeWiggle ) 
                     && ClusterStartEndTime[ClusNum].second < ( ClusterStartEndTime[ClusIt].second + fTimeWiggle ) ) {
                  // Within the same time...now check if same collection channel range...
                  if ( ClusterStartEndColChan[ClusNum].first > ( ClusterStartEndColChan[ClusIt].first + fColChanWiggle ) 
                       && ClusterStartEndColChan[ClusNum].first < ( ClusterStartEndColChan[ClusIt].second + fColChanWiggle ) ) {
                    if ( ClusterStartEndColChan[ClusNum].second > ( ClusterStartEndColChan[ClusIt].first + fColChanWiggle ) 
                         && ClusterStartEndColChan[ClusNum].second < ( ClusterStartEndColChan[ClusIt].second + fColChanWiggle ) ) {
                      // Within same time and collection channel range...Bad Cluster?
                      boolVector[ClusNum] = false;
                    }
                  }
                }
              }
            }
          }
          /*
            std::cout << "\nLooking at Cluster Number " << ClusNum << ".\n" 
            << "It had start time " << ClusterStartEndTime[ClusNum].first << ", and end time " << ClusterStartEndTime[ClusNum].second << ".\n"
            << "It had start ColChan " << ClusterStartEndColChan[ClusNum].first << ", and end ColChan " << ClusterStartEndColChan[ClusNum].second << ".\n"
            << "The bool vector value for this cluster is..." << boolVector[ClusNum]
            << std::endl;
          //*/
        }

        //std::cout << "Now going to look at v hits....." << std::endl;

        //hitstoaddu holds all u hits in the time separated clusters 
        for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
          if (int(fDBScan.fpointId_to_clusterId[j]) == -1 ) continue;
          if (boolVector[fDBScan.fpointId_to_clusterId[j]] == true ) hitstoaddu.push_back(j);
        }

        //Now to do the same for v hits
        fDBScan.InitScan(clockData, detProp, allhitsv[i], chanFilt.SetOfBadChannels(), wireidsv[i]);
        fDBScan.run_cluster();
        if (allhitsv[i].size()!=fDBScan.fpointId_to_clusterId.size())
          throw cet::exception("DisambigAlg35t") <<"DBScan hits do not match input hits"<<allhitsv[i].size()<<" "<<fDBScan.fpointId_to_clusterId.size()<<"\n";

        std::vector<unsigned int> hitstoaddv;
        dbcluhits.clear();
        boolVector.clear();
        ClusterStartEndTime.clear();
        ClusterStartEndColChan.clear();
        dbcluhits.resize(fDBScan.fclusters.size(),0);
        boolVector.resize(fDBScan.fclusters.size(),true);
        for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
	  //for c2: fDBScan.fpointId_to_clusterId[j]>=0 is always true
          //if (fDBScan.fpointId_to_clusterId[j]>=0&&fDBScan.fpointId_to_clusterId[j]<fDBScan.fclusters.size()) {
          if (fDBScan.fpointId_to_clusterId[j]<fDBScan.fclusters.size()) {
            ++dbcluhits[fDBScan.fpointId_to_clusterId[j]];

            CorrectedHitTime = allhitsv[i][j]->PeakTime() - detProp.GetXTicksOffset(allhitsv[i][j]->WireID().Plane, allhitsv[i][j]->WireID().TPC, allhitsv[i][j]->WireID().Cryostat);
            if (ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first==0 || ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first > CorrectedHitTime) 
              ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first = CorrectedHitTime;
            if (ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second==0 || ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second < CorrectedHitTime) 
              ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second = CorrectedHitTime;

            ChannelNumber = allhitsz[i][j]->Channel();
            if (ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first==0 || ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first > ChannelNumber) 
              ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first = ChannelNumber;
            if (ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second==0 || ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second < ChannelNumber) 
              ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second = ChannelNumber;
            /*
              std::cout << "Looking at Cluster " << fDBScan.fpointId_to_clusterId[j] << " in TPC " << (int)i
              << ", It has start time " << ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].first 
              << ", and end time " << ClusterStartEndTime[fDBScan.fpointId_to_clusterId[j]].second 
              << ", It has start ColChan " << ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].first 
              << ", and end ColChan " << ClusterStartEndColChan[fDBScan.fpointId_to_clusterId[j]].second 
              << std::endl;
	
            //*/
          }
        }
        for ( unsigned int ClusNum = 0; ClusNum < fDBScan.fclusters.size(); ++ClusNum ) {
          for ( unsigned int ClusIt = 0; ClusIt < fDBScan.fclusters.size(); ++ClusIt ) {
            if ( dbcluhits[ClusIt] > dbcluhits[ClusNum] ) { 
              // Smaller cluster, so check not in the same time range...
              if ( ClusterStartEndTime[ClusNum].first > ( ClusterStartEndTime[ClusIt].first + fTimeWiggle ) 
                   && ClusterStartEndTime[ClusNum].first < ( ClusterStartEndTime[ClusIt].second + fTimeWiggle ) ) {
                if ( ClusterStartEndTime[ClusNum].second > ( ClusterStartEndTime[ClusIt].first + fTimeWiggle ) 
                     && ClusterStartEndTime[ClusNum].second < ( ClusterStartEndTime[ClusIt].second + fTimeWiggle ) ) {
                  // Within the same time...now check if same collection channel range...
                  if ( ClusterStartEndColChan[ClusNum].first > ( ClusterStartEndColChan[ClusIt].first + fColChanWiggle ) 
                       && ClusterStartEndColChan[ClusNum].first < ( ClusterStartEndColChan[ClusIt].second + fColChanWiggle ) ) {
                    if ( ClusterStartEndColChan[ClusNum].second > ( ClusterStartEndColChan[ClusIt].first + fColChanWiggle ) 
                         && ClusterStartEndColChan[ClusNum].second < ( ClusterStartEndColChan[ClusIt].second + fColChanWiggle ) ) {
                      // Within same time and collection channel range...Bad Cluster?	      
                      boolVector[ClusNum] = false;
                    }
                  }
                }
              }
            }
          }
          /*
            std::cout << "\nLooking at Cluster Number " << ClusNum << ".\n" 
            << "It had start time " << ClusterStartEndTime[ClusNum].first << ", and end time " << ClusterStartEndTime[ClusNum].second << ".\n"
            << "It had start ColChan " << ClusterStartEndColChan[ClusNum].first << ", and end ColChan " << ClusterStartEndColChan[ClusNum].second << ".\n"
            << "The bool vector value for this cluster is..." << boolVector[ClusNum] << std::endl;
          //*/
        }
        //hitstoaddv holds all v hits in the time separated clusters
        for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j){
          if (int(fDBScan.fpointId_to_clusterId[j]) == -1 ) continue;
          if (boolVector[fDBScan.fpointId_to_clusterId[j]] == true ) hitstoaddv.push_back(j);
        }
        //std::cout << "Looking back through all hits to see if can add hits? Loops through " << hitstoaddu.size() << " hits for U, and " << hitstoaddv.size() << " for V." << std::endl;
        //size(allhitsu)=size(allhitsv)
        for (size_t j = 0; j < allhitsu[i].size(); ++j){
          bool foundu = false;
          bool foundv = false;
          for (size_t k = 0; k<hitstoaddu.size(); ++k){
            if (j==hitstoaddu[k]) {
              foundu = true;
            }
          }
          for (size_t k = 0; k<hitstoaddv.size(); ++k){
            if (j==hitstoaddv[k]) {
              foundv = true;
            }
          }
          if (foundu&&foundv){
            //only add each hit once
            if (fHasBeenDisambigedUV[0].find(cluidu[i][j])==fHasBeenDisambigedUV[0].end()){
              fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(allhitsu[i][j],wireidsu[i][j]));
              fHasBeenDisambigedUV[0][cluidu[i][j]] = bestwireidu[i][j];
            }
            if (fHasBeenDisambigedUV[1].find(cluidv[i][j])==fHasBeenDisambigedUV[1].end()){

              fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(allhitsv[i][j],wireidsv[i][j]));
              fHasBeenDisambigedUV[1][cluidv[i][j]] = bestwireidv[i][j];
            }
          }
        }
      }

    }//if (fDoCleanUpHits)
    else{//just take all triplets of hits
      for (size_t i = 0; i<ntpc; ++i){//loop over all TPCs
        for (size_t j = 0; j < allhitsu[i].size(); ++j){
          //only add each hit once
          if (fHasBeenDisambigedUV[0].find(cluidu[i][j])==fHasBeenDisambigedUV[0].end()){
            fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(allhitsu[i][j],wireidsu[i][j]));
            fHasBeenDisambigedUV[0][cluidu[i][j]] = bestwireidu[i][j];
          }
          if (fHasBeenDisambigedUV[1].find(cluidv[i][j])==fHasBeenDisambigedUV[1].end()){
          
            fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(allhitsv[i][j],wireidsv[i][j]));
            fHasBeenDisambigedUV[1][cluidv[i][j]] = bestwireidv[i][j];
          }
        }
      }
    }
    
    //loop over undisambiguated hits, find the nearest channel of disambiguated hits and determine the correct wire segment.
    for (size_t i = 0; i<2; ++i){//loop over U and V hits
      for (size_t hit = 0; hit<hitsUV[i].size(); ++hit){
        if (fHasBeenDisambigedUV[i].find(hit)!=fHasBeenDisambigedUV[i].end()) continue;
        unsigned int apa1(0), cryo1(0);
        fAPAGeo.ChannelToAPA(hitsUV[i][hit]->Channel(), apa1, cryo1);
        //unsigned int channdiff = 100000;
        //geo::WireID nearestwire;
        std::vector<geo::WireID> wires = geo->ChannelToWire(hitsUV[i][hit]->Channel());
        //std::vector<double> disttoallhits(wires.size(),-1);
        std::vector<int> nearbyhits(wires.size(),-1);
        for (auto& u2 : fHasBeenDisambigedUV[i]){
          unsigned int apa2(0), cryo2(0);
          fAPAGeo.ChannelToAPA(hitsUV[i][u2.first]->Channel(), apa2, cryo2);
          if (apa1!=apa2) continue;
          geo::WireID hitwire = geo->ChannelToWire(hitsUV[i][u2.first]->Channel())[u2.second-1];
          double wire_pitch = geo->WirePitch(hitwire.asPlaneID());    //wire pitch in cm
          for (size_t w = 0; w<wires.size(); ++w){
            if (wires[w].TPC!= hitwire.TPC) continue;
            if (nearbyhits[w]<0) nearbyhits[w] = 0;
            double distance = sqrt(pow((double(hitwire.Wire)-double(wires[w].Wire))*wire_pitch,2)+pow(detProp.ConvertTicksToX(hitsUV[i][u2.first]->PeakTime(),hitwire.Plane,hitwire.TPC,hitwire.Cryostat)-detProp.ConvertTicksToX(hitsUV[i][hit]->PeakTime(),hitwire.Plane,hitwire.TPC,hitwire.Cryostat),2));
            if (distance<fDistanceCutClu) ++nearbyhits[w];
          }
        }
        unsigned bestwire = 0;
        int maxnumhits = 0;
        for (size_t w = 0; w<wires.size(); ++w){
          if (nearbyhits[w]<0) continue;
          if (nearbyhits[w]>maxnumhits){
            bestwire = w;
            maxnumhits = nearbyhits[w];
          }
        }
        if (maxnumhits){
          fDisambigHits.push_back(std::pair<art::Ptr<recob::Hit>, geo::WireID>(hitsUV[i][hit],wires[bestwire]));
          fHasBeenDisambigedUV[i][hit] = 1+bestwire;
        }
        else{
          mf::LogInfo("DisambigAlg35t")<<"Could not find disambiguated hit for  "<<std::string(hitsUV[i][hit]->WireID())<<" "<<hitsUV[i][hit]->PeakTime();
        }
      }
    }  
  }

} //end namespace apa
