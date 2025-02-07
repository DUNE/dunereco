////////////////////////////////////////////////////////////////////////
/// \file    TransformerPixelMapProducer.h
/// \brief   TransformerPixelMapProducer for TransformerCVN modified from RegPixelMapProducer.h
/// \author  Alejandro Yankelevich - ayankele@uci.edu
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>

#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunereco/TransformerCVN/art/TransformerPixelMapProducer.h"
#include  "TVector2.h"

namespace cnn
{

  TransformerPixelMapProducer::TransformerPixelMapProducer(unsigned int nWire, unsigned int wRes, unsigned int nTdc, double tRes, int Global):
  fNWire(nWire),
  fWRes(wRes),
  fNTdc(nTdc),
  fTRes(tRes),
  fGlobalWireMethod(Global),
  fOffset{0,0}
  {}

  TransformerPixelMap TransformerPixelMapProducer::CreateMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster,
                                             art::FindManyP<recob::Wire> const& fmwire,
                                             art::FindManyP<recob::Track> const& fmtrkhit,
                                             const std::vector<float> &vtx,
                                             int ProngID)
  {
      // tag prong by track
      // Create pixel maps around the vertex
      hitwireidx.clear();
      tmin_each_wire.clear();
      tmax_each_wire.clear();
      trms_max_each_wire.clear();
      fOffset[0] = 0; fOffset[1] = 0;

      RegCNNBoundary bound = DefineBoundary(detProp, cluster, vtx);

      return CreateMapGivenBoundaryByHit(clockData, detProp, cluster, bound, fmwire, fmtrkhit, ProngID);
  }

  TransformerPixelMap TransformerPixelMapProducer::CreateMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster,
                                             art::FindManyP<recob::Wire> const& fmwire,
                                             art::FindManyP<recob::Shower> const& fmshwhit,
                                             const std::vector<float> &vtx,
                                             int ProngID)
  {
      // tag prong by shower
      // Create pixel maps around the vertex
      hitwireidx.clear();
      tmin_each_wire.clear();
      tmax_each_wire.clear();
      trms_max_each_wire.clear();
      fOffset[0] = 0; fOffset[1] = 0;

      RegCNNBoundary bound = DefineBoundary(detProp, cluster, vtx);

    return CreateMapGivenBoundaryByHit(clockData, detProp, cluster, bound, fmwire, fmshwhit, ProngID);
  }

  TransformerPixelMap TransformerPixelMapProducer::CreateMapGivenBoundaryByHit(detinfo::DetectorClocksData const& clockData,
                                                          detinfo::DetectorPropertiesData const& detProp,
                                                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                          const RegCNNBoundary& bound,
                                                          art::FindManyP<recob::Wire> const& fmwire,
                                                          art::FindManyP<recob::Track> const& fmtrkhit,
                                                          int ProngID)
  {

      TransformerPixelMap pm(fNWire, fWRes, fNTdc, fTRes, bound, ProngID);

      if (!fmwire.isValid()) return pm;

      // Loop over hits
      unsigned int nhits = cluster.size();
      for (unsigned int ihit= 0; ihit< nhits; ++ihit) {
          art::Ptr<recob::Hit> hit = cluster.at(ihit);
          geo::WireID wireid = hit->WireID();

          unsigned int planeid = wireid.Plane;

          double peaktime = hit->PeakTime() ;
          if (wireid.TPC%2 == 0) peaktime = -peaktime;

          unsigned int globalWire = 0;
          unsigned int globalplane = planeid;
          if (fGlobalWireMethod == 1) {
              globalWire = GetGlobalWire(wireid);
              if (wireid.TPC%2 == 1) {
                  if (wireid.Plane == 0) globalWire += fOffset[0];
                  if (wireid.Plane == 1) globalWire += fOffset[1];
              }
          } else if (fGlobalWireMethod == 2) {
              globalWire = wireid.Wire;
              GetDUNEGlobalWireTDC(detProp, wireid, cluster[ihit]->PeakTime(), globalWire, globalplane, peaktime);
          } else {
              std::cout<<"Wrong GlobalWireMethod"<<std::endl;
              abort();
          }

          double correctedHitCharge = (hit->Integral()*TMath::Exp((sampling_rate(clockData)*hit->PeakTime()) / (detProp.ElectronLifetime()*1.e3) ) );

          int hit_prong_tag = -1;
          if (fmtrkhit.isValid() && fmtrkhit.at(ihit).size()!=0) {
              hit_prong_tag = fmtrkhit.at(ihit)[0]->ID();
          }

          pm.Add((int)globalWire, (int)round(peaktime), globalplane, correctedHitCharge, wireid.TPC, hit_prong_tag);
      }

      pm.Finish();

      return pm;

  }

  TransformerPixelMap TransformerPixelMapProducer::CreateMapGivenBoundaryByHit(detinfo::DetectorClocksData const& clockData,
                                                          detinfo::DetectorPropertiesData const& detProp,
                                                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                          const RegCNNBoundary& bound,
                                                          art::FindManyP<recob::Wire> const& fmwire,
                                                          art::FindManyP<recob::Shower> const& fmshwhit,
                                                          int ProngID)
  {

      TransformerPixelMap pm(fNWire, fWRes, fNTdc, fTRes, bound, ProngID);

      if (!fmwire.isValid()) return pm;

      // Loop over hits
      unsigned int nhits = cluster.size();
      for (unsigned int ihit= 0; ihit< nhits; ++ihit) {
          art::Ptr<recob::Hit> hit = cluster.at(ihit);
          geo::WireID wireid = hit->WireID();

          unsigned int planeid = wireid.Plane;

          double peaktime = hit->PeakTime() ;
          if (wireid.TPC%2 == 0) peaktime = -peaktime;

          unsigned int globalWire = 0;
          unsigned int globalplane = planeid;
          if (fGlobalWireMethod == 1) {
              globalWire = GetGlobalWire(wireid);
              if (wireid.TPC%2 == 1) {
                  if (wireid.Plane == 0) globalWire += fOffset[0];
                  if (wireid.Plane == 1) globalWire += fOffset[1];
              }
          } else if (fGlobalWireMethod == 2) {
              globalWire = wireid.Wire;
              GetDUNEGlobalWireTDC(detProp, wireid, cluster[ihit]->PeakTime(), globalWire, globalplane, peaktime);
          } else {
              std::cout<<"Wrong GlobalWireMethod"<<std::endl;
              abort();
          }

          double correctedHitCharge = (hit->Integral()*TMath::Exp((sampling_rate(clockData)*hit->PeakTime()) / (detProp.ElectronLifetime()*1.e3) ) );

          int hit_prong_tag = -1;
          if (fmshwhit.isValid() && fmshwhit.at(ihit).size()!=0) {
              hit_prong_tag = fmshwhit.at(ihit)[0]->ID();
          }

          pm.Add((int)globalWire, (int)round(peaktime), globalplane, correctedHitCharge, wireid.TPC, hit_prong_tag);
      }

      std::cout<<"FIXME: Produce prong pixelmap from fmshwhit" << std::endl;
      pm.Finish();

      return pm;

  }


  std::ostream& operator<<(std::ostream& os, const TransformerPixelMapProducer& p)
  {
    os << "TransformerPixelMapProducer: "
       << p.NTdc()  <<" tdcs X  " <<  p.NWire() << " wires";
    return os;
  }

  RegCNNBoundary TransformerPixelMapProducer::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                                     std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                     const std::vector<float> &vtx)
  {

      unsigned int temp_wire = cluster[0]->WireID().Wire;
      float temp_time_min = 1e5; //99999;
      float temp_time_max = 0;   //-99999;
      float temp_trms_max = 0;   //-99999;
      for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
      {
          geo::WireID wireid = cluster[iHit]->WireID();
          if (temp_wire != wireid.Wire ){ // lost last hit info
              temp_wire = wireid.Wire;
              hitwireidx.push_back(iHit-1);
              tmin_each_wire.push_back(int(temp_time_min));
              tmax_each_wire.push_back(int(temp_time_max));
              trms_max_each_wire.push_back(temp_trms_max);
             // temp_time_min = 99999; temp_time_max = -99999, temp_trms_max = -99999;
		      temp_time_min = 1e5; temp_time_max = 0; temp_trms_max = 0;
          }
          if (temp_time_min > cluster[iHit]->PeakTime()) temp_time_min = cluster[iHit]->PeakTime();
          if (temp_time_max < cluster[iHit]->PeakTime()) temp_time_max = cluster[iHit]->PeakTime();
          if (temp_trms_max < cluster[iHit]->RMS()) temp_trms_max = (float)cluster[iHit]->RMS();
      }

      double center_wire[3] = {-99999, -99999, -99999};
      double center_tick[3] = {-99999, -99999, -99999};
      unsigned int size_vtx = (unsigned int) vtx.size();
      if (size_vtx == 3){
          geo::Point_t const regvtx_loc{(double)vtx[0], (double)vtx[1], (double)vtx[2]};
          int rawcrys = 0;
          bool inTPC = false;
          if (geom->FindTPCAtPosition(regvtx_loc).isValid) inTPC = true;
          if (inTPC){
              for (int iplane = 0; iplane<3; iplane++){
                  int rawtpc = (int) (geom->FindTPCAtPosition(regvtx_loc)).TPC;
                  geo::PlaneGeo const& planegeo_temp = wireReadout->Plane(geo::PlaneID(0, 0, iplane));
                  geo::PlaneID const planeID(rawcrys, rawtpc, iplane);
                  geo::WireID w1;
                  try { 
                      w1 = wireReadout->Plane(planeID).NearestWireID(regvtx_loc);
                  }
                  catch (geo::InvalidWireError const& e){
                      if (!e.hasSuggestedWire()) throw;
                      w1 = planegeo_temp.ClosestWireID(e.suggestedWireID());
                  }
                  double time1 = detProp.ConvertXToTicks(regvtx_loc.X(), planeID);
                  if (fGlobalWireMethod == 1){
                      std::cout << "You Can't use this with GlobalWireMethod = 1" << std::endl;
                      abort();
                  } else if (fGlobalWireMethod == 2){
                      unsigned int globalWire  = (double)w1.Wire;
                      unsigned int globalPlane = w1.Plane;
                      double globalTDC = time1;
                      GetDUNEGlobalWireTDC(detProp, w1, time1, globalWire, globalPlane, globalTDC);
                      center_wire[globalPlane] = globalWire;
                      //center_tick[globalPlane] = (double)globalTDC;
                      center_tick[globalPlane] = (int)globalTDC;
                  } else {
                      std::cout << "Wrong Global Wire Method" << std::endl;
                      abort();
                  }
              } // end of iplane
          } // end of inTPC
      } else if ( size_vtx == 6){
	    for (int ii = 0; ii < 3; ii++){
	        center_wire[ii] = vtx[ii*2+1];
	        center_tick[ii] = vtx[ii*2];
	    }
      } else{
	      std::cout << "Wrong reconstructed vertex" << std::endl;
	      std::cout << "-->" << size_vtx << std::endl;
      }


      int shift = 56*fWRes;
      center_wire[0] += fNWire*fWRes/2 - shift;
      center_wire[1] -= fNWire*fWRes/2 - shift;
      center_wire[2] += fNWire*fWRes/2 - shift;

      RegCNNBoundary bound(fNWire,fNTdc,fWRes,fTRes,round(center_wire[0]),round(center_wire[1]),round(center_wire[2]),round(center_tick[0]),round(center_tick[1]),round(center_tick[2]));
      return bound;
  }


  double TransformerPixelMapProducer::GetGlobalWire(const geo::WireID& wireID){
    // Get Global Wire Coordinate for TransformerCNN
    double globalWire = -9999;
    unsigned int nwires = wireReadout->Nwires(geo::PlaneID(wireID.Cryostat, 0, wireID.Plane));
    // Induction
    if (wireReadout->SignalType(wireID) == geo::kInduction) {
      auto const WireCentre = wireReadout->Wire(wireID).GetCenter();
      int temp_tpc = 0;
      if (wireID.TPC % 2 == 0) { temp_tpc = 0; }
      else { temp_tpc = 1;  }
      geo::PlaneID const p1(wireID.Cryostat, temp_tpc, wireID.Plane);
      globalWire = wireReadout->Plane(p1).WireCoordinate(WireCentre);
    }
    // Collection
    else {
       int block = wireID.TPC / 4;
       globalWire = (double)( ((int)nwires*block) + wireID.Wire );
    }
    return round(globalWire);

  }

  std::vector<unsigned int> getUniques(std::vector<unsigned int> coll)
  {
    std::vector<unsigned int> uniques;
    for (unsigned int tpc : coll) 
    {
      if (std::find(uniques.begin(), uniques.end(), tpc) == uniques.end())
          uniques.push_back(tpc);
    }
   return uniques;
  }

  void TransformerPixelMapProducer::ShiftGlobalWire(std::vector< art::Ptr< recob::Hit > > const& cluster)
  {
    // find hits passing through an APA
    std::vector<unsigned int> list_TPCs;
    std::map<unsigned int, double> map_mean_U, map_mean_V, map_Twei_U, map_Twei_V;
    int nhits_z[24] = {0};
    for (unsigned int iHit = 0; iHit < cluster.size(); ++iHit)
    {
	geo::WireID wireid = cluster[iHit]->WireID();
        if (wireid.Plane == 2){
		nhits_z[wireid.TPC] += 1;
	}	
	//if (cluster[iHit]->PeakTime() < 30){ // select 30 ticks
	if (cluster[iHit]->PeakTime() < 150){ // select 150 ticks
		double Twei = cluster[iHit]->PeakTime()+1.; //
		list_TPCs.push_back(wireid.TPC);
		if (wireid.Plane == 0){
		    map_mean_U[wireid.TPC] += GetGlobalWire(wireid)/Twei;
		    map_Twei_U[wireid.TPC] += 1.0/Twei;
		}
		if (wireid.Plane == 1){
		    map_mean_V[wireid.TPC] += GetGlobalWire(wireid)/Twei;
		    map_Twei_V[wireid.TPC] += 1.0/Twei;
		}
	}
    }
    std::vector<unsigned int> unique_TPCs = getUniques(list_TPCs);
    for (unsigned int tmp_idx : unique_TPCs)
    {
      if (map_Twei_U[tmp_idx]>0) map_mean_U[tmp_idx] /= map_Twei_U[tmp_idx];
      if (map_Twei_V[tmp_idx]>0) map_mean_V[tmp_idx] /= map_Twei_V[tmp_idx];
    }
    // obtain offsets for U and V planes
    if (unique_TPCs.size() > 1){
      // shift global wire for an event passing through an APA
      int maxsum = 0;
      unsigned int nmax_tpcs[2] = {10000,10000};
      for (unsigned int ii = 0; ii < unique_TPCs.size()-1; ++ii)
      {
        unsigned int tpc1 = unique_TPCs[ii];
        unsigned int tpc2 = unique_TPCs[ii+1];
        if (tpc1/4 != tpc2/4) continue;
        if (tpc1%2 == tpc2%2) continue;
	int sum = nhits_z[tpc1] + nhits_z[tpc2];
	if (sum > maxsum){
		maxsum = sum;
		nmax_tpcs[0] = tpc1;
		nmax_tpcs[1] = tpc2;
	}
      }
      if (nmax_tpcs[0] < 10000 && nmax_tpcs[1] < 10000){
        double offsetU = 0, offsetV = 0;
        if (map_Twei_U[nmax_tpcs[0]] > 0 && map_Twei_U[nmax_tpcs[1]] > 0) offsetU = map_mean_U[nmax_tpcs[0]] - map_mean_U[nmax_tpcs[1]];
        if (map_Twei_V[nmax_tpcs[0]] > 0 && map_Twei_V[nmax_tpcs[1]] > 0) offsetV = map_mean_V[nmax_tpcs[0]] - map_mean_V[nmax_tpcs[1]];
        fOffset[0] = round(offsetU);
        fOffset[1] = round(offsetV);
      } // end of nmax_tpc
    } // ned of unique_TPCs
  }

  void TransformerPixelMapProducer::GetDUNEGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                                                 const geo::WireID& wireID, double localTDC,
                                                 unsigned int& globalWire, unsigned int& globalPlane, 
                                                 double& globalTDC){
    unsigned int localWire = wireID.Wire;
    unsigned int plane = wireID.Plane;
    unsigned int tpc   = wireID.TPC;

    unsigned int nWiresTPC = 400;
    unsigned int wireGap = 4;
    auto const& tpcg = geom->TPC(geo::TPCID{0, tpc});
    double driftLen = tpcg.DriftDistance();
    double apaLen = tpcg.Width() - tpcg.ActiveWidth();
    double driftVel = detProp.DriftVelocity();
    //unsigned int drift_size = (driftLen / driftVel) * 2;
    //unsigned int apa_size   = 4*(apaLen / driftVel) * 2;
    double drift_size = (driftLen / driftVel) * 2;
    double apa_size   = 4*(apaLen / driftVel) * 2;

    globalWire = 0;
    globalPlane = 0;
    // Collection plane has more wires
    if(plane == 2){
        nWiresTPC=480;
        wireGap = 5;
        globalPlane = 2;
    }
    bool includeZGap = true;
    if(includeZGap) nWiresTPC += wireGap;

  // Workspace geometry has two drift regions
  //                  |-----|-----| /  /
  //      y ^         |  3  |  2  |/  /
  //        | -| z    |-----|-----|  /
  //        | /       |  1  |  0  | /
  //  x <---|/        |-----|-----|/
  //
  
    int tpcMod4 = tpc%4;
    int offset = 0;
    // Induction views depend on the drift direction
    if(plane < 2){
        // For TPCs 0 and 3 keep U and V as defined.
        if(tpcMod4 == 0 || tpcMod4 == 3){
          globalPlane = plane;
      // But reverse the TDCs
        }
    // For TPCs 1 and 2, swap U and V.
        else{
            if(plane == 0) globalPlane = 1;
            else globalPlane = 0;
        }
    }
    if(globalPlane != 1){
        globalWire += (tpc/4)*nWiresTPC + (tpcMod4>1)*offset + localWire;
    }
    else{
        globalWire += ((23-tpc)/4)*nWiresTPC + (tpcMod4>1)*offset + localWire;
    }
    
    if(tpcMod4 == 0 || tpcMod4 == 2){
        //globalTDC = (float)( drift_size - (double)localTDC );
        globalTDC =  drift_size - localTDC;
    }
    else {
        //globalTDC = (float)( (double)localTDC + drift_size + apa_size );
        globalTDC = localTDC + drift_size + apa_size;
    }


  } // end of GetDUNEGlobalWireTDC

}
