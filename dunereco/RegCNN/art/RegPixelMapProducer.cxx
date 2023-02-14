////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMapProducer.h
/// \brief   RegPixelMapProducer for RegCNN modified from PixelMapProducer.h
/// \author  Ilsoo Seong - iseong@uci.edu
//           Wenjie Wu - wenjieww@uci.edu
//
//  Modifications to allow unwrapped collection view
//   - Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>

#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunereco/RegCNN/art/RegPixelMapProducer.h"
#include  "TVector2.h"

namespace cnn
{

  RegPixelMapProducer::RegPixelMapProducer(unsigned int nWire, unsigned int wRes, unsigned int nTdc, double tRes, int Global, bool ProngOnly, bool ByHit):
  fNWire(nWire),
  fWRes(wRes),
  fNTdc(nTdc),
  fTRes(tRes),
  fGlobalWireMethod(Global),
  fProngOnly(ProngOnly),
  fByHit(ByHit),
  fOffset{0,0}
  {}

  RegPixelMap RegPixelMapProducer::CreateMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster, 
                                             art::FindManyP<recob::Wire> const& fmwire)
  {

    hitwireidx.clear();
    tmin_each_wire.clear();
    tmax_each_wire.clear();
    trms_max_each_wire.clear();
    fOffset[0] = 0; fOffset[1] = 0;

    RegCNNBoundary bound = DefineBoundary(detProp, cluster);

    //return CreateMapGivenBoundaryByHit(clockData, detProp, cluster, bound, fmwire, fProngOnly);
    return CreateMapGivenBoundary(clockData, detProp, cluster, bound, fmwire);

  }

  RegPixelMap RegPixelMapProducer::CreateMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster,
                                             art::FindManyP<recob::Wire> const& fmwire,
                                             const std::vector<float> &vtx)
  {
      // Create pixel maps around the vertex
      hitwireidx.clear();
      tmin_each_wire.clear();
      tmax_each_wire.clear();
      trms_max_each_wire.clear();
      fOffset[0] = 0; fOffset[1] = 0;

      RegCNNBoundary bound = DefineBoundary(detProp, cluster, vtx);

      return CreateMapGivenBoundary(clockData, detProp, cluster, bound, fmwire);
  }

  RegPixelMap RegPixelMapProducer::CreateMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster,
                                             art::FindManyP<recob::Wire> const& fmwire,
                                             art::FindManyP<recob::Track> const& fmtrkhit,
                                             const std::vector<float> &vtx)
  {
      // tag prong by track
      // Create pixel maps around the vertex
      hitwireidx.clear();
      tmin_each_wire.clear();
      tmax_each_wire.clear();
      trms_max_each_wire.clear();
      fOffset[0] = 0; fOffset[1] = 0;

      RegCNNBoundary bound = DefineBoundary(detProp, cluster, vtx);

      if (fByHit) {
        return CreateMapGivenBoundaryByHit(clockData, detProp, cluster, bound, fmwire, fmtrkhit, fProngOnly);
      } else {
        return CreateMapGivenBoundary(clockData, detProp, cluster, bound, fmwire);
      }
  }

  RegPixelMap RegPixelMapProducer::CreateMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster,
                                             art::FindManyP<recob::Wire> const& fmwire,
                                             art::FindManyP<recob::Shower> const& fmshwhit,
                                             const std::vector<float> &vtx)
  {
      // tag prong by shower
      // Create pixel maps around the vertex
      hitwireidx.clear();
      tmin_each_wire.clear();
      tmax_each_wire.clear();
      trms_max_each_wire.clear();
      fOffset[0] = 0; fOffset[1] = 0;

      RegCNNBoundary bound = DefineBoundary(detProp, cluster, vtx);

      if (fByHit) {
        return CreateMapGivenBoundaryByHit(clockData, detProp, cluster, bound, fmwire, fmshwhit, fProngOnly);
      } else {
        return CreateMapGivenBoundary(clockData, detProp, cluster, bound, fmwire);
      }
  }


  RegPixelMap RegPixelMapProducer::CreateMapGivenBoundary(detinfo::DetectorClocksData const& clockData,
          detinfo::DetectorPropertiesData const& detProp,
          std::vector< art::Ptr< recob::Hit > > const& cluster,
          const RegCNNBoundary& bound,
          art::FindManyP<recob::Wire> const& fmwire)
  {
      RegPixelMap pm(fNWire, fWRes, fNTdc, fTRes, bound, 0);

      if (!fmwire.isValid()) return pm;

      // get all raw adc of every hit wire
      for (size_t iwire = 0; iwire < hitwireidx.size(); ++iwire)
      {
          unsigned int iHit = hitwireidx[iwire];
          std::vector< art::Ptr<recob::Wire> > wireptr = fmwire.at(iHit);
          geo::WireID wireid = cluster[iHit]->WireID();

          for (size_t iwireptr = 0; iwireptr < wireptr.size(); ++iwireptr){
              std::vector<geo::WireID> wireids = geom->ChannelToWire(wireptr[iwireptr]->Channel());
              bool goodWID = false; 
              for (auto const & wid:wireids){ 
                  if (wid.Plane == wireid.Plane &&   
                          wid.Wire  == wireid.Wire &&
                          wid.TPC   == wireid.TPC &&
                          wid.Cryostat == wireid.Cryostat) goodWID = true;
              }
              if (!goodWID) continue;


              //int t0_hit = (int)( tmin_each_wire[iwire] - 3 * (trms_max_each_wire[iwire]) );
              //int t1_hit = (int)( tmax_each_wire[iwire] + 3 * (trms_max_each_wire[iwire]) );
              float hit_first_time = tmin_each_wire[iwire] - 3 * (trms_max_each_wire[iwire]);
              float hit_end_time = tmax_each_wire[iwire] + 3 * (trms_max_each_wire[iwire]);
              int t0_hit = (hit_first_time < 0) ? 0 : (int)hit_first_time; 
              int t1_hit = (hit_end_time > 4491) ? 4491 : (int)hit_end_time;

              const std::vector<float>& signal = wireptr[0]->Signal();
              unsigned int globalWire = 0;
              if (fGlobalWireMethod == 1){
                  globalWire = (unsigned int)GetGlobalWire(wireid);
                  if (wireid.TPC%2 == 1) {
                      if (wireid.Plane == 0) globalWire += fOffset[0];
                      if (wireid.Plane == 1) globalWire += fOffset[1];
                  }
              }

              unsigned int globalplane = wireid.Plane;
              for (int tt = t0_hit; tt <= t1_hit; ++tt)
              {
                  //double correctedadc = (double) signal[tt];
                  double correctedadc = ( signal[tt] * TMath::Exp( (sampling_rate(clockData) * tt) / (detProp.ElectronLifetime()*1.e3) ) );
                  int tdc = tt;
                  if (wireid.TPC%2 == 0) tdc = -tdc;
                  if (fGlobalWireMethod == 2){
                      double globaltick = double(tt);
                      GetDUNEGlobalWireTDC(detProp, wireid, (double)tt, globalWire, globalplane, globaltick);
                      //tdc = (int)round(globaltick);
                      tdc = (int)globaltick;
                      // FIXIT
                      //if (globalplane==0 && correctedadc) {
                      //  std::cout<<tt<<", "<<globalplane<<" | ";
                      //  std::cout<<globalWire<<", "<<globaltick<<", "<<tdc<<", "<<correctedadc<<std::endl;
                      //}
                  }
                  pm.Add((int)globalWire, tdc, globalplane, correctedadc, wireid.TPC, 0);

              } // end of tt
          } // end of iwireptr
      } // end of iwire

      //std::cout<< "===============> Offsets: " << fOffset[0] << " " << fOffset[1] << std::endl;
      return pm;

  }

  RegPixelMap RegPixelMapProducer::CreateMapGivenBoundaryByHit(detinfo::DetectorClocksData const& clockData,
                                                          detinfo::DetectorPropertiesData const& detProp,
                                                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                          const RegCNNBoundary& bound,
                                                          art::FindManyP<recob::Wire> const& fmwire,
                                                          art::FindManyP<recob::Track> const& fmtrkhit,
                                                          const bool& ProngOnly)
  {

      RegPixelMap pm(fNWire, fWRes, fNTdc, fTRes, bound, ProngOnly);

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

  RegPixelMap RegPixelMapProducer::CreateMapGivenBoundaryByHit(detinfo::DetectorClocksData const& clockData,
                                                          detinfo::DetectorPropertiesData const& detProp,
                                                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                          const RegCNNBoundary& bound,
                                                          art::FindManyP<recob::Wire> const& fmwire,
                                                          art::FindManyP<recob::Shower> const& fmshwhit,
                                                          const bool& ProngOnly)
  {

      RegPixelMap pm(fNWire, fWRes, fNTdc, fTRes, bound, ProngOnly);

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

      std::cout<<"FIXME: Produce pixelmap from fmshwhit and prongonly = "<<ProngOnly<<std::endl;
      pm.Finish();

      return pm;

  }


  std::ostream& operator<<(std::ostream& os, const RegPixelMapProducer& p)
  {
    os << "RegPixelMapProducer: "
       << p.NTdc()  <<" tdcs X  " <<  p.NWire() << " wires";
    return os;
  }



  RegCNNBoundary RegPixelMapProducer::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                                     std::vector< art::Ptr< recob::Hit > > const& cluster)
  {
    if (fGlobalWireMethod == 1) { ShiftGlobalWire(cluster); }

    std::vector<float> time_0;
    std::vector<float> time_1;
    std::vector<float> time_2;

    std::vector<int> wire_0;
    std::vector<int> wire_1;
    std::vector<int> wire_2;

    unsigned int temp_wire = cluster[0]->WireID().Wire;
    float temp_time_min = 1e5; //99999;
    float temp_time_max = 0;   //-99999;
    float temp_trms_max = 0;   //-99999;
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
        geo::WireID wireid = cluster[iHit]->WireID();
	//if (temp_wire != wireid.Wire || iHit == cluster.size()-1){ // lost last hit info
	if (temp_wire != wireid.Wire ){ // lost last hit info
		temp_wire = wireid.Wire;
		hitwireidx.push_back(iHit-1);
		tmin_each_wire.push_back(int(temp_time_min));
		tmax_each_wire.push_back(int(temp_time_max));
		trms_max_each_wire.push_back(temp_trms_max);
		//temp_time_min = 99999; temp_time_max = -99999, temp_trms_max = -99999;
		temp_time_min = 1e5; temp_time_max = 0; temp_trms_max = 0;
	}
	if (temp_time_min > cluster[iHit]->PeakTime()) temp_time_min = cluster[iHit]->PeakTime();
	if (temp_time_max < cluster[iHit]->PeakTime()) temp_time_max = cluster[iHit]->PeakTime();
	if (temp_trms_max < cluster[iHit]->RMS()) temp_trms_max = (float)cluster[iHit]->RMS();

	unsigned int planeid = wireid.Plane;
	double peaktime = cluster[iHit]->PeakTime() ;
	if (wireid.TPC%2 == 0) peaktime = -peaktime;

        unsigned int globalWire = 0;
        unsigned int globalplane = planeid;
        if (fGlobalWireMethod == 1){
            globalWire = GetGlobalWire(wireid);
        } else if (fGlobalWireMethod == 2 ){
            GetDUNEGlobalWireTDC(detProp, wireid, cluster[iHit]->PeakTime(), globalWire, globalplane, peaktime);
        } else {
            std::cout << "Wrong GlobalWireMethod" << std::endl;
            abort();
        }

	if (wireid.TPC%2 == 1 && fGlobalWireMethod == 1) {
	  if (wireid.Plane == 0) globalWire += fOffset[0];
	  if (wireid.Plane == 1) globalWire += fOffset[1];
	}

        if(globalplane==0){
          time_0.push_back(peaktime);
          wire_0.push_back((int)globalWire);
        }
        if(globalplane==1){
          time_1.push_back(peaktime);
          wire_1.push_back((int)globalWire);
        }
        if(globalplane==2){
          time_2.push_back(peaktime);
          wire_2.push_back((int)globalWire);
        }
    }
  

    double tsum_0 = std::accumulate(time_0.begin(), time_0.end(), 0.0);
    double tmean_0 = tsum_0 / time_0.size();

    double tsum_1 = std::accumulate(time_1.begin(), time_1.end(), 0.0);
    double tmean_1 = tsum_1 / time_1.size();

    double tsum_2 = std::accumulate(time_2.begin(), time_2.end(), 0.0);
    double tmean_2 = tsum_2 / time_2.size();


    double wiresum_0 = std::accumulate(wire_0.begin(), wire_0.end(), 0.0);
    double wiremean_0 = wiresum_0 / wire_0.size();

    double wiresum_1 = std::accumulate(wire_1.begin(), wire_1.end(), 0.0);
    double wiremean_1 = wiresum_1 / wire_1.size();

    double wiresum_2 = std::accumulate(wire_2.begin(), wire_2.end(), 0.0);
    double wiremean_2 = wiresum_2 / wire_2.size();

    //std::cout << "TDC ===> " << round(tmean_0) << " " <<   round(tmean_1) << " " << round(tmean_2) << std::endl;
    //std::cout << "Wire ==> " << round(wiremean_0) << " " << round(wiremean_1) << " " << round(wiremean_2) << std::endl;
    //std::cout << "Offset ==> " << fOffset[0] << " " << fOffset[1] << std::endl;

    //auto minwireelement_0= std::min_element(wire_0.begin(), wire_0.end());
    //std::cout<<"minwire 0: "<<*minwireelement_0<<std::endl;
    //auto minwireelement_1= std::min_element(wire_1.begin(), wire_1.end());
    //std::cout<<"minwire 1: "<<*minwireelement_1<<std::endl;
    //auto minwireelement_2= std::min_element(wire_2.begin(), wire_2.end());
    //std::cout<<"minwire 2: "<<*minwireelement_2<<std::endl;

    //auto maxwireelement_0= std::max_element(wire_0.begin(), wire_0.end());
    //std::cout<<"maxwire 0: "<<*maxwireelement_0<<std::endl;
    //auto maxwireelement_1= std::max_element(wire_1.begin(), wire_1.end());
    //std::cout<<"maxwire 1: "<<*maxwireelement_1<<std::endl;
    //auto maxwireelement_2= std::max_element(wire_2.begin(), wire_2.end());
    //std::cout<<"maxwire 2: "<<*maxwireelement_2<<std::endl;


    //int minwire_0 = *minwireelement_0-1;
    //int minwire_1 = *minwireelement_1-1;
    //int minwire_2 = *minwireelement_2-1;

    RegCNNBoundary bound(fNWire,fNTdc,fWRes,fTRes,round(wiremean_0),round(wiremean_1),round(wiremean_2),round(tmean_0),round(tmean_1),round(tmean_2));
    std::cout<<bound<<std::endl;

    return bound;
  }

  RegCNNBoundary RegPixelMapProducer::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
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
                  geo::PlaneGeo const& planegeo_temp = geom->Plane(geo::PlaneID(0, 0, iplane));
                  geo::PlaneID const planeID(rawcrys, rawtpc, iplane);
                  geo::WireID w1;
                  try { 
                      w1 = geom->NearestWireID(regvtx_loc, planeID);
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


  double RegPixelMapProducer::GetGlobalWire(const geo::WireID& wireID){
    // Get Global Wire Coordinate for RegCNN
    double globalWire = -9999;
    unsigned int nwires = geom->Nwires(geo::PlaneID(wireID.Cryostat, 0, wireID.Plane));
    // Induction
    if (geom->SignalType(wireID) == geo::kInduction) {
      auto const WireCentre = geom->WireIDToWireGeo(wireID).GetCenter();
      geo::PlaneID p1;
      int temp_tpc = 0;
      if (wireID.TPC % 2 == 0) { temp_tpc = 0; }
      else { temp_tpc = 1;  }
      p1 = geo::PlaneID(wireID.Cryostat, temp_tpc, wireID.Plane);
      globalWire = geom->WireCoordinate(WireCentre, p1);
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

  void RegPixelMapProducer::ShiftGlobalWire(std::vector< art::Ptr< recob::Hit > > const& cluster)
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

  void RegPixelMapProducer::GetDUNEGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
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
