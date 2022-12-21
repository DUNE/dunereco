////////////////////////////////////////////////////////////////////////
/// \file    PixelMapSimProducer.h
/// \brief   PixelMapSimProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
//
//  Modifications to allow unwrapped collection view
//   - Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>

#include "dunereco/CVN/art/PixelMapSimProducer.h"
#include "dunereco/CVN/func/AssignLabels.h"
#include "TVector2.h"
#include "TH2D.h"
// #include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cvn
{

  PixelMapSimProducer::PixelMapSimProducer(unsigned int nWire, unsigned int nTdc, double tRes, double threshold):
    fNWire(nWire),
    fNTdc(nTdc),
    fTRes(tRes),
    fThreshold(threshold),
    fUnwrapped(2),
    fProtoDUNE(false),
    fTotHits(0)
  {

    fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos)
      _cacheIntercepts();
    
  }

  PixelMapSimProducer::PixelMapSimProducer()
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos)
      _cacheIntercepts();
  }

  PixelMap PixelMapSimProducer::CreateMap(detinfo::DetectorPropertiesData const& detProp,
                                       const std::vector< art::Ptr< sim::SimChannel > >& cluster)
  {
    std::vector<const sim::SimChannel*> newCluster;
    for(const art::Ptr<sim::SimChannel> hit : cluster){
      newCluster.push_back(hit.get());
    }
    return CreateMap(detProp, newCluster);
  }

  PixelMap PixelMapSimProducer::CreateMap(detinfo::DetectorPropertiesData const& detProp,
                                       const std::vector<const sim::SimChannel* >& cluster)
  {
    Boundary bound = DefineBoundary(detProp, cluster);
    return CreateMapGivenBoundary(detProp, cluster, bound);
  }

  PixelMap PixelMapSimProducer::CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,
                                                    const std::vector<const sim::SimChannel*>& cluster,
      const Boundary& bound)
  {

    PixelMap pm(fNWire, fNTdc, bound);
    
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
      
      const sim::SimChannel* reco_wire = cluster[iHit];
      auto& ROIs = reco_wire->TDCIDEMap();
      if(!(ROIs.size())) continue;
      // auto ROIs = reco_wire->SignalROI();

      std::vector<geo::WireID> wireids = fGeometry->ChannelToWire(reco_wire->Channel());
      if(!wireids.size()) continue;
      geo::WireID wireid = wireids[0];
      
      // if(wireids.size() > 1){
      //   for(auto iwire : wireids)
      //     if(iwire.Plane == reco_wire->View()) wireid = iwire;
      // }
      unsigned int tempWire  = wireid.Wire;
      unsigned int tempPlane = wireid.Plane;

      if(!fProtoDUNE){
        if(fUnwrapped == 1){
          if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos){
            GetDUNEVertDrift3ViewGlobalWire(wireid.Wire, wireid.Plane,wireid.TPC,tempWire,tempPlane);
          }
        }
        else if(fUnwrapped == 2){
          // Old method that has problems with the APA crossers, kept for old times' sake
          GetDUNEGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,tempWire,tempPlane);
        }
      }
      else{
        GetProtoDUNEGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,tempWire,tempPlane);
      }

      for(auto iROI = ROIs.begin(); iROI != ROIs.end(); ++iROI){
        auto& ROI = *iROI;
        auto tick = ROI.first;
        double temptdc  = (double)tick;
        double charge =  0.005*reco_wire->Charge(tick); 
        if(!(charge > fThreshold)) continue;   
        // Leigh: Simple modification to unwrap the collection view wire plane
        if(!fProtoDUNE){
          if(fUnwrapped == 1){
            // Jeremy: Autodetect geometry for DUNE 10kt module. Is this a bad idea??
            if (fGeometry->DetectorName() == "dune10kt_v1") {
              if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue; // Skip dummy TPCs in 10kt module
              GetDUNE10ktGlobalWireTDC(detProp, wireid.Wire,(double)tick,
                wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
            }
            else if (fGeometry->DetectorName().find("dunevd10kt_3view") == std::string::npos){
              GetDUNEGlobalWireTDC(detProp, wireid.Wire,(double)tick,
                wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
            }
          }
        }
        const double pe  = charge;
        const unsigned int wire = tempWire;
        const unsigned int wirePlane = tempPlane;
        const double tdc = temptdc;
        pm.Add(wire, tdc, wirePlane, pe);
      }

    }
    pm.SetTotHits(fTotHits);
    return pm;
  }

  std::ostream& operator<<(std::ostream& os, const PixelMapSimProducer& p)
  {
    os << "PixelMapSimProducer: "
      << p.NTdc()  <<" tdcs X  " <<  p.NWire() << " wires";
    return os;
  }

  double PixelMapSimProducer::_getIntercept(geo::WireID wireid) const
  {
    const geo::WireGeo* pwire = fGeometry->WirePtr(wireid);
    geo::Point_t center = pwire->GetCenter();
    double slope = 0.;
    if(!pwire->isVertical()) slope = pwire->TanThetaZ();
    
    double intercept = center.Y() - slope*center.Z();
    if(wireid.Plane == 2) intercept = 0.;
    
    return intercept;
  }

  void PixelMapSimProducer::_cacheIntercepts(){
   
    // double spacing = 0.847;
    for(int plane = 0; plane < 2; plane++){
      
      int nCRM_row = 6;
      int nCRM_col = 6;
      bool is8x6 = fGeometry->DetectorName().find("8x6") != std::string::npos;
      if(is8x6)
        nCRM_col = 8;
      
      geo::WireID wstart = geo::WireID(0, 0, plane,0);
      geo::WireID wstartplus1 = geo::WireID(0, 0, plane, 1);
      double wstart_intercept = _getIntercept(wstart);
      double wstartplus1_intercept = _getIntercept(wstartplus1);
      if(plane == 0)
        fSpacing0 = std::abs(wstartplus1_intercept - wstart_intercept);
      else
        fSpacing1 = std::abs(wstartplus1_intercept - wstart_intercept);

      bool is3view30deg = fGeometry->DetectorName().find("30deg") != std::string::npos;
      
      for(int diag_tpc = 0; diag_tpc < nCRM_row; diag_tpc++){
        
        int tpc = (plane == 0 || !is3view30deg) ? (nCRM_col+1)*diag_tpc : (nCRM_col-1)*(nCRM_row-diag_tpc);
        geo::PlaneID const planeID(0, tpc, plane);
        unsigned int nWiresTPC = fGeometry->Nwires(planeID);
        
        geo::WireID start = geo::WireID(planeID, 0);
        geo::WireID end = geo::WireID(planeID, nWiresTPC-1);

        double start_intercept = _getIntercept(start);
        double end_intercept = _getIntercept(end);
        if(plane == 0){
          fVDPlane0.push_back(start_intercept);
          fVDPlane0.push_back(end_intercept);
        }
        else if (is3view30deg){
          fVDPlane1.push_back(end_intercept);
          fVDPlane1.push_back(start_intercept);
        }
        else{
          fVDPlane1.push_back(start_intercept);
          fVDPlane1.push_back(end_intercept);
        }
      }

    }
  }

  Boundary PixelMapSimProducer::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                            const std::vector< const sim::SimChannel*>& cluster)
  {

    std::vector<int> wire_0;
    std::vector<int> wire_1;
    std::vector<int> wire_2;
    
    std::vector<double> twire_0;
    std::vector<double> twire_1;
    std::vector<double> twire_2;
   
    double tsum_0 = 0., tsum_1 = 0., tsum_2 = 0.;
    int total_t0 = 0, total_t1 = 0, total_t2 = 0;

    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
      const sim::SimChannel* reco_wire = cluster[iHit];
      auto& ROIs = reco_wire->TDCIDEMap();
      if(!(ROIs.size())) continue;

      std::vector<geo::WireID> wireids = fGeometry->ChannelToWire(reco_wire->Channel());
      if(!wireids.size()) continue;
      geo::WireID wireid = wireids[0];
      
      // if(wireids.size() > 1){
      //   for(auto iwire : wireids)
      //     if(iwire.Plane == reco_wire->View()) wireid = iwire;
      // }
      unsigned int globalWire  = wireid.Wire;
      unsigned int globalPlane = wireid.Plane;
     
      if(!fProtoDUNE){
        if(fUnwrapped == 1){
          if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos){
            GetDUNEVertDrift3ViewGlobalWire(wireid.Wire, wireid.Plane,wireid.TPC,globalWire,globalPlane);
          }
        }
        else if(fUnwrapped == 2){
          // Old method that has problems with the APA crossers, kept for old times' sake
          GetDUNEGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,globalWire,globalPlane);
        }
      }
      else{
        GetProtoDUNEGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,globalWire,globalPlane);
      }

      for(auto iROI = ROIs.begin(); iROI != ROIs.end(); ++iROI){
        auto& ROI = *iROI;
        auto tick = ROI.first;
        // for(int tick = ROI.begin_index(); tick < (int)ROI.end_index(); tick++){
          
        double globalTime  = (double)tick;
        double charge = 0.005*reco_wire->Charge(tick);
        if(!(charge > fThreshold)) continue;  
        // Leigh: Simple modification to unwrap the collection view wire plane
        if(!fProtoDUNE){
          if(fUnwrapped == 1){
            // Jeremy: Autodetect geometry for DUNE 10kt module. Is this a bad idea??
            if (fGeometry->DetectorName() == "dune10kt_v1") {
              if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue; // Skip dummy TPCs in 10kt module
              GetDUNE10ktGlobalWireTDC(detProp, wireid.Wire,(double)tick,
                wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
            }
            else if (fGeometry->DetectorName().find("dunevd10kt_3view") == std::string::npos){
              GetDUNEGlobalWireTDC(detProp, wireid.Wire,(double)tick,
                wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
            }
          }
        }

        if(globalPlane==0){
          tsum_0 += globalTime;
          total_t0++;
          wire_0.push_back(globalWire);
          twire_0.push_back((double)tick);
        }
        if(globalPlane==1){
          tsum_1 += globalTime;
          total_t1++;
          wire_1.push_back(globalWire);
          twire_1.push_back((double)tick);
        }
        if(globalPlane==2){
          tsum_2 += globalTime;
          total_t2++;
          wire_2.push_back(globalWire);
          twire_2.push_back((double)tick);
        }
        
        // }
      }
   
      // if(!none_threshold){
      //   if(globalPlane==0){
      //     wire_0.push_back(globalWire);
      //     twire_0.push_back((double)min_tick);
      //   }
      //   if(globalPlane==1){
      //     wire_1.push_back(globalWire);
      //     twire_1.push_back((double)min_tick);
      //   }
      //   if(globalPlane==2){
      //     wire_2.push_back(globalWire);
      //     twire_2.push_back((double)min_tick);
      //   }
      // }

    }
    double tmean_0 = tsum_0/total_t0; 
    double tmean_1 = tsum_1/total_t1; 
    double tmean_2 = tsum_2/total_t2; 
    

    //auto maxwireelement_0= std::max_element(wire_0.begin(), wire_0.end());
    //std::cout<<"maxwire 0: "<<*maxwireelement_0<<std::endl;
    //auto maxwireelement_1= std::max_element(wire_1.begin(), wire_1.end());
    //std::cout<<"maxwire 1: "<<*maxwireelement_1<<std::endl;
    //auto maxwireelement_2= std::max_element(wire_2.begin(), wire_2.end());
    //std::cout<<"maxwire 2: "<<*maxwireelement_2<<std::endl;


    std::vector<int> bwire_0;
    std::vector<int> bwire_1;
    std::vector<int> bwire_2;
    for(int i = 0; i < (int)wire_0.size(); i++){
      double t = twire_0[i];
      if(std::abs(t-tmean_0) < (double)fTRes)
        bwire_0.push_back(wire_0[i]);
    }
    for(int i = 0; i < (int)wire_1.size(); i++){
      double t = twire_1[i];
      if(std::abs(t-tmean_1) < (double)fTRes)
        bwire_1.push_back(wire_1[i]);
    }
    for(int i = 0; i < (int)wire_2.size(); i++){
      double t = twire_2[i];
      if(std::abs(t-tmean_2) < (double)fTRes)
        bwire_2.push_back(wire_2[i]);
    }
    
    std::cout << "Boundary wire vector sizes: " << bwire_0.size() << ", " << bwire_1.size() << ", " << bwire_2.size() << std::endl;
   
    int minwire_0 = 0;
    int minwire_1 = 0;
    int minwire_2 = 0;
    auto minwireelement_0 = std::min_element(bwire_0.begin(), bwire_0.end());
    auto minwireelement_1 = std::min_element(bwire_1.begin(), bwire_1.end());
    auto minwireelement_2 = std::min_element(bwire_2.begin(), bwire_2.end());
    
    if(bwire_0.size() > 0) { minwire_0 = *minwireelement_0-1; std::cout<<"minwire 0: "<<(*minwireelement_0) <<std::endl;}
    if(bwire_1.size() > 0) { minwire_1 = *minwireelement_1-1; std::cout<<"minwire 1: "<<(*minwireelement_1) <<std::endl;}
    if(bwire_2.size() > 0) { minwire_2 = *minwireelement_2-1; std::cout<<"minwire 2: "<<(*minwireelement_2) <<std::endl;}
    
    fTotHits = bwire_0.size() + bwire_1.size() + bwire_2.size();

    Boundary bound(fNWire,fTRes,minwire_0,minwire_1,minwire_2,tmean_0,tmean_1,tmean_2);

    return bound;
  }

  void PixelMapSimProducer::GetDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
  {
    unsigned int nWiresTPC = 400;

    globalWire = localWire;
    globalPlane = 0;

    // Collection plane has more wires
    if(plane == 2){
      nWiresTPC=480;
      globalPlane = 2;
    }

    // Workspace geometry has two drift regions
    //                  |-----|-----| /  /
    //      y ^         |  3  |  2  |/  /
    //        | -| z    |-----|-----|  /
    //        | /       |  1  |  0  | /
    //  x <---|/        |-----|-----|/
    //

    int tpcMod4 = tpc%4;

    // Induction views depend on the drift direction
    if(plane < 2){
      // For drift in negative x direction keep U and V as defined.
      if(tpcMod4 == 0 || tpcMod4 == 3){
        globalPlane = plane;
      }
      // For drift in positive x direction, swap U and V.
      else{
        if(plane == 0) globalPlane = 1;
        else globalPlane = 0;
      }
    }

    if(globalPlane != 1){
      globalWire += (tpc/4)*nWiresTPC;
    }
    else{
      globalWire += ((23-tpc)/4)*nWiresTPC;
    }

  }

  // Based on Robert's code in adcutils
  void PixelMapSimProducer::GetDUNEGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                                              unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                             unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const
  {

    unsigned int nWiresTPC = 400;
    unsigned int wireGap = 4;
    auto const& tpcgeom = fGeometry->TPC(geo::TPCID{0, tpc});
    double driftLen = tpcgeom.DriftDistance();
    double apaLen = tpcgeom.Width() - tpcgeom.ActiveWidth();
    double driftVel = detProp.DriftVelocity();
    unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC 
    unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC

    globalWire = 0;
    globalPlane = 0;
//    int dir = tpcgeom.DetectDriftDirection();

    // Collection plane has more wires
    if(plane == 2){
      nWiresTPC = 480;
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
    // Induction views depend on the drift direction
    if (plane < 2 and tpc%2 == 1) globalPlane = !plane;
    else globalPlane = plane;

    int offset = 752; // Offset between upper and lower modules in induction views, from Robert & Dorota's code
    // Second induction plane gets offset from the back of the TPC
    if (globalPlane != 1) globalWire += (tpc/4)*nWiresTPC;
    else globalWire += ((23-tpc)/4)*nWiresTPC;
    // Reverse wires and add offset for upper modules in induction views
    // Nitish : what's the difference between Nwires here and nWiresTPC?
    if (tpcMod4 > 1 and globalPlane < 2) globalWire += fGeometry->Nwires(geo::PlaneID{0, tpc, globalPlane}) + offset - localWire;
    else globalWire += localWire;

    if(tpcMod4 == 0 || tpcMod4 == 2){
      globalTDC = drift_size - localTDC;
    }
    else{
      globalTDC = localTDC + drift_size + apa_size;
    }
  }

  void PixelMapSimProducer::GetDUNE10ktGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
                                                  unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
                                                  unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const
  {
    unsigned int nWiresTPC = 400;
    unsigned int wireGap = 4;
    auto const& tpcgeom = fGeometry->TPC(geo::TPCID{0, tpc});
    double driftLen = tpcgeom.DriftDistance();
    double apaLen = tpcgeom.Width() - tpcgeom.ActiveWidth();
    double driftVel = detProp.DriftVelocity();
    unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC 
    unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC


    globalWire = 0;
    globalPlane = 0;

    // Collection plane has more wires
    if(plane == 2){
      nWiresTPC = 480;
      wireGap = 5;
      globalPlane = 2;
    }

    bool includeZGap = true;
    if(includeZGap) nWiresTPC += wireGap;

    // 10kt has four real TPCs and two dummies in each slice
    //
    //                 |--|-----|-----|-----|-----|--| /  /
    //      y ^        |11| 10  |  9  |  8  |  7  | 6|/  /
    //        | -| z   |--|-----|-----|-----|-----|--|  /
    //        | /      | 5|  4  |  3  |  2  |  1  | 0| /
    //  x <---|/       |--|-----|-----|-----|-----|--|/
    //                     ^  wires  ^ ^  wires  ^
    //
    // We already filtered out the dummies, so we can assume 0->3 as follows:
    //
    //                 |-----|-----|-----|-----| /  /
    //      y ^        |  7  |  6  |  5  |  4  |/  /
    //        | -| z   |-----|-----|-----|-----|  /
    //        | /      |  3  |  2  |  1  |  0  | /
    //  x <---|/       |-----|-----|-----|-----|/
    //                  ^  wires  ^ ^  wires  ^
    //

    size_t tpc_x = (tpc%6) - 1;   // x coordinate in 0->4 range
    size_t tpc_xy = (tpc%12) - 1; // xy coordinate as 0->3 & 6->9 (converted from 1->4, 7->10)
    if (tpc_xy > 3) tpc_xy -= 2;  // now subtract 2 so it's in the 0->7 range

    // Induction views depend on the drift direction
    if (plane < 2 and tpc%2 == 1) globalPlane = !plane;
    else globalPlane = plane;

    int offset = 752; // Offset between upper and lower modules in induction views, from Robert & Dorota's code
    // Second induction plane gets offset from the back of the TPC
    if (globalPlane != 1) globalWire += (tpc/12)*nWiresTPC;
    else globalWire += ((300-tpc)/12)*nWiresTPC;
    // Reverse wires and add offset for upper modules in induction views
    if (tpc_xy > 3 and globalPlane < 2) globalWire += fGeometry->Nwires(geo::PlaneID{0, tpc, globalPlane}) + offset - localWire;
    else globalWire += localWire;

    if (tpc_x % 2 == 0) globalTDC = localTDC;
    else globalTDC = (2*drift_size) - localTDC;
    if (tpc_x > 1) globalTDC += 2 * (drift_size + apa_size);

  } // function PixelMapSimProducer::GetDUNE10ktGlobalWireTDC

  // Special case for ProtoDUNE where we want to extract single particles to mimic CCQE interactions. The output pixel maps should be the same as the workspace
  // but we need different treatment of the wire numbering
  void PixelMapSimProducer::GetProtoDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
  { 
    unsigned int nWiresTPC = 400;
    
    globalWire = localWire;
    globalPlane = 0;
    
    // Collection plane has more wires
    if(plane == 2){
      nWiresTPC=480;
      globalPlane = 2;
    }
    
    // ProtoDUNE has a central CPA so times are fine
    // It (annoyingly) has two dummy TPCs on the sides
    //                  
    //      y ^       |-|-----|-----|-|   /
    //        | -| z  | |     |     | |  /
    //        | /     |3|  2  |  1  |0| /
    //  x <---|/      |-|-----|-----|-|/
    //
    
    int tpcMod4 = tpc%4;
    // tpcMod4: 1 for -ve drift, 2 for +ve drift
    // Induction views depend on the drift direction
    if(plane < 2){
      // For drift in negative x direction keep U and V as defined.
      if(tpcMod4 == 1){
        globalPlane = plane;
      }
      // For drift in positive x direction, swap U and V.
      else{
        if(plane == 0) globalPlane = 1;
        else globalPlane = 0;
      }
    }
    
    if(globalPlane != 1){
      globalWire += (tpc/4)*nWiresTPC;
    }
    else{
      globalWire += ((12-tpc)/4)*nWiresTPC;
    }
  
  } // function PixelMapSimProducer::GetProtoDUNEGlobalWire

  // Special case for ProtoDUNE where we want to extract single particles to mimic CCQE interactions. The output pixel maps should be the same as the workspace
  // but we need different treatment of the wire numbering
  void PixelMapSimProducer::GetProtoDUNEGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
    unsigned int& globalWire, double& globalTDC, unsigned int& globalPlane) const
  {
    // We can just use the existing function to get the global wire & plane
 //   GetProtoDUNEGlobalWire(localWire, plane, tpc, globalWire, globalPlane);
    GetDUNEGlobalWire(localWire, plane, tpc, globalWire, globalPlane);
    // Implement additional mirroring here?

  } // function GetProtoDUNEGlobalWireTDC
  
  void PixelMapSimProducer::GetDUNEVertDrift3ViewGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
  {
    // Preliminary function for VD Geometries
    
    int nCRM_row = 6;
    int nCRM_col = 6;
    bool is8x6 = fGeometry->DetectorName().find("8x6") != std::string::npos;
    if(is8x6)
      nCRM_col = 8;
    // spacing between y-intercepts of parallel wires in a given plane. 
    double spacing = 0.847; 
    
    globalPlane = plane;
    geo::PlaneID const planeID{0, tpc, globalPlane};
    unsigned int nWiresTPC = fGeometry->Nwires(planeID);
    bool is3view30deg = fGeometry->DetectorName().find("30deg") != std::string::npos;
    
    if(globalPlane < 2){
      
      geo::WireID wire_id = geo::WireID(planeID, localWire);
      double wire_intercept = _getIntercept(wire_id);
   
      double low_bound = 0., upper_bound = 0.; 
      int start, end, diag_tpc;
      // get wires on diagonal CRMs and their intercepts which bound the current wire's intercept 
      if(globalPlane == 0){
        start = std::lower_bound(fVDPlane0.begin(), fVDPlane0.end(), wire_intercept) - fVDPlane0.begin() - 1;
        end = std::upper_bound(fVDPlane0.begin(), fVDPlane0.end(), wire_intercept) - fVDPlane0.begin();
        
        low_bound = fVDPlane0[start];
        if(end < (int) fVDPlane0.size())
          upper_bound = fVDPlane0[end];
        diag_tpc = (start/2);
        spacing = fSpacing0;
      }
      else if (is3view30deg){
        end = std::lower_bound(fVDPlane1.begin(), fVDPlane1.end(), wire_intercept) - fVDPlane1.begin() - 1;
        start = std::upper_bound(fVDPlane1.begin(), fVDPlane1.end(), wire_intercept) - fVDPlane1.begin();
       
        if(end > 0)
          low_bound = fVDPlane1[end];
        upper_bound = fVDPlane1[start];
        diag_tpc = (nCRM_row-(end/2) - 1);
        if(end < 0) 
          diag_tpc = nCRM_row - 1;
        spacing = fSpacing1;
      }
      else{
        int tpc_y = tpc % nCRM_col;
        globalWire = localWire + tpc_y*nWiresTPC;
        
        start = std::lower_bound(fVDPlane1.begin(), fVDPlane1.end(), wire_intercept) - fVDPlane1.begin() - 1;
        end = std::upper_bound(fVDPlane1.begin(), fVDPlane1.end(), wire_intercept) - fVDPlane1.begin();
        
        low_bound = fVDPlane1[start];
        if(end < (int) fVDPlane1.size())
          upper_bound = fVDPlane1[end];
        diag_tpc = (start/2);
      }
      // if the intercept of the wire is in between two diagonal CRMs, assign it to the diagonal CRM its closest to 
      if((((start % 2)^globalPlane && is3view30deg) || (globalPlane == 0 && !is3view30deg && (start % 2 == 1))) && (diag_tpc < nCRM_row-1)){
        int diag_idx = diag_tpc + !globalPlane;
        globalWire = (wire_intercept > (low_bound+upper_bound)*0.5) ? (nWiresTPC-1)*diag_idx + !globalPlane : (nWiresTPC-1)*diag_idx + globalPlane;
      }
      // otherwise assign it to the closest wire within the same CRM
      else if(is3view30deg || (globalPlane == 0 && !is3view30deg && (start %2 == 0))){
        int diag_idx = diag_tpc;
        int offset = globalPlane ? std::round((upper_bound - wire_intercept)/spacing) : std::round((wire_intercept-low_bound)/spacing);
        globalWire = (nWiresTPC-1)*diag_idx + offset + 1;
        
      }
    }
    else{
      int tpc_z = tpc/nCRM_col;
      globalWire = localWire + tpc_z*nWiresTPC;
    }
  }

} // namespace cvn
