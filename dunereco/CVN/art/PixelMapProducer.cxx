////////////////////////////////////////////////////////////////////////
/// \file    PixelMapProducer.h
/// \brief   PixelMapProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
//
//  Modifications to allow unwrapped collection view
//   - Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>
#include <numeric>

#include "dunereco/CVN/art/PixelMapProducer.h"
#include "dunereco/CVN/func/AssignLabels.h"
#include "TVector2.h"
#include "TH2D.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace cvn
{

  PixelMapProducer::PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes):
    fNWire(nWire),
    fNTdc(nTdc),
    fTRes(tRes),
    fUnwrapped(2),
    fProtoDUNE(false)
  {

    fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos)
      _cacheIntercepts();
  }

  PixelMapProducer::PixelMapProducer()
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos)
      _cacheIntercepts();
  }

  PixelMap PixelMapProducer::CreateMap(detinfo::DetectorPropertiesData const& detProp,
                                       const std::vector< art::Ptr< recob::Hit > >& cluster)
  {
    std::vector<const recob::Hit*> newCluster;
    for(const art::Ptr<recob::Hit> hit : cluster){
      newCluster.push_back(hit.get());
    }
    return CreateMap(detProp, newCluster);
  }

  PixelMap PixelMapProducer::CreateMap(detinfo::DetectorPropertiesData const& detProp,
                                       const std::vector<const recob::Hit* >& cluster)
  {
    if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos){
      fVDPlane0.clear();
      fVDPlane1.clear();
      _cacheIntercepts();
    }
    
    Boundary bound = DefineBoundary(detProp, cluster);
    return CreateMapGivenBoundary(detProp, cluster, bound);
  }

  PixelMap PixelMapProducer::CreateMapGivenBoundary(detinfo::DetectorPropertiesData const& detProp,
                                                    const std::vector<const recob::Hit*>& cluster,
      const Boundary& bound)
  {

    PixelMap pm(fNWire, fNTdc, bound);

    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {

      geo::WireID wireid     = cluster[iHit]->WireID();
      double temptdc         = cluster[iHit]->PeakTime();
      unsigned int tempWire  = wireid.Wire;
      unsigned int tempPlane = wireid.Plane;
      // Leigh: Simple modification to unwrap the collection view wire plane
      if(!fProtoDUNE){
        if(fUnwrapped == 1){
          // Jeremy: Autodetect geometry for DUNE 10kt module. Is this a bad idea??
          if (fGeometry->DetectorName() == "dune10kt_v1") {
            if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue; // Skip dummy TPCs in 10kt module
            GetDUNE10ktGlobalWireTDC(detProp, wireid.Wire,cluster[iHit]->PeakTime(),
              wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
          }
          else if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos){
            GetDUNEVertDrift3ViewGlobalWire(wireid.Wire, wireid.Plane,wireid.TPC,tempWire,tempPlane);
          }
          // Default to 1x2x6. Should probably specifically name this function as such
          else {
            GetDUNEGlobalWireTDC(detProp, wireid.Wire,cluster[iHit]->PeakTime(),
              wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
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
      const double pe  = cluster[iHit]->Integral();
      const unsigned int wire = tempWire;
      const unsigned int wirePlane = tempPlane;
      const double tdc = temptdc;
      pm.Add(wire, tdc, wirePlane, pe);

    }
    return pm;
  }

  std::ostream& operator<<(std::ostream& os, const PixelMapProducer& p)
  {
    os << "PixelMapProducer: "
      << p.NTdc()  <<" tdcs X  " <<  p.NWire() << " wires";
    return os;
  }

  double PixelMapProducer::_getIntercept(geo::WireID wireid) const
  {
    const geo::WireGeo* pwire = fGeometry->WirePtr(wireid);
    geo::Point_t center = pwire->GetCenter();
    double slope = 0.;
    if(!pwire->isVertical()) slope = pwire->TanThetaZ();
    
    double intercept = center.Y() - slope*center.Z();
    if(wireid.Plane == 2) intercept = 0.;
    
    return intercept;
  }

  void PixelMapProducer::_cacheIntercepts(){
   
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

  Boundary PixelMapProducer::DefineBoundary(detinfo::DetectorPropertiesData const& detProp,
                                            const std::vector< const recob::Hit*>& cluster)
  {

    std::vector<double> time_0;
    std::vector<double> time_1;
    std::vector<double> time_2;
    //std::vector<int> wire;

    std::vector<int> wire_0;
    std::vector<int> wire_1;
    std::vector<int> wire_2;
    

    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
      geo::WireID wireid = cluster[iHit]->WireID();
      
      unsigned int globalWire = wireid.Wire;
      unsigned int globalPlane = wireid.Plane;
      double globalTime = cluster[iHit]->PeakTime();
      if(!fProtoDUNE){
        if(fUnwrapped == 1) {
          if (fGeometry->DetectorName() == "dune10kt_v1") {
            if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue; // Skip dummy TPCs in 10kt module
            GetDUNE10ktGlobalWireTDC(detProp, wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
          }
          else if (fGeometry->DetectorName().find("dunevd10kt_3view") != std::string::npos){
            GetDUNEVertDrift3ViewGlobalWire(wireid.Wire, wireid.Plane,wireid.TPC,globalWire,globalPlane);
          }
          else {
            GetDUNEGlobalWireTDC(detProp, wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
          }
        }
        else if(fUnwrapped == 2){
          GetDUNEGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,globalWire,globalPlane);
        }
      }
      else{
        GetProtoDUNEGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,globalWire,globalPlane);
      }

      if(globalPlane==0){
        time_0.push_back(globalTime);
        wire_0.push_back(globalWire);
      }
      if(globalPlane==1){
        time_1.push_back(globalTime);
        wire_1.push_back(globalWire);
      }
      if(globalPlane==2){
        time_2.push_back(globalTime);
        wire_2.push_back(globalWire);
      }
    }

    double tsum_0 = std::accumulate(time_0.begin(), time_0.end(), 0.0);
    double tmean_0 = tsum_0 / time_0.size();

    double tsum_1 = std::accumulate(time_1.begin(), time_1.end(), 0.0);
    double tmean_1 = tsum_1 / time_1.size();

    double tsum_2 = std::accumulate(time_2.begin(), time_2.end(), 0.0);
    double tmean_2 = tsum_2 / time_2.size();

    std::cout << "Boundary wire vector sizes: " << wire_0.size() << ", " << wire_1.size() << ", " << wire_2.size() << std::endl;

    auto minwireelement_0= std::min_element(wire_0.begin(), wire_0.end());
    if(wire_0.size() > 0) std::cout<<"minwire 0: "<<*minwireelement_0<<std::endl;
    auto minwireelement_1= std::min_element(wire_1.begin(), wire_1.end());
    if(wire_1.size() > 0) std::cout<<"minwire 1: "<<*minwireelement_1<<std::endl;
    auto minwireelement_2= std::min_element(wire_2.begin(), wire_2.end());
    if(wire_2.size() > 0) std::cout<<"minwire 2: "<<*minwireelement_2<<std::endl;

    //auto maxwireelement_0= std::max_element(wire_0.begin(), wire_0.end());
    //std::cout<<"maxwire 0: "<<*maxwireelement_0<<std::endl;
    //auto maxwireelement_1= std::max_element(wire_1.begin(), wire_1.end());
    //std::cout<<"maxwire 1: "<<*maxwireelement_1<<std::endl;
    //auto maxwireelement_2= std::max_element(wire_2.begin(), wire_2.end());
    //std::cout<<"maxwire 2: "<<*maxwireelement_2<<std::endl;


    int minwire_0 = 0;
    int minwire_1 = 0;
    int minwire_2 = 0;
    if(wire_0.size() > 0) minwire_0 = *minwireelement_0-1;
    if(wire_1.size() > 0) minwire_1 = *minwireelement_1-1;
    if(wire_2.size() > 0) minwire_2 = *minwireelement_2-1;

    Boundary bound(fNWire,fTRes,minwire_0,minwire_1,minwire_2,tmean_0,tmean_1,tmean_2);

    return bound;
  }

  void PixelMapProducer::GetDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
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
  void PixelMapProducer::GetDUNEGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
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
//    int dir = fGeometry->TPC(tpc,0).DetectDriftDirection();

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
    // if (plane < 2 and tpc%2 == 1) globalPlane = !plane;
    if (plane < 2 and tpcMod4 > 0 and tpcMod4 < 3) globalPlane = !plane;
    else globalPlane = plane;

    // int offset = 752; // Offset between upper and lower modules in induction views, from Robert & Dorota's code
    int offset = 0; // Offset between upper and lower modules in induction views, from Robert & Dorota's code
    // Second induction plane gets offset from the back of the TPC
    // if (globalPlane != 1) globalWire += (tpc/4)*nWiresTPC;
    // else globalWire += ((23-tpc)/4)*nWiresTPC;
    if (globalPlane != 1) globalWire += (tpc/4)*nWiresTPC + (tpcMod4 > 1)*offset + localWire;
    else globalWire += ((23-tpc)/4)*nWiresTPC + (tpcMod4 > 1)*offset + localWire;
    // Reverse wires and add offset for upper modules in induction views
    // Nitish : what's the difference between Nwires here and nWiresTPC?
    // if (tpcMod4 > 1 and globalPlane < 2) globalWire += fGeometry->Nwires(globalPlane, tpc, 0) + offset - localWire;
    // else globalWire += localWire;

    if(tpcMod4 == 0 || tpcMod4 == 2){
      globalTDC = drift_size - localTDC;
    }
    else{
      globalTDC = localTDC + drift_size + apa_size;
    }
  }

  void PixelMapProducer::GetDUNE10ktGlobalWireTDC(detinfo::DetectorPropertiesData const& detProp,
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
    if (tpc_xy > 3 and globalPlane < 2) globalWire += fGeometry->Nwires(geo::PlaneID{tpcgeom.ID(), globalPlane}) + offset - localWire;
    else globalWire += localWire;

    if (tpc_x % 2 == 0) globalTDC = localTDC;
    else globalTDC = (2*drift_size) - localTDC;
    if (tpc_x > 1) globalTDC += 2 * (drift_size + apa_size);

  } // function PixelMapProducer::GetDUNE10ktGlobalWireTDC

  // Special case for ProtoDUNE where we want to extract single particles to mimic CCQE interactions. The output pixel maps should be the same as the workspace
  // but we need different treatment of the wire numbering
  void PixelMapProducer::GetProtoDUNEGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
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
  
  } // function PixelMapProducer::GetProtoDUNEGlobalWire

  // Special case for ProtoDUNE where we want to extract single particles to mimic CCQE interactions. The output pixel maps should be the same as the workspace
  // but we need different treatment of the wire numbering
  void PixelMapProducer::GetProtoDUNEGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc,
    unsigned int& globalWire, double& globalTDC, unsigned int& globalPlane) const
  {
    // We can just use the existing function to get the global wire & plane
 //   GetProtoDUNEGlobalWire(localWire, plane, tpc, globalWire, globalPlane);
    GetDUNEGlobalWire(localWire, plane, tpc, globalWire, globalPlane);
    // Implement additional mirroring here?

  } // function GetProtoDUNEGlobalWireTDC
  
  void PixelMapProducer::GetDUNEVertDrift3ViewGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
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

  void PixelMapProducer::GetHitTruth(detinfo::DetectorClocksData const& clockData,
                                     art::Ptr<recob::Hit>& hit, std::vector<int>& pdgs,
    std::vector<int>& tracks, std::vector<float>& energy, std::vector<std::string>& process) {

    // BackTracker and ParticleInventory
    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi;

    // Get true particle and PDG responsible for this hit
    std::vector<sim::TrackIDE> IDEs = bt->HitToTrackIDEs(clockData, hit);
    for (sim::TrackIDE & k : IDEs) {
      tracks.push_back(k.trackID); // add track ID
      simb::MCParticle p = pi->TrackIdToParticle(k.trackID);
      process.push_back(p.Process()); // add G4 process string
    
      // Manually check to see if we have a Michel electron here
      int pdg = p.PdgCode();
      if (abs(pdg) == 11 && p.Mother() != 0) {
        int pdgParent = pi->TrackIdToParticle(p.Mother()).PdgCode();
        if (abs(pdgParent) == 13 && ((pdg>0)==(pdgParent>0))) {
          // If it's a Michel electron, tag it with an unphysical PDG code             
          if (p.Process() == "muMinusCaptureAtRest" || p.Process() == "muPlusCaptureAtRest") {
            pdg = 99;
          }
        } // If electron parent is muon
      } // If non-primary electron
      pdgs.push_back(pdg);
      energy.push_back(k.energy);
    }  // for trackIDE
  } // function PixelMapProducer::GetHitTruth

  SparsePixelMap PixelMapProducer::CreateSparseMap2D(
    detinfo::DetectorClocksData const& clockData,
    detinfo::DetectorPropertiesData const& detProp,
    std::vector< art::Ptr< recob::Hit> >& cluster,
    bool usePixelTruth) {

    // 3-dimensional coordinates (wire, time TPC) and 3 views
    SparsePixelMap map(3, 3, usePixelTruth);

    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi;
    
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit) {

      geo::WireID wireid       = cluster[iHit]->WireID();
      double globalTime        = cluster[iHit]->PeakTime();
      unsigned int globalWire  = wireid.Wire;
      unsigned int globalPlane = wireid.Plane;

      if (fGeometry->DetectorName().find("1x2x6") != std::string::npos) {
        GetDUNEGlobalWireTDC(detProp, wireid.Wire, cluster[iHit]->PeakTime(),
          wireid.Plane, wireid.TPC, globalWire, globalPlane, globalTime);
      }
      else if (fGeometry->DetectorName() == "dune10kt_v1") {
        if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue;
        GetDUNE10ktGlobalWireTDC(detProp, wireid.Wire, cluster[iHit]->PeakTime(),
          wireid.Plane, wireid.TPC, globalWire, globalPlane, globalTime);
      }
      else if (fGeometry->DetectorName().find("protodune") != std::string::npos) {
        GetProtoDUNEGlobalWireTDC(wireid.Wire, cluster[iHit]->PeakTime(),
          wireid.Plane, wireid.TPC, globalWire, globalTime, globalPlane);
      }
      else throw art::Exception(art::errors::UnimplementedFeature)
        << "Geometry " << fGeometry->DetectorName() << " not implemented "
        << "in CreateSparseMap." << std::endl;

      std::vector<float> coordinates = { (float)globalWire, (float)globalTime, (float)wireid.TPC };

      if (usePixelTruth) {
        // Get true information for this hit
        std::vector<int> pdgs, tracks;
        std::vector<float> energy;
        std::vector<std::string> process;
        GetHitTruth(clockData, cluster[iHit], pdgs, tracks, energy, process);
        map.AddHit(globalPlane, coordinates, {cluster[iHit]->Integral()}, pdgs, tracks, energy, process); 
      } // if PixelTuth 

      else {
        map.AddHit(0, coordinates, {cluster[iHit]->Integral()});
      }
    } // for iHit

    return map;

  } // function PixelMapProducer::CreateSparseMap

  SparsePixelMap PixelMapProducer::CreateSparseMap3D(
    detinfo::DetectorClocksData const& clockData,
    std::vector<art::Ptr<recob::SpacePoint>>& sp,
    std::vector<std::vector<art::Ptr<recob::Hit>>>& hit) {

    // 3D coordinates (x,y,z) and a single 3D view
    SparsePixelMap map(3, 1, true);

    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi;

    for (size_t iSP = 0; iSP < sp.size(); ++iSP) { // Loop over spacepoints

      std::vector<float> features(3, 0.); // charge on each plane
      std::vector<float> coordinates(3);
      const double *pos = sp[iSP]->XYZ();
      for (size_t p = 0; p < 3; ++p) coordinates[p] = pos[p];

      std::vector<int> pdgs, tracks;
      std::vector<float> energy;
      std::vector<std::string> process;

      for (size_t iH = 0; iH < hit[iSP].size(); ++iH) { // Loop over this spacepoint's hits
        features[hit[iSP][iH]->View()] += hit[iSP][iH]->Integral(); // Add hit integral to corresponding view's features
        GetHitTruth(clockData, hit[iSP][iH], pdgs, tracks, energy, process);
        map.AddHit(0, coordinates, features, pdgs, tracks, energy, process); 
      } // for hit iH
    } // for spacepoint iSP

    return map;

  } // function PixelMapProducer::CreateSparseMap

} // namespace cvn
