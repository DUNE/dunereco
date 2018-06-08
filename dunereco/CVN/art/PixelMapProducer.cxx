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

#include "dune/CVN/art/PixelMapProducer.h"
#include "dune/CVN/func/AssignLabels.h"
#include "TVector2.h"
#include "TH2D.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

namespace cvn
{

  PixelMapProducer::PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes):
    fNWire(nWire),
    fNTdc(nTdc),
    fTRes(tRes),
    fUnwrapped(false)
  {

    fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }


  PixelMap PixelMapProducer::CreateMap(std::vector< art::Ptr< recob::Hit > >& cluster)
  {

    Boundary bound = DefineBoundary(cluster);

    return CreateMapGivenBoundary(cluster, bound);



  }

  PixelMap PixelMapProducer::CreateMapGivenBoundary(std::vector< art::Ptr< recob::Hit > >& cluster,
      const Boundary& bound)
  {

    PixelMap pm(fNWire, fNTdc, bound);

    TH2D *hWiresX = new TH2D("hWireX","",3000,0,3000,400,0,400);
    TH2D *hWiresY = new TH2D("hWireY","",3000,0,3000,400,0,400);
    TH2D *hWiresZ = new TH2D("hWireZ","",3000,0,3000,480,0,480);

    bool hasLeftTPC = false;
    bool hasRightTPC = false;
  
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {
      geo::WireID wireid = cluster[iHit]->WireID();
      double temptdc                   = cluster[iHit]->PeakTime();
      unsigned int tempWire                  = wireid.Wire;
      unsigned int tempPlane              = wireid.Plane;
      // Leigh: Simple modification to unwrap the collection view wire plane
      if(fUnwrapped){
        unsigned int oldWire = 0;
        unsigned int oldPlane = 0;
        GetBetterGlobalWire(wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
        GetGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,oldWire,oldPlane);
        if(wireid.Plane == 0) hWiresX->SetBinContent(tempWire+1,wireid.Wire+1,wireid.TPC);
        if(wireid.Plane == 1) hWiresY->SetBinContent(tempWire+1,wireid.Wire+1,wireid.TPC);
        if(wireid.Plane == 2) hWiresZ->SetBinContent(tempWire+1,wireid.Wire+1,wireid.TPC);
        // Test the new and old methods
//        std::cout << "Local: " << wireid.Wire << ", " << wireid.Plane << ", " << wireid.TPC << " :: New: " << tempWire << ", " << tempPlane << " :: Old " << oldWire << ", " << oldPlane << std::endl;  
//        std::cout << "TPC = " << wireid.TPC << std::endl;
      }
//      const double pe  = cluster[iHit]->Integral();
      const double pe  = static_cast<double>(wireid.TPC);
      const unsigned int wire = tempWire;
      const unsigned int wirePlane = tempPlane;
      const double tdc = temptdc;
      pm.Add(wire, tdc, wirePlane, pe);

      if((wireid.TPC % 4) == 0 || (wireid.TPC % 4) == 2) hasRightTPC = true;
      else hasLeftTPC = true;

    }

    if(hasLeftTPC && hasRightTPC) std::cout << "*** Event spans the APA" << std::endl;

    hWiresX->Write();
    hWiresY->Write();
    hWiresZ->Write();

    return pm;

  }




  std::ostream& operator<<(std::ostream& os, const PixelMapProducer& p)
  {
    os << "PixelMapProducer: "
      << p.NTdc()  <<" tdcs X  " <<  p.NWire() << " wires";
    return os;
  }



  Boundary PixelMapProducer::DefineBoundary(std::vector< art::Ptr< recob::Hit > >& cluster)
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
      //wire.push_back(wireid.Wire);
      unsigned int globalWire = wireid.Wire;
      unsigned int globalPlane = wireid.Plane;
      double globalTime = cluster[iHit]->PeakTime();
      if(fUnwrapped){
        GetBetterGlobalWire(wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
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

    auto minwireelement_0= std::min_element(wire_0.begin(), wire_0.end());
    //std::cout<<"minwire 0: "<<*minwireelement_0<<std::endl;
    auto minwireelement_1= std::min_element(wire_1.begin(), wire_1.end());
    //std::cout<<"minwire 1: "<<*minwireelement_1<<std::endl;
    auto minwireelement_2= std::min_element(wire_2.begin(), wire_2.end());
    //std::cout<<"minwire 2: "<<*minwireelement_2<<std::endl;

    //auto maxwireelement_0= std::max_element(wire_0.begin(), wire_0.end());
    //std::cout<<"maxwire 0: "<<*maxwireelement_0<<std::endl;
    //auto maxwireelement_1= std::max_element(wire_1.begin(), wire_1.end());
    //std::cout<<"maxwire 1: "<<*maxwireelement_1<<std::endl;
    //auto maxwireelement_2= std::max_element(wire_2.begin(), wire_2.end());
    //std::cout<<"maxwire 2: "<<*maxwireelement_2<<std::endl;


    int minwire_0 = *minwireelement_0-1;
    int minwire_1 = *minwireelement_1-1;
    int minwire_2 = *minwireelement_2-1;

    Boundary bound(fNWire,fTRes,minwire_0,minwire_1,minwire_2,tmean_0,tmean_1,tmean_2);

    return bound;
  }

  void PixelMapProducer::GetGlobalWire(unsigned int localWire, unsigned int plane, unsigned int tpc, unsigned int& globalWire, unsigned int& globalPlane) const
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
    //      y ^         |  2  |  3  |/  /
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
  void PixelMapProducer::GetBetterGlobalWire(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc, 
                                             unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const
  {

    unsigned int nWiresTPC = 400;
    unsigned int wireGap = 4;
    double driftLen = fGeometry->TPC(tpc,0).DriftDistance();
    double apaLen = fGeometry->TPC(tpc,0).Width() - fGeometry->TPC(tpc,0).ActiveWidth();
    double driftVel = fDetProp->DriftVelocity();
    unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC 
    unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC

//    std::cout << "MEASUREMENTS: " << drift_size << ", " << apa_size << std::endl;

    globalWire = 0;
    globalPlane = 0;
//    int dir = fGeometry->TPC(tpc,0).DetectDriftDirection();

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

      // Need to apply an offset for the upper TPCs
//      if(tpcMod4 > 1){ 
//        offset = (globalPlane == 0) ? 48 : 752;     // top-bottom offsets
//        offset = -1*dir*48;
//      }
    }

    if(globalPlane != 1){
      globalWire += (tpc/4)*nWiresTPC + (tpcMod4>1)*offset + localWire;
    }
    else{
      globalWire += ((23-tpc)/4)*nWiresTPC + (tpcMod4>1)*offset + localWire;
    }

    if(tpcMod4 == 0 || tpcMod4 == 2){
      globalTDC = drift_size - localTDC;
    }
    else{
      globalTDC = localTDC + drift_size + apa_size;
    }
/*
      //size_t weff[3] = { 404, 404, 485 }; // effective offset in wires between consecutive tpc's, incl. dead spaces, per plane
      size_t weff[3] = { 400, 400, 480 }; // effective offset in wires between consecutive tpc's... no dead space
      size_t tpcs[3] = { 2, 2, 6 };       // FD workspace dimension in #tpc's

      size_t tpc_z = tpc / (tpcs[0] * tpcs[1]);
      size_t tpc_y = (tpc / tpcs[0]) % tpcs[1];
//      size_t tpc_x = tpc % tpcs[0];

      // Always zero for now
      unsigned short cryo = 0;

      int dir = fGeometry->TPC(tpc, cryo).DetectDriftDirection();
      std::cout << "TPC = " << tpc << " :: " << "Drift direction = " << dir << std::endl;
      bool flip_w = false;     // need the flip in wire direction?
      globalPlane = plane;        // global plane
      int offset = 0;
//      int p_align = 0;

      if (plane < 2) // no wire flip nor top-bottom alignments in Z (collection) plane
      {
        // For TPCs with negative drift direction
        if ((tpc % 4 == 0) || (tpc % 4 == 3))
        //if (dir == -1)
        {
          globalPlane = (plane == 1) ? 0 : 1;     // swap U and V
        }

        if((tpc % 4) > 1){
          if(dir > 0) flip_w = (globalPlane == 1);
          else flip_w = (globalPlane == 0);
        }

//        if ((tpc % 4 == 0) || (tpc % 4 == 1)) // bottom tpcs
//        {
//          flip_w = (dir > 0) ? (globalPlane == 1) : (globalPlane == 0);
//        }
//        else                              // top tpc's
//        {
//          flip_w = (dir < 0) ? (globalPlane == 1) : (globalPlane == 0);
//        }

        offset = (plane == 0) ? 48 : 752;     // top-bottom offsets
      }

      globalWire = tpc_z * weff[plane] + tpc_y*offset + localWire;
      if(flip_w){
          globalWire = tpc_z*weff[plane] + tpc_y*offset + weff[plane] - localWire - 1;
      }

      if ((t % 4 == 0) || (t % 4 == 2))
      {
        p_align = plane_align0[p];
      }
      else
      {
        p_align = plane_align1[p];
      }

      size_t gw = tpc_z * weff[p] + tpc_y * offset;                          // global wire
      size_t gd = tpc_x * drift_size + apa_gap * (1 + tpc_x) / 2 + p_align;  // global drift
*/

  }
}



