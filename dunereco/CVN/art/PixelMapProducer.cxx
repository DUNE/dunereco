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
    fUnwrapped(2),
    fProtoDUNE(false)
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

//    bool hasLeftTPC = false;
//    bool hasRightTPC = false;

    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
    {

      geo::WireID wireid = cluster[iHit]->WireID();
      double temptdc                   = cluster[iHit]->PeakTime();
      unsigned int tempWire                  = wireid.Wire;
      unsigned int tempPlane              = wireid.Plane;
      // Leigh: Simple modification to unwrap the collection view wire plane
      if(!fProtoDUNE){
        if(fUnwrapped == 1){
          GetDUNEGlobalWireTDC(wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
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
      if(!fProtoDUNE){
        if(fUnwrapped == 1){
          GetDUNEGlobalWireTDC(wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
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
  void PixelMapProducer::GetDUNEGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc, 
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
  }

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
  
  }

}
