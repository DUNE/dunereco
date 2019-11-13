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
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }

  PixelMapProducer::PixelMapProducer()
  {
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());  
    fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  }

  PixelMap PixelMapProducer::CreateMap(const std::vector< art::Ptr< recob::Hit > >& cluster)
  {
    std::vector<const recob::Hit*> newCluster;
    for(const art::Ptr<recob::Hit> hit : cluster){
      newCluster.push_back(hit.get());
    }
    return CreateMap(newCluster);
  }

  PixelMap PixelMapProducer::CreateMap(const std::vector<const recob::Hit* >& cluster)
  {
    Boundary bound = DefineBoundary(cluster);
    return CreateMapGivenBoundary(cluster, bound);
  }

  PixelMap PixelMapProducer::CreateMapGivenBoundary(const std::vector<const recob::Hit*>& cluster,
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
            GetDUNE10ktGlobalWireTDC(wireid.Wire,cluster[iHit]->PeakTime(),
              wireid.Plane,wireid.TPC,tempWire,tempPlane,temptdc);
          }
          // Default to 1x2x6. Should probably specifically name this function as such
          else {
            GetDUNEGlobalWireTDC(wireid.Wire,cluster[iHit]->PeakTime(),
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



  Boundary PixelMapProducer::DefineBoundary(const std::vector< const recob::Hit*>& cluster)
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
        if(fUnwrapped == 1) {
          if (fGeometry->DetectorName() == "dune10kt_v1") {
            if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue; // Skip dummy TPCs in 10kt module
            GetDUNE10ktGlobalWireTDC(wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
          }
          else {
            GetDUNEGlobalWireTDC(wireid.Wire,cluster[iHit]->PeakTime(),wireid.Plane,wireid.TPC,globalWire,globalPlane,globalTime);
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
    std::cout<<"minwire 0: "<<*minwireelement_0<<std::endl;
    auto minwireelement_1= std::min_element(wire_1.begin(), wire_1.end());
    std::cout<<"minwire 1: "<<*minwireelement_1<<std::endl;
    auto minwireelement_2= std::min_element(wire_2.begin(), wire_2.end());
    std::cout<<"minwire 2: "<<*minwireelement_2<<std::endl;

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
    if (tpcMod4 > 1 and globalPlane < 2) globalWire += fGeometry->Nwires(globalPlane, tpc, 0) + offset - localWire;
    else globalWire += localWire;

    if(tpcMod4 == 0 || tpcMod4 == 2){
      globalTDC = drift_size - localTDC;
    }
    else{
      globalTDC = localTDC + drift_size + apa_size;
    }
  }

  void PixelMapProducer::GetDUNE10ktGlobalWireTDC(unsigned int localWire, double localTDC, unsigned int plane, unsigned int tpc, 
                                                  unsigned int& globalWire, unsigned int& globalPlane, double& globalTDC) const
  {
    unsigned int nWiresTPC = 400;
    unsigned int wireGap = 4;
    double driftLen = fGeometry->TPC(tpc).DriftDistance();
    double apaLen = fGeometry->TPC(tpc).Width() - fGeometry->TPC(tpc).ActiveWidth();
    double driftVel = fDetProp->DriftVelocity();
    unsigned int drift_size = (driftLen / driftVel) * 2; // Time in ticks to cross a TPC 
    unsigned int apa_size   = 4*(apaLen / driftVel) * 2; // Width of the whole APA in TDC

    // std::cout << "Drift size is " << drift_size << ", APA size is " << apa_size << std::endl;

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
    if (tpc_xy > 3 and globalPlane < 2) globalWire += fGeometry->Nwires(globalPlane, tpc, 0) + offset - localWire;
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
    GetProtoDUNEGlobalWire(localWire, plane, tpc, globalWire, globalPlane);

    // Implement additional mirroring here?

  } // function GetProtoDUNEGlobalWireTDC

  SparsePixelMap PixelMapProducer::CreateSparseMap(std::vector< art::Ptr< recob::Hit> >& cluster,
    bool usePixelTruth) {

    SparsePixelMap map(3, 3, usePixelTruth);

    art::ServiceHandle<cheat::BackTrackerService> bt;
    art::ServiceHandle<cheat::ParticleInventoryService> pi;
     
    //int count =0;
    
    for(size_t iHit = 0; iHit < cluster.size(); ++iHit) {

      geo::WireID wireid       = cluster[iHit]->WireID();
      double globalTime        = cluster[iHit]->PeakTime();
      unsigned int globalWire  = wireid.Wire;
      unsigned int globalPlane = wireid.Plane;

      if (fGeometry->DetectorName().find("1x2x6") != std::string::npos) {
        GetDUNEGlobalWireTDC(wireid.Wire, cluster[iHit]->PeakTime(),
          wireid.Plane, wireid.TPC, globalWire, globalPlane, globalTime);
      }
      else if (fGeometry->DetectorName() == "dune10kt_v1") {
        if (wireid.TPC%6 == 0 or wireid.TPC%6 == 5) continue;
        GetDUNE10ktGlobalWireTDC(wireid.Wire, cluster[iHit]->PeakTime(),
          wireid.Plane, wireid.TPC, globalWire, globalPlane, globalTime);
      }
      else if (fGeometry->DetectorName() == "protodune") {
        GetProtoDUNEGlobalWireTDC(wireid.Wire, cluster[iHit]->PeakTime(),
          wireid.Plane, wireid.TPC, globalWire, globalTime, globalPlane);
      }
      else throw art::Exception(art::errors::UnimplementedFeature)
        << "Geometry " << fGeometry->DetectorName() << " not implemented "
        << "in CreateSparseMap." << std::endl;


      if (usePixelTruth) {
        // Get true particle and PDG responsible for this hit
        std::vector<sim::TrackIDE> IDEs = bt->HitToTrackIDEs(cluster[iHit]);
           // Get track ids, PDG and energy responsibles for this hit
        std::vector<int> tracks, pdgs, pdgParent;
        std::vector<float> energy;
        std::vector<std::string> process;
        for (auto k : IDEs){
            tracks.push_back(k.trackID);
            simb::MCParticle p = pi->TrackIdToParticle(k.trackID);
            process.push_back(p.Process());
            if ( p.Mother()!=0){
                  int pdg_parent = pi->TrackIdToParticle(p.Mother()).PdgCode();
                  //std::cout << " track ID  " << k.trackID << " PDG "<< p.PdgCode() << " Parent "<< pdg_parent <<"process "  << p.Process() << std::endl;
                  pdgParent.push_back(pdg_parent);

                 // std::cout<< "pdg delta parent: "<< pdgDelta_Parent << std::endl;
            } 
          
            // Manually check to see if we have a Michel electron here
            int pdg = p.PdgCode();
            if (pdg == 11 && p.Mother() != 0) {
               int pdgParent = pi->TrackIdToParticle(p.Mother()).PdgCode();
             //  std::cout <<" Michel e-  parent "<< pdgParent << std::endl; 
               if (pdgParent == 13 && p.Process() == "muMinusCaptureAtRest") {
                 // If it's a Michel electron, tag it with an unphysical PDG code
                 pdg = 99;
               }
             } // If non-primary electron
             pdgs.push_back(pdg);
             energy.push_back(k.energy);
           } 
           for (auto k : process)
          {
              std:: cout << "size " <<  process.size() << " process: " <<  k << std::endl;
          }
           map.AddHit(globalPlane, {globalWire, (unsigned int)globalTime, wireid.TPC},
             cluster[iHit]->Integral(), pdgs, tracks, energy, pdgParent, process) ; 
      } // if PixelTuth 
      

      else {
        map.AddHit(globalPlane, {globalWire, (unsigned int)globalTime},
          cluster[iHit]->Integral());
      }
    } // for iHit

    return map;

  } // function PixelMapProducer::CreateSparseMap

} // namespace cvn
