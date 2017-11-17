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
#include  "TVector2.h"
#include "lardataobj/RecoBase/Hit.h"

namespace cvn
{

  PixelMapProducer::PixelMapProducer(unsigned int nWire, unsigned int nTdc, double tRes):
  fNWire(nWire),
  fNTdc(nTdc),
  fTRes(tRes),
  fUnwrapped(false)
  {}


  PixelMap PixelMapProducer::CreateMap(std::vector< art::Ptr< recob::Hit > >& cluster)
  {

    Boundary bound = DefineBoundary(cluster);

    return CreateMapGivenBoundary(cluster, bound);



  }

  PixelMap PixelMapProducer::CreateMapGivenBoundary(std::vector< art::Ptr< recob::Hit > >& cluster,
                                                    const Boundary& bound)
  {

    PixelMap pm(fNWire, fNTdc, bound);

    for(size_t iHit = 0; iHit < cluster.size(); ++iHit)
      {
        geo::WireID wireid = cluster[iHit]->WireID();
        const double tdc                   = cluster[iHit]->PeakTime();
        unsigned int tempwire                  = wireid.Wire;
        unsigned int tempPlane              = wireid.Plane;
        // Leigh: Simple modification to unwrap the collection view wire plane
        if(fUnwrapped){
          GetGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,tempwire,tempPlane);
        }
        const double pe  = cluster[iHit]->Integral();
        const unsigned int wire = tempwire;
        const unsigned int wirePlane = tempPlane;
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
        unsigned int globalWire;
        unsigned int globalPlane;
        GetGlobalWire(wireid.Wire,wireid.Plane,wireid.TPC,globalWire,globalPlane);
        if(globalPlane==0){
          time_0.push_back(cluster[iHit]->PeakTime());
          wire_0.push_back(globalWire);
        }
        if(globalPlane==1){
          time_1.push_back(cluster[iHit]->PeakTime());
          wire_1.push_back(globalWire);
        }
        if(globalPlane==2){
          time_2.push_back(cluster[iHit]->PeakTime());
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


}
