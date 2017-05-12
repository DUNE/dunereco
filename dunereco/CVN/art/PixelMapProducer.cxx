////////////////////////////////////////////////////////////////////////
/// \file    PixelMapProducer.h
/// \brief   PixelMapProducer for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
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
  fTRes(tRes)
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
        const unsigned int wire                  = wireid.Wire;
        const unsigned int wirePlane              = wireid.Plane;

        const double pe  = cluster[iHit]->Integral();
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
        if(wireid.Plane==0){
          time_0.push_back(cluster[iHit]->PeakTime());
          wire_0.push_back(wireid.Wire);
        }
        if(wireid.Plane==1){
          time_1.push_back(cluster[iHit]->PeakTime());
          wire_1.push_back(wireid.Wire);
        }
        if(wireid.Plane==2){
          time_2.push_back(cluster[iHit]->PeakTime());
          wire_2.push_back(wireid.Wire);
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




}
