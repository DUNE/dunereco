////////////////////////////////////////////////////////////////////////
/// \file    RegPixelMap3DProducer.h
/// \brief   RegPixelMap3DProducer for RegCNN modified from PixelMapProducer.h
/// \author  Wenjie Wu - wenjieww@uci.edu
////////////////////////////////////////////////////////////////////////

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>

#include "larcorealg/Geometry/Exceptions.h" // geo::InvalidWireError
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dune/RegCNN/art/RegPixelMap3DProducer.h"
#include  "TVector2.h"

namespace cnn
{
  RegPixelMap3DProducer::RegPixelMap3DProducer(unsigned int nbinsX, float xmin, float xmax, 
                        unsigned int nbinsY, float ymin, float ymax,
                        unsigned int nbinsZ, float zmin, float zmax):
    fNBinsX(nbinsX),
    fNBinsY(nbinsY),
    fNBinsZ(nbinsZ),
    fLengthX(xmax-xmin),
    fLengthY(ymax-ymin),
    fLengthZ(zmax-zmin)
    {
    }


  RegPixelMap3D RegPixelMap3DProducer::Create3DMap(detinfo::DetectorClocksData const& clockData,
                                             detinfo::DetectorPropertiesData const& detProp,
                                             std::vector< art::Ptr< recob::Hit > > const& cluster,
                                             art::FindManyP<recob::SpacePoint> const& fmSPFromHits,
                                             const std::vector<float> &vtx)
  {

      RegCNNBoundary3D bound = Define3DBoundary(detProp, cluster, vtx); 

      return Create3DMapGivenBoundaryBySP(clockData, detProp, cluster, bound, fmSPFromHits);
  }


  RegPixelMap3D RegPixelMap3DProducer::Create3DMapGivenBoundaryBySP(detinfo::DetectorClocksData const& clockData,
                                                          detinfo::DetectorPropertiesData const& detProp,
                                                          std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                          const RegCNNBoundary3D& bound,
                                                          art::FindManyP<recob::SpacePoint> const& fmSPFromHits)
  {
      std::cout<<"create 3D pixel maps"<<std::endl;

      RegPixelMap3D pm(bound, false);

      if (!fmSPFromHits.isValid()) return pm;

      // Loop over hit
      unsigned int nhits = cluster.size();
      for (unsigned int ihit= 0; ihit< nhits; ++ihit) {
          // Get hit
          art::Ptr<recob::Hit> hit = cluster.at(ihit);
          std::vector<art::Ptr<recob::SpacePoint> > sp = fmSPFromHits.at(ihit);

          std::vector<float> coordinates(3);
          if (!sp.empty()) {
              coordinates[0] = sp[0]->XYZ()[0];
              coordinates[1] = sp[0]->XYZ()[1];
              coordinates[2] = sp[0]->XYZ()[2];
              if (coordinates[0] == -999) continue;
          }

          std::vector<float> shifted_coordinates(3);
          shifted_coordinates[0] = coordinates[0] - bound.Center(0) + bound.Length(0)/2;
          shifted_coordinates[1] = coordinates[1] - bound.Center(1) + bound.Length(1)/2;
          shifted_coordinates[2] = coordinates[2] - bound.Center(2) + 1000./8;
          //std::cout<<"HELLO ===> "<<shifted_coordinates[0];
          //std::cout<<", "<<shifted_coordinates[1];
          //std::cout<<", "<<shifted_coordinates[2]<<std::endl;

          float correctedHitCharge = (hit->Integral()*TMath::Exp((sampling_rate(clockData)*hit->PeakTime()) / (detProp.ElectronLifetime()*1.e3) ) );

          pm.AddHit(shifted_coordinates[0], shifted_coordinates[1], shifted_coordinates[2], correctedHitCharge);
      }

      return pm;
  }


  RegCNNBoundary3D RegPixelMap3DProducer::Define3DBoundary(detinfo::DetectorPropertiesData const& detProp,
                                                     std::vector< art::Ptr< recob::Hit > > const& cluster,
                                                     const std::vector<float> &vtx)
  {
      RegCNNBoundary3D bound(fNBinsX, fNBinsY, fNBinsZ, fLengthX, fLengthY, fLengthZ, vtx[0], vtx[1], vtx[2]);

      return bound;
  }

  std::ostream& operator<<(std::ostream& os, const RegPixelMap3DProducer& p)
  {
    os << "RegPixelMapProducer...... ";
    return os;
  }

}
