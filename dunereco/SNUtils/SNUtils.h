#ifndef SNUTILS_H
#define SNUTILS_H

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TMath.h"

#include <vector>

#include "lardataobj/RecoBase/OpHit.h"

#include "dune/SNSlicer/SNSlice.h"

namespace{
  double mysqr(double x){return x*x;}
}

namespace sn
{
  //---------------------------------------------------------------------------
  double TimeCluster(std::vector<std::pair<double, double>>& pts,
                     double& clustPE)
  {
    // "Clustering by Fast Search and Find of Density Peaks"
    // Rodriguez and Laio, Science, 2014

    const double sigma = 5; // characteristic time
    const double sigsq = mysqr(sigma);

    //    const double norm = 1/(sqrt(2*TMath::Pi())*sigma);

    std::sort(pts.begin(), pts.end()); // Sort points by time

    std::vector<double> rhos(pts.size(), 0); // densities

    for(unsigned int i = 0; i < pts.size(); ++i){
      double& rho = rhos[i];

      const double ti = pts[i].first;

      for(unsigned int j = 0; j < pts.size(); ++j){
        //        if(j == i) continue;

        const double tj = pts[j].first;
        const double qj = pts[j].second;

        rho += qj*exp(-mysqr(ti-tj)/(2*sigsq));
      } // end for j
    } // end for i

    double bestScore = 0;
    double bestT = 0;
    for(unsigned int i = 0; i < pts.size(); ++i){
      if(rhos[i] > bestScore){
        bestScore = rhos[i];
        bestT = pts[i].first;
        clustPE = rhos[i];
      }
    }
    return bestT;
  }

  //---------------------------------------------------------------------------
  // Pass a vector of OpHits and fields from the SNSlice. Returns best estimate
  // of the slice time, and the total number of photons in the ophits it based
  // that on (as flashPhots).
  double FindSliceTZero(const std::vector<recob::OpHit>& ophits,
                        double sliceMeanT,
                        double sliceMeanZ,
                        double& flashPhots)
  {
    const geo::GeometryCore& geom(*lar::providerFrom<geo::Geometry>());
    const detinfo::DetectorProperties& dp(*lar::providerFrom<detinfo::DetectorPropertiesService>());

    const double driftLen = geom.Cryostat(0).TPC(0).DriftDistance();
    const double driftT = driftLen / dp.DriftVelocity();

    double ophitTotQ = 0;

    // Times and PE weights, to calculate median
    std::vector<std::pair<double, double>> ophitTs;

    for(const recob::OpHit& ophit: ophits){
      double xyz[3];
      geom.OpDetGeoFromOpChannel(ophit.OpChannel()).GetCenter(xyz);

      // Only use hits that are close enough to be consistent
      if(ophit.PeakTime() < sliceMeanT &&
         ophit.PeakTime() > sliceMeanT-driftT &&
         fabs(xyz[2]-sliceMeanZ) < 600){ // 2sigma?*/

        ophitTotQ += ophit.PE();
        ophitTs.emplace_back(ophit.PeakTime(), ophit.PE());
      }
    } // end for ophit

    return TimeCluster(ophitTs, flashPhots);
  }

  //---------------------------------------------------------------------------
  // Corrects for PE->MeV as well as lifetime of the drifting charge
  double AttenFactor(double dt)
  {
    // From an expo fit of totQ/trueE vs dt
    //    return exp(4.83617-0.00037413*dt);
    return 128.734*exp(-0.000357316*dt);
  }

  //---------------------------------------------------------------------------
  double SliceEnergy(const sn::SNSlice& slice,
                     const std::vector<recob::OpHit>& ophits)
  {
    double junk;
    const double t0 = FindSliceTZero(ophits, slice.meanT, slice.meanZ, junk);
    const double dt = slice.meanT - t0;

    return slice.totQ / AttenFactor(dt);
  }
}

#endif
