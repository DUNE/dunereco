#ifndef HITLINEFITALG_H
#define HITLINEFITALG_H

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RobustHitFinderSupport.h"

#include "TVector3.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"

#include <algorithm>
#include <memory>
#include <vector>

namespace dune {

  class HitLineFitAlg {
public:

    struct HitLineFitResults {
      float fitconstant;
      float fitconstanterr;
      float fitlinear;
      float fitlinearerr;
      float fitquadratic;
      float fitquadraticerr;
      float fitchi2;
      float fitsumsqrresidual;
      float fitmle;
      float fitndf;
      bool fitsuccess;
    };

    HitLineFitAlg(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& p);

    int FindTrack(std::vector<dune::HitInformation> & data, HitLineFitResults & bestfit);
    void SetCounterPositions(float c1vert, float c1horiz, float c2vert, float c2horiz);
    void SetRanges(float hmin, float hmax, float vmin, float vmax);
    void SetSeed(UInt_t seed) {
      fSeedValue = seed;
    }

private:
    float hitGeomDist(TVector3 hitloc, TVector3 trigloc1, TVector3 trigloc2);
    void DeterministicShuffle(std::vector<unsigned int> & vec);

    float fC1Vert;
    float fC1Horiz;
    float fC2Vert;
    float fC2Horiz;
    float fVertRangeMin;
    float fVertRangeMax;
    float fHorizRangeMin;
    float fHorizRangeMax;

    UInt_t fSeedValue;

    int fFitPolN;
    int fMinStartPoints;
    int fMinAlsoPoints;
    float fIterationsMultiplier;
    float fInclusionThreshold;
    int fLogLevel;
  };

}

#endif
