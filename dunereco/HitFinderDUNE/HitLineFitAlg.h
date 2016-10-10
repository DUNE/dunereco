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
    
    struct HitLineFitData {
      double hitHoriz;
      double hitHorizErrLo;
      double hitHorizErrHi;
      double hitVert;
      double hitVertErrLo;
      double hitVertErrHi;
      bool hitREAL;
    };

    struct HitLineFitResults {
      std::map<int,float> bestVal;
      std::map<int,float> bestValError;
      double chi2;
      double sum2resid;
      double mle;
      int ndf;
      bool fitsuccess;
    };

    struct ParVals {
      double start;
      double min;
      double max;
    };

    HitLineFitAlg(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& p);

    int FitLine(std::vector<HitLineFitData> & data, HitLineFitResults & bestfit);
    void SetParameter(int i, double startValue, double minValue, double maxValue);
    void SetHorizVertRanges(float hmin, float hmax, float vmin, float vmax);
    void SetSeed(UInt_t seed) {
      fSeedValue = seed;
    }

private:
    float PointToLineDist(TVector3 ptloc, TVector3 linept1, TVector3 linept2);
    void DeterministicShuffle(std::vector<unsigned int> & vec);
    bool CheckModelParameters();

    float fVertRangeMin;
    float fVertRangeMax;
    float fHorizRangeMin;
    float fHorizRangeMax;

    std::map<int,ParVals> fParIVal; 

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
