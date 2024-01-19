#ifndef FDSENSOPT_ANGULARRECOOUTPUT_H
#define FDSENSOPT_ANGULARRECOOUTPUT_H

#include "Math/GenVector/PositionVector3D.h" 
#include "Math/GenVector/PxPyPzE4D.h" 
#include "Math/GenVector/LorentzVector.h" 

using Point_t = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>;
using Direction_t = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>;
using Momentum_t = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>>;

namespace dune
{
  class AngularRecoOutput
  {
  public:

    int recoMethodUsed;         // 1 = longest reco direction, 2 = shower with highest charge direction, -1 = not set 

    Point_t fRecoVertex;
    Direction_t fRecoDirection;
  };
}

#endif
