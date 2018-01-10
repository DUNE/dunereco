#ifndef FDSENSOPT_ENERGYRECOOUTPUT_H
#define FDSENSOPT_ENERGYRECOOUTPUT_H

namespace dune
{
  class EnergyRecoOutput
  {
  public:
    int recoMethodUsed;         // 1 = longest reco track + hadronic, 2 = reco shower with highest charge + hadronic, 3 = all hit charges, -1 = not set 
    double nuEnergy;
    double lepEnergy;
    double hadEnergy;
    int longestTrackContained;  // 1 = contained, 0 = exiting, -1 = not set
    int trackMomMethod;         // 1 = range, 0 = MCS, -1 = not set
  };
}

#endif
