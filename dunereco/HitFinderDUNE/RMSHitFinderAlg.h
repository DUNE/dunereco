#ifndef RMSHITFINDERALG_H
#define RMSHITFINDERALG_H

#include "fhiclcpp/ParameterSet.h"

#include "TMath.h"

#include "RobustHitFinderSupport.h"

#include <memory>
#include <algorithm>
#include <map>

namespace dune {
  class RMSHitFinderAlg {
public:
    RMSHitFinderAlg(fhicl::ParameterSet const& p);

    void reconfigure(fhicl::ParameterSet const& p);

    void FindHits(dune::ChannelInformation & chan);

private:
    void FilterWaveform(std::vector<float> wf, std::vector<float> & fwf);
    void RobustRMSBase(std::vector<float> wf, float & bl, float & r);
    void FindPulses(std::vector<float> wf, float bl, float r, std::vector<std::pair<int,int> > & pulse_ends);
    void MergeHits(std::vector<std::pair<int,int> > & pulse_ends);

    int fWindowWidth;
    float fFilterWidth;
    float fSigmaRiseThreshold;
    float fSigmaFallThreshold;

  };

}

#endif
