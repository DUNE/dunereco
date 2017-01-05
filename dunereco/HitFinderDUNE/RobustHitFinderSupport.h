#ifndef ROBUSTHITFINDERSUPPORT_H
#define ROBUSTHITFINDERSUPPORT_H

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "canvas/Persistency/Common/Ptr.h"

#include <memory>
#include <algorithm>
#include <vector>
#include <map>
#include <utility>

namespace dune {

  class HitInformation {
public:
    float pitchgamma;

    float hitIntegral;
    float hitIntegralFilter;
    float hitSigmaIntegral;
    float hitSigmaIntegralFilter;
    float hitAmplitude;
    float hitAmplitudeFilter;
    float hitPeakTick;
    float hitPeakTickFilter;
    float hitPeakTime;
    float hitPeakTimeFilter;
    int hitBeginTick;
    int hitEndTick;
    int hitWidth;     //in units of ticks
    float hitSumADC;

    float hitx;
    float hity;
    float hitz;
    float hiterrxlo;
    float hiterrxhi;
    float hiterrylo;
    float hiterryhi;
    float hiterrzlo;
    float hiterrzhi;

    float perpdist;
    float hitt;
    float driftdist;

    bool fitrealhit;
    bool countercut;
    bool assumedhit;

    int channelID;

    recob::Hit artHit;
  };

  typedef std::vector<HitInformation> HitVec_t;

  class ChannelInformation {
public:
    int signalSize;
    std::vector<float> signalVec;
    std::vector<float> signalFilterVec;
    art::Ptr<recob::Wire> artWire;
    art::Ptr<raw::RawDigit> artRawDigit;
    float baseline;
    float rms;
    float baselineFilter;
    float rmsFilter;
    int channelID;
    int wireID;
    int tpcNum;
    float chanz;
    int nGoodHits;
    float goodHitStartTick;
    float goodHitEndTick;
    std::vector<std::pair<int,int> > pulse_ends;
  };

  typedef std::map<int,ChannelInformation> ChanMap_t;

}

#endif
