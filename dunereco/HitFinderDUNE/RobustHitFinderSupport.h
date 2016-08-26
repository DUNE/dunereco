#ifndef ROBUSTHITFINDERSUPPORT_H
#define ROBUSTHITFINDERSUPPORT_H

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "canvas/Persistency/Common/Ptr.h"

#include <memory>
#include <algorithm>
#include <vector>
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

    float hitx;
    float hity;
    float hitz;
    float hiterrxlo;
    float hiterrxhi;
    float hiterrylo;
    float hiterryhi;
    float hiterrzlo;
    float hiterrzhi;
    float hithoriz;
    float hitvert;
    float hithorizerrlo;
    float hithorizerrhi;
    float hitverterrlo;
    float hitverterrhi;
    float perpdist;
    float hitt;
    float driftdist;

    bool fitrealhit;
    bool countercut;

    int channelID;

    recob::Hit artHit;
  };

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
    int tpcNum;
    std::vector<std::pair<int,int> > pulse_ends;
  };

}

#endif
