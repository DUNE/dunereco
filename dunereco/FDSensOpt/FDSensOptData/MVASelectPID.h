#ifndef FDSENSOPT_MVASELECTPID_H
#define FDSENSOPT_MVASELECTPID_H

namespace dunemva
{
  class MVASelectPID
  {
  public:

    int selectMode; ///< What neutrino species are we IDing? Use PDG code to identify numu or nue selection.
    double pid; ///< How confident are we?

    // Input variables
    float evtcharge;
    float rawcharge;
    float wirecharge;

    float ntrack;
    float avgtrklength;
    float maxtrklength;
    float trkdedx;
    float trkrch;
    float trkrt;
    float trkfr;
    float trkpida_save;
    float nshower;
    float showerdedx;
    float eshower;
    float frshower;
    float nhitspershw;
    float shwlength;
    float shwmax;
    float fract_5_wires;
    float fract_10_wires;
    float fract_50_wires;
    float fract_100_wires;
    float shwdis;
    float shwdisx;
    float shwdisy;
    float shwdisz;
    float shwcosx;
    float shwcosy;
    float shwcosz;
    float trkcosx;
    float trkcosy;
    float trkcosz;
    float et;

  };
}

#endif
