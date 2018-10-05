////////////////////////////////////////////////////////////////////////
// \file    Assignlabels.h
///\brief   Utility class for ruth labels
///
// \author  Leigh Whitehead leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////
#ifndef CVN_ASSIGNLABELS_H
#define CVN_ASSIGNLABELS_H


#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/PixelMap.h"
#include "art/Framework/Principal/Handle.h"
#include "nusimdata/SimulationBase/MCTruth.h"


namespace cvn
{

  class AssignLabels{

    public:

    AssignLabels();
    ~AssignLabels(){};

    InteractionType GetInteractionType(simb::MCNeutrino& truth);
    InteractionType GetInteractionTypeFromSlice(int nuPDG, bool nuCCNC,
                                                int nuMode );

    // Use the topology information
    void GetTopology(const art::Ptr<simb::MCTruth> truth, unsigned int nTopologyHits);
    void PrintTopology();
    unsigned short GetNProtons()  { return nProton;  };
    unsigned short GetNPions()    { return nPion;    };
    unsigned short GetNPizeros()  { return nPizero;  };
    unsigned short GetNNeutrons() { return nNeutron; };
    short GetPDG()       { return pdgCode;  };
    unsigned short TauMode()      { return tauMode;  };
    bool IsAntineutrino();
    unsigned short GetTopologyType();
    unsigned short GetTopologyTypeAlt();
    
    private:

    // unused unsigned int fTopologyHitsCut;

    unsigned short nProton;
    unsigned short nPion;
    unsigned short nPizero;
    unsigned short nNeutron;
    short pdgCode;
    unsigned short tauMode;
  
  };
}

#endif // CVN_ASSIGNLABELS_H
