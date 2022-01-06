////////////////////////////////////////////////////////////////////////
// \file    Assignlabels.h
///\brief   Utility class for truth labels
///
// \author  Leigh Whitehead leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////
#ifndef CVN_ASSIGNLABELS_H
#define CVN_ASSIGNLABELS_H

#include "art/Framework/Principal/Handle.h"

#include "dunereco/CVN/func/InteractionType.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"


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

    // Get the pion interaction mode for ProtoDUNE specific code
    unsigned short GetProtoDUNEBeamInteractionType(const simb::MCParticle &particle) const;
   
    private:

    // Recursive function to get all hits from daughters of a neutral particle
    unsigned int GetNeutralDaughterHitsRecursive(const simb::MCParticle &particle) const;

    int GetProcessKey(std::string process) const;
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
