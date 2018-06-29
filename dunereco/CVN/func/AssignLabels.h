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

#include "nusimdata/SimulationBase/MCTruth.h"


namespace cvn
{

  class AssignLabels{

    public:

    AssignLabels(){};
    ~AssignLabels(){};

    InteractionType GetInteractionType(simb::MCNeutrino& truth);
    InteractionType GetInteractionTypeFromSlice(int nuPDG, bool nuCCNC,
                                                int nuMode );

    // Use the topology information
    TopologyType GetTopology(const simb::MCTruth &truth);
    void PrintTopology(TopologyType &top);
    unsigned short GetNProtons(TopologyType top);
    unsigned short GetNPions(TopologyType top);
    unsigned short GetNPizeros(TopologyType top);
    unsigned short GetNNeutrons(TopologyType top);
    short GetPDGFromTopology(TopologyType top);
    
  };
}

#endif // CVN_ASSIGNLABELS_H
