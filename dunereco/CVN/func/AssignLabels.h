////////////////////////////////////////////////////////////////////////
// \file    Assignlabels.h
///\brief   Defines an enumeration for prong classification
///
// \author psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef CVN_ASSIGNLABELS_H
#define CVN_ASSIGNLABELS_H


#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/PixelMap.h"

#include "nusimdata/SimulationBase/MCTruth.h"


namespace cvn
{

  InteractionType GetInteractionType(simb::MCNeutrino& truth);
  InteractionType GetInteractionTypeFromSlice(int nuPDG, bool nuCCNC,
                                              int nuMode );
}

#endif // CVN_ASSIGNLABELS_H
