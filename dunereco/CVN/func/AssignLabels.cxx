#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/TrainingData.h"
#include "larsim/MCCheater/BackTracker.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <iostream>
#include <iomanip>

namespace cvn
{
  /// Get Interaction_t from pdg, mode and iscc.
  /// Setting pdg and mode to zero triggers cosmic ray
  InteractionType GetInteractionType(simb::MCNeutrino& truth)
  {

    //if(truth->NeutrinoSet())
    //{
     int pdg = truth.Nu().PdgCode();
     bool iscc = truth.CCNC() == simb::kCC;
     int trueMode = truth.Mode();

     if(iscc)
     {
       if(abs(pdg) == 14)
       {
        switch(trueMode){
          case simb::kQE: return kNumuQE; break;
          case simb::kRes: return kNumuRes; break;
          case simb::kDIS: return kNumuDIS; break;
          default: return kNumuOther;
        }
      }
      else if(abs(pdg) == 12)
      {
        switch(trueMode){
          case simb::kQE: return kNueQE; break;
          case simb::kRes: return kNueRes; break;
          case simb::kDIS: return kNueDIS; break;
          default: return kNueOther;
        }
      }
      else if(abs(pdg) == 16)
      {
        switch(trueMode){
          case simb::kQE: return kNutauQE; break;
          case simb::kRes: return kNutauRes; break;
          case simb::kDIS: return kNutauDIS; break;
          default: return kNutauOther;
        }
      }
    }
    else if(trueMode==simb::kNuElectronElastic){
     return kNuElectronElastic;
   }
   else return kNC;
     //}
     //else return kCosmic;

 return kOther;
}

InteractionType GetInteractionTypeFromSlice(int pdg, bool iscc, int trueMode)
{

  if(iscc){

    if(abs(pdg) == 14){
      switch(trueMode){
        case simb::kQE: return kNumuQE; break;
        case simb::kRes: return kNumuRes; break;
        case simb::kDIS: return kNumuDIS; break;
        default: return kNumuOther;
      }
    }
    else if(abs(pdg) == 12){
      switch(trueMode){
        case simb::kQE: return kNueQE; break;
        case simb::kRes: return kNueRes; break;
        case simb::kDIS: return kNueDIS; break;
        default: return kNueOther;
      }
    }
    else if(abs(pdg) == 16){
      switch(trueMode){
        case simb::kQE: return kNutauQE; break;
        case simb::kRes: return kNutauRes; break;
        case simb::kDIS: return kNutauDIS; break;
        default: return kNutauOther;
      }
    }

  }
  else if(trueMode==simb::kNuElectronElastic){
    return kNuElectronElastic;
  }
  else return kNC;

  return kOther;
}



}
