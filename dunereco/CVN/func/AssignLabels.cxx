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

// This function uses purely the information from the neutrino generator to
// find all of the final-state particles that contribute to the event.
TopologyType GetTopology(const simb::MCTruth &truth){

  TopologyType top = cvn::kTopUnset;

  const simb::MCNeutrino &nu = truth.GetNeutrino();

  // Find out the neutrino flavour and CC/NC
  if(nu.CCNC() == simb::kCC){
    if(abs(nu.Nu().PdgCode()) == 12){
      top = static_cast<cvn::TopologyType>(top | kTopNue);
    }
    else if(abs(nu.Nu().PdgCode()) == 14){
      top = static_cast<cvn::TopologyType>(top | kTopNumu);
    }
    else if(abs(nu.Nu().PdgCode()) == 16){
      top = static_cast<cvn::TopologyType>(top | kTopNue);
    }   
  }
  else{
    top = static_cast<cvn::TopologyType>(top | kTopNC);
  }

  std::cout << "Topology after neutrino flavour = " << top << " :: " << nu.Nu().PdgCode() << std::endl;

  // Now we need to do some final state particle counting.
  unsigned int nParticle = truth.NParticles();

  unsigned short nProton = 0;
  unsigned short nPion = 0; // Charged pions, that is
  unsigned short nPizero = 0;
  unsigned short nNeutron = 0;
  
  // Loop over all of the particles
  for(unsigned int p = 0; p < nParticle; ++p){

    int pdg = truth.GetParticle(p).PdgCode();

    // Make sure this is a final state particle
    if(truth.GetParticle(p).StatusCode() != 1){
      continue;
    }

    // GENIE has some fake particles for energy conservation - eg nuclear binding energy. Ignore these
    if(pdg > 2000000000){
      continue;
    }

    // Check if we have more than 100 MeV of kinetic energy
    float ke = truth.GetParticle(p).E() - truth.GetParticle(p).Mass();
//    if( ke < 0.0){
//      continue;
//    }

    std::cout << "Final state particle " << pdg << " with ke " << ke << std::endl;

    switch(abs(pdg)){
      case 111 : ++nPizero;  break;
      case 211 : ++nPion;    break;
      case 2112: ++nNeutron; break;
      case 2212: ++nProton;  break;
      default  : break;
    }

  }

  // Assign the enums based on the counters
  switch(nProton){
    case 0  : top = static_cast<cvn::TopologyType>(top | kTop0proton); break;
    case 1  : top = static_cast<cvn::TopologyType>(top | kTop1proton); break;
    case 2  : top = static_cast<cvn::TopologyType>(top | kTop2proton); break;
    default : top = static_cast<cvn::TopologyType>(top | kTopNproton); break;
  }

  switch(nPion){
    case 0  : top = static_cast<cvn::TopologyType>(top | kTop0pion); break;
    case 1  : top = static_cast<cvn::TopologyType>(top | kTop1pion); break;
    case 2  : top = static_cast<cvn::TopologyType>(top | kTop2pion); break;
    default : top = static_cast<cvn::TopologyType>(top | kTopNpion); break;
  }

  switch(nPizero){
    case 0  : top = static_cast<cvn::TopologyType>(top | kTop0pizero); break;
    case 1  : top = static_cast<cvn::TopologyType>(top | kTop1pizero); break;
    case 2  : top = static_cast<cvn::TopologyType>(top | kTop2pizero); break;
    default : top = static_cast<cvn::TopologyType>(top | kTopNpizero); break;
  }

  switch(nNeutron){
    case 0  : top = static_cast<cvn::TopologyType>(top | kTop0neutron); break;
    case 1  : top = static_cast<cvn::TopologyType>(top | kTop1neutron); break;
    case 2  : top = static_cast<cvn::TopologyType>(top | kTop2neutron); break;
    default : top = static_cast<cvn::TopologyType>(top | kTopNneutron); break;
  }

  std::cout << "Final topology : " << top << " with particle counts " << nProton << ", " << nPion << ", " << nPizero << ", " << nNeutron << std::endl;

  return top;

}


}


