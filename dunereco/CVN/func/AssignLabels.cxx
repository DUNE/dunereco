#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/TrainingData.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <iostream>
#include <iomanip>

namespace cvn
{
  /// Get Interaction_t from pdg, mode and iscc.
  /// Setting pdg and mode to zero triggers cosmic ray
  InteractionType AssignLabels::GetInteractionType(simb::MCNeutrino& truth)
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

  InteractionType AssignLabels::GetInteractionTypeFromSlice(int pdg, bool iscc, int trueMode)
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
  TopologyType AssignLabels::GetTopology(const simb::MCTruth &truth){

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
        top = static_cast<cvn::TopologyType>(top | kTopNutau);
      }   
    }
    else{
      top = static_cast<cvn::TopologyType>(top | kTopNC);
    }

    if(nu.Nu().PdgCode() < 0){
      top = static_cast<cvn::TopologyType>(top | kTopIsAntiNeutrino);
    }

    std::cout << "Topology after neutrino flavour = " << top << " :: " << nu.Nu().PdgCode() << std::endl;

    // Now we need to do some final state particle counting.
    unsigned int nParticle = truth.NParticles();

    unsigned short nProton = 0;
    unsigned short nPion = 0; // Charged pions, that is
    unsigned short nPizero = 0;
    unsigned short nNeutron = 0;

    // We need an instance of the backtracker to find the number of simulated hits for each track
    art::ServiceHandle<cheat::BackTrackerService> backTrack;

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

      // Find how many SimIDEs the track has
      unsigned int nSimIDE = backTrack->TrackIdToSimIDEs_Ps(truth.GetParticle(p).TrackId()).size();

      // Check if we have more than 100 MeV of kinetic energy
      float ke = truth.GetParticle(p).E() - truth.GetParticle(p).Mass();
      //    if( ke < 0.0){
      //      continue;
      //    }

      std::cout << "Final state particle " << pdg << " with ke " << ke << " and " << nSimIDE << " true hits" << std::endl;

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

  void AssignLabels::PrintTopology(TopologyType &top){
    unsigned short nu = GetPDGFromTopology(top);
    std::cout << "== Topology Information ==" << std::endl;
    std::cout << " - Raw Topology = " << top << std::endl;

    std::cout << " - Neutrino flavour (1 = NC) = " << nu << std::endl;

    unsigned short nProton = GetNProtons(top);
    std::cout << " - Number of protons (3 means >2) = " << nProton << std::endl;

    unsigned short nPion = GetNPions(top);
    std::cout << " - Number of charged pions (3 means >2) = " << nPion << std::endl;

    unsigned short nPizero = GetNPizeros(top);
    std::cout << " - Number of pizeros (3 means >2) = " << nPizero << std::endl;

    unsigned short nNeutron = GetNNeutrons(top);
    std::cout << " - Number of neutrons (3 means >2) = " << nNeutron << std::endl;

  }

  unsigned short AssignLabels::GetNProtons(TopologyType top){

    unsigned short nProton = 0;

    if(top & cvn::kTop0proton) nProton = 0;
    if(top & cvn::kTop1proton) nProton = 1;
    if(top & cvn::kTop2proton) nProton = 2;
    if(top & cvn::kTopNproton) nProton = 3;

    return nProton;
  }

  unsigned short AssignLabels::GetNPions(TopologyType top){

    unsigned short nPion = 0;

    if(top & cvn::kTop0pion) nPion = 0;
    if(top & cvn::kTop1pion) nPion = 1;
    if(top & cvn::kTop2pion) nPion = 2;
    if(top & cvn::kTopNpion) nPion = 3;

    return nPion;
  }

  unsigned short AssignLabels::GetNPizeros(TopologyType top){

    unsigned short nPizero = 0;

    if(top & cvn::kTop0pizero) nPizero = 0;
    if(top & cvn::kTop1pizero) nPizero = 1;
    if(top & cvn::kTop2pizero) nPizero = 2;
    if(top & cvn::kTopNpizero) nPizero = 3;

    return nPizero;
  }
  
  unsigned short AssignLabels::GetNNeutrons(TopologyType top){

    unsigned short nNeutron = 0;

    if(top & cvn::kTop0neutron) nNeutron = 0;
    if(top & cvn::kTop1neutron) nNeutron = 1;
    if(top & cvn::kTop2neutron) nNeutron = 2;
    if(top & cvn::kTopNneutron) nNeutron = 3;

    return nNeutron;
  }

  short AssignLabels::GetPDGFromTopology(TopologyType top){

    short nu = 0;

    if(top & cvn::kTopNue)   nu = 12;
    if(top & cvn::kTopNumu)  nu = 14;
    if(top & cvn::kTopNutau) nu = 16;
    if(top & cvn::kTopNC) nu = 1;

    if(top & cvn::kTopIsAntiNeutrino) nu = nu * -1;

    return nu;
  }

}


