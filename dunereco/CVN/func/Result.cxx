////////////////////////////////////////////////////////////////////////
/// \file    Result.h
/// \brief   Result for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead - leigh.howard.whitehead@cern.ch
///          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <ostream>
#include <algorithm>

#include "dune/CVN/func/Result.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace cvn
{

  Result::Result(const float* output, unsigned int& nOutputs):
  fOutput(1)
  {
    fOutput[0].resize(nOutputs);
    for(size_t i = 0; i < nOutputs; ++i)
    {
     fOutput[0][i] = output[i];
    }
  }

  Result::Result(const std::vector< std::vector<float> > output){
    fOutput = output; 
  }

  Result::Result():
  fOutput()
  {}

  unsigned int Result::ArgMax(int output_n) const
  {
    // Get the max element iterator and convert to vector index

    // single-output
    if(fOutput.size() == 1)
        return std::distance(fOutput[0].begin(),std::max_element(fOutput[0].begin(),fOutput[0].end()));
    // multi-output
    return std::distance(fOutput[output_n].begin(),std::max_element(fOutput[output_n].begin(),fOutput[output_n].end()));
  }

  /*
  float Result::Max(){
    // Get the maximum value by dereferencing the iterator
    return *std::max_element(fOutput.begin(),fOutput.end());
  }

  unsigned int Result::NOutput(){
    return fOutput.size();
  }
  */

  /// Return the predicted is_antineutrion
  TFIsAntineutrino Result::PredictedIsAntineutrino() const{
    return static_cast<TFIsAntineutrino>((int)round(this->GetIsAntineutrinoProbability()));
  }

  /// Return the predicted flavour
  TFFlavour Result::PredictedFlavour() const
  {
    return static_cast<TFFlavour>(this->ArgMax(TFMultioutputs::flavour));
  }

  /// Return the predicted interaction
  TFInteraction Result::PredictedInteraction() const{
    return static_cast<TFInteraction>(this->ArgMax(TFMultioutputs::interaction));
  }

  /// Return the predicted protons
  TFTopologyProtons Result::PredictedProtons() const{
    return static_cast<TFTopologyProtons>(this->ArgMax(TFMultioutputs::protons));
  }

  /// Return the predicted pions
  TFTopologyPions Result::PredictedPions() const{
    return static_cast<TFTopologyPions>(this->ArgMax(TFMultioutputs::pions));
  }

  /// Return the predicted pizeros
  TFTopologyPizeros Result::PredictedPizeros() const{
    return static_cast<TFTopologyPizeros>(this->ArgMax(TFMultioutputs::pizeros));
  }

  /// Return the predicted neutrons
  TFTopologyNeutrons Result::PredictedNeutrons() const{
    return static_cast<TFTopologyNeutrons>(this->ArgMax(TFMultioutputs::neutrons));
  }

  /// Return the is_antineutrino probability
  float Result::GetIsAntineutrinoProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no is_antineutrino output
    // multi-output
    return fOutput[TFMultioutputs::is_antineutrino][0]; 
  }

  /// Return the numu flavour probability
  float Result::GetNumuProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return fOutput[0][TFResultType::kTFNumuQE] + fOutput[0][TFResultType::kTFNumuRes]
           + fOutput[0][TFResultType::kTFNumuDIS] + fOutput[0][TFResultType::kTFNumuOther];
    // multi-output
    return fOutput[TFMultioutputs::flavour][TFFlavour::kFlavNumuCC];
  }
  
  /// Return the nue flavour probability
  float Result::GetNueProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return fOutput[0][TFResultType::kTFNueQE] + fOutput[0][TFResultType::kTFNueRes]
           + fOutput[0][TFResultType::kTFNueDIS] + fOutput[0][TFResultType::kTFNueOther];
    // multi-output
    return fOutput[TFMultioutputs::flavour][TFFlavour::kFlavNueCC];
  }

  /// Return the nutau flavour probability
  float Result::GetNutauProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return fOutput[0][TFResultType::kTFNutauQE] + fOutput[0][TFResultType::kTFNutauRes]
           + fOutput[0][TFResultType::kTFNutauDIS] + fOutput[0][TFResultType::kTFNutauOther];
    // multi-output
    return fOutput[TFMultioutputs::flavour][TFFlavour::kFlavNutauCC];
  }

  /// Return the NC probability
  float Result::GetNCProbability() const
  {

    // The old caffe network didn't give us an NC probability
    // So make sure we have enough values to grab it
    float result = -999;

    // single-output
    if(fOutput.size() == 1){
      if(fOutput[0].size() > static_cast<unsigned int>(TFResultType::kTFNC)){
        result = fOutput[0][TFResultType::kTFNC];
      }
      else{
        mf::LogError("cvn::Result") << "Output vector too short to include an NC probability" << std::endl;
      }
      return result;
    }
    // multi-output
    return fOutput[TFMultioutputs::flavour][TFFlavour::kFlavNC];
  }

  /// Return the CC QE interaction probability
  float Result::GetQEProbability() const
  { 
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no interaction probability
    // multi-output
    return fOutput[TFMultioutputs::interaction][TFInteraction::kInteQECC];
  }

  /// Return the CC Res interaction probability
  float Result::GetResProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no interaction probability
    // multi-output
    return fOutput[TFMultioutputs::interaction][TFInteraction::kInteResCC];
  }

  /// Return the CC DIS interaction probability
  float Result::GetDISProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no interaction probability
    // multi-output
    return fOutput[TFMultioutputs::interaction][TFInteraction::kInteDISCC];
  }
  
  /// Return the CC Other interaction probability
  float Result::GetOtherProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no interaction probability
    // multi-output
    return fOutput[TFMultioutputs::interaction][TFInteraction::kInteOtherCC];
  }

  /// Return the 0 protons topology probability
  float Result::Get0protonsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no protons probability
    // multi-output
    return fOutput[TFMultioutputs::protons][TFTopologyProtons::kTop0proton]; 
  }

  /// Return the 1 protons topology probability
  float Result::Get1protonsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no protons probability
    // multi-output
    return fOutput[TFMultioutputs::protons][TFTopologyProtons::kTop1proton]; 
  }

  /// Return the 2 protons topology probability
  float Result::Get2protonsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no protons probability
    // multi-output
    return fOutput[TFMultioutputs::protons][TFTopologyProtons::kTop2proton]; 
  }

  /// Return the >2 protons topology probability
  float Result::GetNprotonsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no protons probability
    // multi-output
    return fOutput[TFMultioutputs::protons][TFTopologyProtons::kTopNproton]; 
  }

  /// Return the 0 pions topology probability
  float Result::Get0pionsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pions probability
    // multi-output
    return fOutput[TFMultioutputs::pions][TFTopologyPions::kTop0pion]; 
  }

  /// Return the 1 pions topology probability
  float Result::Get1pionsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pions probability
    // multi-output
    return fOutput[TFMultioutputs::pions][TFTopologyPions::kTop1pion]; 
  }

  /// Return the 2 pions topology probability
  float Result::Get2pionsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pions probability
    // multi-output
    return fOutput[TFMultioutputs::pions][TFTopologyPions::kTop2pion]; 
  }

  /// Return the >2 pions topology probability
  float Result::GetNpionsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pions probability
    // multi-output
    return fOutput[TFMultioutputs::pions][TFTopologyPions::kTopNpion]; 
  }

  /// Return the 0 pizeros topology probability
  float Result::Get0pizerosProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pizeros probability
    // multi-output
    return fOutput[TFMultioutputs::pizeros][TFTopologyPizeros::kTop0pizero]; 
  }

  /// Return the 1 pizeros topology probability
  float Result::Get1pizerosProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pizeros probability
    // multi-output
    return fOutput[TFMultioutputs::pizeros][TFTopologyPizeros::kTop1pizero]; 
  }

  /// Return the 2 pizeros topology probability
  float Result::Get2pizerosProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pizeros probability
    // multi-output
    return fOutput[TFMultioutputs::pizeros][TFTopologyPizeros::kTop2pizero]; 
  }

  /// Return the >2 pizeros topology probability
  float Result::GetNpizerosProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no pizeros probability
    // multi-output
    return fOutput[TFMultioutputs::pizeros][TFTopologyPizeros::kTopNpizero]; 
  }

  /// Return the 0 neutrons topology probability
  float Result::Get0neutronsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no neutrons probability
    // multi-output
    return fOutput[TFMultioutputs::neutrons][TFTopologyNeutrons::kTop0neutron]; 
  }

  /// Return the 1 neutrons topology probability
  float Result::Get1neutronsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no neutrons probability
    // multi-output
    return fOutput[TFMultioutputs::neutrons][TFTopologyNeutrons::kTop1neutron]; 
  }

  /// Return the 2 neutrons topology probability
  float Result::Get2neutronsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no neutrons probability
    // multi-output
    return fOutput[TFMultioutputs::neutrons][TFTopologyNeutrons::kTop2neutron]; 
  }

  /// Return the >2 neutrons topology probability
  float Result::GetNneutronsProbability() const
  {
    // single-output
    if(fOutput.size() == 1)
      return -1; // There is no neutrons probability
    // multi-output
    return fOutput[TFMultioutputs::neutrons][TFTopologyNeutrons::kTopNneutron]; 
  }
}
