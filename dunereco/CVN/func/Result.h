////////////////////////////////////////////////////////////////////////
/// \file    Result.h
/// \brief   Result for CVN
/// \author  Alexander Radovic - a.radovic@gmail.com
///          Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CVN_RESULT_H
#define CVN_RESULT_H

#include <math.h>
#include <vector>
#include "dune/CVN/func/InteractionType.h"

namespace cvn
{
  /// Result, basic output of CVN neural net
  class Result
  {
  public:
    Result(const float* output, unsigned int& nOutputs);
    // Vector version of the constructor
    Result(const std::vector< std::vector<float> > output);
    //Result(const float* output, unsigned int& nOutputs, const float* features, unsigned int& nFeatures);
    Result();

    /// Index of maximum value in vector
    unsigned int ArgMax(int output_n) const;

    /// Maximum value in vector
    float Max();

    /// Return the predicted is_antineutrino
    TFIsAntineutrino PredictedIsAntineutrino() const;
 
    /// Return the predicted flavour
    TFFlavour PredictedFlavour() const;
    
    /// Return the predicted interaction
    TFInteraction PredictedInteraction() const;

    /// Return the predicted protons
    TFTopologyProtons PredictedProtons() const;

    /// Return the predicted pions
    TFTopologyPions PredictedPions() const;
 
    /// Return the predicted pizeros
    TFTopologyPizeros PredictedPizeros() const;

    /// Return the predicted neutrons
    TFTopologyNeutrons PredictedNeutrons() const;
 
    /// Return the is_antineutrino probability
    float GetIsAntineutrinoProbability() const;

    /// Return the numu flavour probability
    float GetNumuProbability() const;

    /// Return the nue flavour probability
    float GetNueProbability() const;

    /// Return the nutau flavour probability
    float GetNutauProbability() const;

    /// Return the NC probability
    float GetNCProbability() const;

    /// Return the CC QE interaction probability
    float GetQEProbability() const;

    /// Return the CC Res interaction probability
    float GetResProbability() const;

    /// Return the CC DIS interaction probability
    float GetDISProbability() const;

    /// Return the CC Other interaction probability
    float GetOtherProbability() const;

    /// Return the 0 protons topology probability
    float Get0protonsProbability() const;

    /// Return the 1 protons topology probability
    float Get1protonsProbability() const;

    /// Return the 2 protons topology probability
    float Get2protonsProbability() const;

    /// Return the >2 protons topology probability
    float GetNprotonsProbability() const;

    /// Return the 0 pions topology probability
    float Get0pionsProbability() const;

    /// Return the 1 pions topology probability
    float Get1pionsProbability() const;

    /// Return the 2 pions topology probability
    float Get2pionsProbability() const;

    /// Return the >2 pions topology probability
    float GetNpionsProbability() const;

    /// Return the 0 pizeros topology probability
    float Get0pizerosProbability() const;

    /// Return the 1 pizeros topology probability
    float Get1pizerosProbability() const;

    /// Return the 2 pizeros topology probability
    float Get2pizerosProbability() const;

    /// Return the >2 pizeros topology probability
    float GetNpizerosProbability() const;

    /// Return the 0 neutrons topology probability
    float Get0neutronsProbability() const;

    /// Return the 1 neutrons topology probability
    float Get1neutronsProbability() const;

    /// Return the 2 neutrons topology probability
    float Get2neutronsProbability() const;

    /// Return the >2 neutrons topology probability
    float GetNneutronsProbability() const;

    /// Number of outputs, i.e. size of vector
    unsigned int NOutput();

    std::vector< std::vector<float> > fOutput;  ///< Vector of outputs from neural net

  };
}

#endif  // CVN_RESULT_H

