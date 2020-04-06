////////////////////////////////////////////////////////////////////////
/// \file    CTPHelper.h
/// \brief   Functions to help use the convolutional track PID
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#ifndef CTPHELPER_H
#define CTPHELPER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcorealg/GeoAlgo/GeoAlgo.h"

#include "dune/TrackPID/CTPResult.h"

namespace ctp
{

  /// Class containing some utility functions for all things CVN
  class CTPHelper
  {
  public:
    CTPHelper(const fhicl::ParameterSet& pset);
    ~CTPHelper();

    // Function to calculate the PID for a given track
    const ctp::CTPResult RunConvolutionalTrackPID(const art::Ptr<recob::PFParticle> particle, const art::Event &evt) const;

    // Calculate the features for the track PID
    const std::vector<std::vector<float>> GetNetworkInputs(const art::Ptr<recob::PFParticle>, const art::Event &evt) const;
    const std::vector<float> GetDeDxVector(const art::Ptr<recob::PFParticle>, const art::Event &evt) const;
    const std::vector<float> GetVariableVector(const art::Ptr<recob::PFParticle>, const art::Event &evt) const;

    // Get true PDG code for training
    const std::pair<const simb::MCParticle*,float> GetTrueParticle(const art::Ptr<recob::PFParticle>, const art::Event &evt) const;
    const int GetTruePDGCode (const art::Ptr<recob::PFParticle>, const art::Event &evt) const;

  private:

    void SmoothDedxVector(std::vector<float> &dedx) const;
    void PadDedxVector(std::vector<float> &dedx, const float mean, const float sigma) const;
    void GetDedxMeanAndSigma(const std::vector<float> &dedx, float &mean, float &sigma) const;
    void GetDeflectionMeanAndSigma(const art::Ptr<recob::Track> track, float &mean, float &sigma) const;

    void GetChildParticles(const art::Ptr<recob::PFParticle> part, const art::Event &evt, float &nTrack, float &nShower, float &nGrand) const;

    void NormaliseInputs(std::vector<std::vector<float>> &netInputs) const;

    // Variables for accessing the network architecture
    std::string fNetDir;
    std::string fNetName;

    // Module names
    std::string fParticleLabel; 
    std::string fTrackLabel;
    std::string fShowerLabel;
    std::string fCalorimetryLabel;

    // Parameters for variable extraction
    unsigned int fMinTrackPoints;
    unsigned int fDedxLength;
    float fQMax; // Maximum allowed charge
    float fQJump; // Maximum difference between consequetive dEdx values

    bool fNormalise; // Normalise the inputs for the network
  };

}

#endif  // CTPHELPER_H
