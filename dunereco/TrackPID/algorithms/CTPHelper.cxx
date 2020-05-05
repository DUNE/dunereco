////////////////////////////////////////////////////////////////////////
/// \file    CTPHelper.cxx
/// \brief   Functions to help use the convolutional track PID
/// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <random>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "dune/TrackPID/algorithms/CTPHelper.h"
#include "dune/TrackPID/products/CTPResult.h"
#include "dune/TrackPID/tf/CTPGraph.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"

#include "cetlib/getenv.h"

namespace ctp
{

  // Constructor
  CTPHelper::CTPHelper(const fhicl::ParameterSet& pset){
    fNetDir  = pset.get<std::string>("NetworkPath");
    fNetName = pset.get<std::string>("NetworkName");
    fParticleLabel = pset.get<std::string>("ParticleLabel");
    fTrackLabel = pset.get<std::string>("TrackLabel");
    fShowerLabel = pset.get<std::string>("ShowerLabel");
    fCalorimetryLabel = pset.get<std::string>("CalorimetryLabel");
    fMinTrackPoints = pset.get<unsigned int>("MinHits",50);
    fDedxLength = pset.get<unsigned int>("DedxLength",100);
    fQMax = pset.get<float>("MaxCharge",1000);
    fQJump = pset.get<float>("MaxChargeJump",500);
    fNormalise = pset.get<bool>("NormaliseInputs",true);
  }

  CTPHelper::~CTPHelper(){
  }

  // Function to calculate the PID for a given track
  const ctp::CTPResult CTPHelper::RunConvolutionalTrackPID(const art::Ptr<recob::PFParticle> part, const art::Event &evt) const{

    // Get the inputs to the network
    std::vector< std::vector< std::vector<float> > > finalInputs;

    std::vector< std::vector<float> > twoVecs = GetNetworkInputs(part,evt);
    if(twoVecs.empty()) return ctp::CTPResult();

    finalInputs.push_back(twoVecs);

//    std::cout << "We have sensible inputs... building TF object" << std::endl;
    // Load the network and run it
    const std::string fullPath = cet::getenv(fNetDir) + "/" + fNetName;
    std::unique_ptr<tf::CTPGraph> convNet = tf::CTPGraph::create(fullPath.c_str(),std::vector<std::string>(),2,1);
//    std::cout << "Calling TF interface" << std::endl;
    std::vector< std::vector< std::vector<float> > > convNetOutput = convNet->run(finalInputs);

//    std::cout << "Got result from TF" << std::endl;
    ctp::CTPResult result(convNetOutput.at(0).at(0));
//    std::cout << "Returning CTPResult" << std::endl;
    return result;
  }

  // Calculate the features for the track PID
  const std::vector<std::vector<float>> CTPHelper::GetNetworkInputs(const art::Ptr<recob::PFParticle> part, const art::Event &evt) const{
  
    std::vector<std::vector<float>> netInputs;

    if(!dune_ana::DUNEAnaPFParticleUtils::IsTrack(part,evt,fParticleLabel,fTrackLabel)){
//      std::cout << "CTPHelper: this PFParticle is not track-like... returning empty vector." << std::endl;
      return std::vector<std::vector<float>>();
    }

    // Use the analysis utilities to simplify finding products and associations
    art::Ptr<recob::Track> thisTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(part,evt,fParticleLabel,fTrackLabel);
    art::Ptr<anab::Calorimetry> thisCalo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(thisTrack,evt,fTrackLabel,fCalorimetryLabel);

    if(thisCalo->dEdx().size() < fMinTrackPoints){
//      std::cout << "CTPHelper: this track has too few points for PID (" << thisCalo->dEdx().size() << ")... returning empty vector." << std::endl;
      return std::vector<std::vector<float>>();
    }

    std::vector<float> dedxVector = thisCalo->dEdx();
    this->SmoothDedxVector(dedxVector);
    float dedxMean = 0.;
    float dedxSigma = 0.;

    // We want to use the middle third of the dedx vector (with max length 100)
    std::vector<float> dedxTrunc;
    unsigned int pointsForAverage = (fDedxLength - fMinTrackPoints) / 3;
    unsigned int avStart = dedxVector.size() - 1 - pointsForAverage;
    unsigned int avEnd   = dedxVector.size() - 1 - (2*pointsForAverage);
    for(unsigned int e = avStart; e > avEnd; --e) dedxTrunc.push_back(dedxVector.at(e));
    this->GetDedxMeanAndSigma(dedxTrunc,dedxMean,dedxSigma);

    // If our dedx vector is between fMinTrackPoints and fDedxLength in size then we need to pad it
    if(dedxVector.size() < fDedxLength){
      this->PadDedxVector(dedxVector,dedxMean,dedxSigma);
    }
    
    std::vector<float> finalInputDedx;
    finalInputDedx.insert(finalInputDedx.begin(),dedxVector.end() - fDedxLength,dedxVector.end());

    std::vector<float> finalInputVariables;
    // Get the number of child particles
    float nTrack, nShower, nGrand;
    this->GetChildParticles(part,evt,nTrack,nShower,nGrand);
    finalInputVariables.push_back(nTrack);
    finalInputVariables.push_back(nShower);
    finalInputVariables.push_back(nGrand);
    // Now add the dedx mean and sigma
    finalInputVariables.push_back(dedxMean);
    finalInputVariables.push_back(dedxSigma);
    // Finally, get the angular deflection mean and sigma
    float deflectionMean, deflectionSigma;
    this->GetDeflectionMeanAndSigma(thisTrack,deflectionMean,deflectionSigma);
    finalInputVariables.push_back(deflectionMean);
    finalInputVariables.push_back(deflectionSigma);
  
    netInputs.push_back(finalInputDedx);
    netInputs.push_back(finalInputVariables);
   
    if(fNormalise) this->NormaliseInputs(netInputs);

    return netInputs;
  }

  const std::vector<float> CTPHelper::GetDeDxVector(const art::Ptr<recob::PFParticle> part, const art::Event &evt) const{
    return GetNetworkInputs(part,evt).at(0);
  }

  const std::vector<float> CTPHelper::GetVariableVector(const art::Ptr<recob::PFParticle> part, const art::Event &evt) const{
    return GetNetworkInputs(part,evt).at(1);
  }

  const std::pair<const simb::MCParticle*,float> CTPHelper::GetTrueParticle(const art::Ptr<recob::PFParticle> part, const art::Event &evt) const{
    // Get hits
    const std::vector<art::Ptr<recob::Hit>> collectionHits  = dune_ana::DUNEAnaPFParticleUtils::GetHits(part,evt,fParticleLabel);
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Function to find a weighted and sorted vector of the true particles matched to a particle
    using weightedMCPair = std::pair<const simb::MCParticle*, float>;
    std::vector<weightedMCPair> outVecHits;

    // Loop over all hits in the input vector and record the contributing MCParticles.
    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    std::unordered_map<const simb::MCParticle*, float> mcHitMap;
    float hitTotal = 0;
    for(const art::Ptr<recob::Hit> &hit : collectionHits) {
      for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(clockData, hit)) {
        const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
        mcHitMap[curr_part] += 1.;
        ++hitTotal;
      }
    }

    for (weightedMCPair const &p : mcHitMap) outVecHits.push_back(p);
    // Can't continue without a truth match
    if(outVecHits.size() == 0) assert(0);

    std::sort(outVecHits.begin(),outVecHits.end(),[](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});
    for(weightedMCPair &p : outVecHits) p.second /= hitTotal;

//    std::cout << "Returning truth match..." << std::endl;
    return outVecHits.at(0);
  }

  const int CTPHelper::GetTruePDGCode(const art::Ptr<recob::PFParticle> part, const art::Event &evt) const{
    return this->GetTrueParticle(part,evt).first->PdgCode();
  }
 
  void CTPHelper::SmoothDedxVector(std::vector<float> &dedx) const{

    // Firstly, get rid of all very high values > fQMax
    for(float &val : dedx){
//      std::cout << "Checking de/dx = " << val;
      if(val > fQMax){
//        std::cout << "CTPHelper: large value " << val << " set to " << fQMax << std::endl; 
        val = fQMax;
      }
      if(val < 1e-3){
        std::cout << "CTPHelper: small value " << val << " set to 1.0e-3 " << std::endl;
        val = 1e-3;
      }
//      std::cout << " value becomes " << val << std::endl;
    }

    // Now try to smooth over jumps
    unsigned int nQ = dedx.size();
    // Now do the rest of the points
    for(unsigned int q = 1; q < nQ - 1; ++q){
      if((dedx[q] - dedx[q-1]) > fQJump)
      {
//        std::cout << "Updating dedx point " << q << " from " << dedx[q] << " to " << 0.5 * (dedx[q-1]+dedx[q+1]) << std::endl;
        dedx[q] = 0.5 * (dedx[q-1]+dedx[q+1]);
      }
    }
    // First and last points are special cases
//    std::cout << "Initial first and last points: " << dedx[0] << ", " << dedx[nQ-1] << std::endl;
//    if((dedx[0] - dedx[1]) > fQJump) dedx[0] = dedx[1] + (dedx[1] - dedx[2]);
//    if((dedx[nQ-1] - dedx[nQ-2]) > fQJump) dedx[nQ-1] = dedx[nQ-2] + (dedx[nQ-2] - dedx[nQ-3]);
//    std::cout << "Final first and last points: " << dedx[0] << ", " << dedx[nQ-1] << std::endl;
//    std::cout << "Using... " << dedx[1] << ", " << dedx[nQ-3] << ", " << dedx[nQ-2] << std::endl;

  }

  void CTPHelper::PadDedxVector(std::vector<float> &dedx, const float mean, const float sigma) const{

    std::default_random_engine generator;
    std::normal_distribution<float> gaussDist(mean,sigma);

    unsigned int originalSize = dedx.size();

    for (unsigned int h = 0; h + originalSize < fDedxLength; ++h)
    {
      // Pick a random Gaussian value but ensure we don't go negative
      float randVal = -1;
      do
      {
        randVal = gaussDist(generator);
      }
      while (randVal < 0);
      // Pad from beginning to keep the real track part at the end
      dedx.insert(dedx.begin(),randVal);
    }

  } 

  void CTPHelper::GetDedxMeanAndSigma(const std::vector<float> &dedx, float &mean, float &sigma) const{

    float averageDedx = 0;
    float sigmaDedx = 0;
    for(const float &q : dedx) averageDedx += q;
    mean = averageDedx / static_cast<float>(dedx.size());
    for(const float &q : dedx) sigmaDedx += (mean -q)*(mean-q);
    sigma = std::sqrt(sigmaDedx / static_cast<float>(dedx.size()));
  }

  void CTPHelper::GetDeflectionMeanAndSigma(const art::Ptr<recob::Track> track, float &mean, float &sigma) const{

    std::vector<float> trajAngle;
    for(unsigned int p = 1; p < track->Trajectory().NPoints(); ++p){
      TVector3 thisDir = track->Trajectory().DirectionAtPoint<TVector3>(p);
      TVector3 prevDir = track->Trajectory().DirectionAtPoint<TVector3>(p-1);
      trajAngle.push_back(thisDir.Angle(prevDir));

//      if(p < 11 || p > track->Trajectory().NPoints() - 11) std::cout << "Deflection at point " << p << " = " << trajAngle.at(p-1) << std::endl;
    }

    // Average and sigma of the angular deflection between trajectory points (wobble)
    float averageAngle = 0;
    float standardDevAngle = 0;
    for(const float &a : trajAngle) averageAngle += a;
    mean = averageAngle / static_cast<float>(trajAngle.size());
    for(const float &a : trajAngle) standardDevAngle += (a - mean)*(a - mean);
    sigma = sqrt(standardDevAngle / static_cast<float>(trajAngle.size()));
  }

  void CTPHelper::GetChildParticles(const art::Ptr<recob::PFParticle> part, const art::Event &evt, float &nTrack, float &nShower, float &nGrand) const{

    nTrack = 0.;
    nShower = 0.;
    nGrand = 0.;

    std::vector<art::Ptr<recob::PFParticle>> children = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(part,evt,fParticleLabel);

    for(const art::Ptr<recob::PFParticle> child : children){
      nTrack += dune_ana::DUNEAnaPFParticleUtils::IsTrack(child,evt,fParticleLabel,fTrackLabel);
      nShower += dune_ana::DUNEAnaPFParticleUtils::IsShower(child,evt,fParticleLabel,fShowerLabel);
      nGrand += child->NumDaughters();
    }
//    std::cout << "Children = " << children.size() << "( " << nTrack << ", " << nShower << ") and grand children = " << nGrand << std::endl;
  }

  void CTPHelper::NormaliseInputs(std::vector<std::vector<float>> &inputs) const{

    // Firstly, we need to normalise the dE/dx values
    for(float &dedx : inputs.at(0)){
      // Protect against negative values
      if(dedx < 1.e-3){
//        std::cout << "Very low dedx value = " << dedx << std::endl;
        dedx = 1.e-3;
      }
      float newVal = std::log(2*dedx) + 1;
      if(newVal > 5.) newVal = 5.;
      dedx = (newVal / 2.5) - 1;
    }

    // For the three child values
    for(unsigned int v = 0; v < 3; ++v){
      float val = inputs.at(1).at(v);
      if(val > 5) val = 5;
      inputs.at(1).at(v) = (val / 2.5) - 1.0;
    }

    // The charge mean
    float val = inputs.at(1).at(3);
    val = std::log(2*val) + 1;
    if(val > 5.) val = 5.;
    inputs.at(1).at(3) = (val / 2.5) - 1;
    
    // The charge sigma
    val = inputs.at(1).at(4);
    if(val > 3) val = 3.;
    inputs.at(1).at(4) = (val/1.5) - 1; 

    // The angle mean
    val = inputs.at(1).at(5);
    if(val > 0.05) val = 0.05;
    inputs.at(1).at(5) = (val / 0.025) - 1.0;

    // The angle sigma
    val = inputs.at(1).at(6);
    if(val > 0.3) val = 0.3;
    inputs.at(1).at(6) = (val / 0.15) - 1.0;

  }

}
