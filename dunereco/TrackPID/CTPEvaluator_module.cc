/**
 *  @file   dune/trackPID/CTPEvaluator_module.cc
 *
 *  @brief  This module performs track PID using the convolutionsl 
 *          track PID network
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "TTree.h"
#include "TVector3.h"

#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"

#include "dune/TrackPID/CTPHelper.h"
#include "dune/TrackPID/CTPResult.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <fstream>
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace ctp
{

/**
 *  @brief  CTPEvaluator class
 */
class CTPEvaluator : public art::EDAnalyzer
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
     CTPEvaluator(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
     virtual ~CTPEvaluator();

     void beginJob();
     void endJob();
     void analyze(const art::Event &evt);

private:

  fhicl::ParameterSet fHelperPars;
  CTPHelper fConvTrackPID;
  std::string fParticleLabel;
  bool fWriteTree;

  std::vector<float> fMuonScoreVector;
  std::vector<float> fPionScoreVector;
  std::vector<float> fProtonScoreVector;
  std::vector<int>   fPDGVector;
  std::vector<unsigned int> fCaloPoints;


  TTree *fPIDTree;
};

DEFINE_ART_MODULE(CTPEvaluator)

} // namespace ctp

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "TRandom3.h"

#include <iostream>
#include <random>

namespace ctp
{

CTPEvaluator::CTPEvaluator(fhicl::ParameterSet const &pset) : art::EDAnalyzer(pset),
fHelperPars(pset.get<fhicl::ParameterSet>("ctpHelper")),
fConvTrackPID(fHelperPars),
fParticleLabel(pset.get<std::string>("particleLabel")),
fWriteTree(pset.get<bool>("writeTree"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CTPEvaluator::~CTPEvaluator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPEvaluator::beginJob()
{
    if (!fWriteTree) return;
  
    art::ServiceHandle<art::TFileService const> tfs;
  
    fPIDTree = tfs->make<TTree>("pidTree","pidTree");
    fPIDTree->Branch("muonScores",&fMuonScoreVector);
    fPIDTree->Branch("pionScores",&fPionScoreVector);
    fPIDTree->Branch("protonScores",&fProtonScoreVector);
    fPIDTree->Branch("pdgCodes",&fPDGVector);
    fPIDTree->Branch("caloPoints",&fCaloPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPEvaluator::endJob()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CTPEvaluator::analyze(const art::Event &evt)
{
    if (fWriteTree)
    {
        fMuonScoreVector.clear();
        fPionScoreVector.clear();
        fProtonScoreVector.clear();
        fPDGVector.clear();
        fCaloPoints.clear();
    }

    // Get all of the PFParticles
    const std::vector<art::Ptr<recob::PFParticle>> particles = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt,fParticleLabel);

    for (const art::Ptr<recob::PFParticle> &particle : particles)
    {
        // Get the track if this particle is track-like
        unsigned int nCaloPoints = 0;
        const std::string trkLabel = fHelperPars.get<std::string>("TrackLabel");
        if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(particle,evt,fParticleLabel,trkLabel))
        {
            const art::Ptr<recob::Track> trk = dune_ana::DUNEAnaPFParticleUtils::GetTrack(particle,evt,fParticleLabel,trkLabel);
            const std::string caloLabel = fHelperPars.get<std::string>("CalorimetryLabel");
            const art::Ptr<anab::Calorimetry> calo = dune_ana::DUNEAnaTrackUtils::GetCalorimetry(trk,evt,trkLabel,caloLabel);

            nCaloPoints = calo->dEdx().size();
        }
        else continue;

        // Returns a dummy value if not a track or not suitable
        CTPResult thisPID = fConvTrackPID.RunConvolutionalTrackPID(particle,evt);
        if (!thisPID.IsValid()) continue;

        const int pdg = fConvTrackPID.GetTruePDGCode(particle,evt);

        std::cout << "Got a track PID for particle of type " << pdg << ": " << thisPID.GetMuonScore() << ", " << thisPID.GetPionScore() << ", " << thisPID.GetProtonScore() << std::endl;       

        if (fWriteTree)
        {
            fMuonScoreVector.push_back(thisPID.GetMuonScore());
            fPionScoreVector.push_back(thisPID.GetPionScore());
            fProtonScoreVector.push_back(thisPID.GetProtonScore());
            fPDGVector.push_back(pdg);
            fCaloPoints.push_back(nCaloPoints);
        }
    }
    if (fWriteTree) fPIDTree->Fill();
}

} //namespace ctp

