////////////////////////////////////////////////////////////////////////
// \file    CVNValidation_module.cc
// \brief   Analyzer module to make some standard validation plots
//          of the CVN performance
// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include "TTree.h"
#include "TH2F.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/CVN/func/Result.h"

#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"


namespace cvn {
  class CVNValidation : public art::EDAnalyzer {
  public:

    explicit CVNValidation(fhicl::ParameterSet const& pset);
    ~CVNValidation();

    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);
    void beginJob() override;
    void endJob() override;

  private:

    std::string fResultLabel;
    std::string fTruthLabel;
    std::string fNumuEnergyLabel;
    std::string fNueEnergyLabel;

    TTree*        fValidTree;

    /// Tree branch variables
    /// Neutrino flavour probabilities
    float fNumuProb;
    float fNueProb;
    float fNutauProb;
    float fNCProb;
    /// Reco energy for numu and nue probabilities
    float fNumuEnergy;
    float fNueEnergy;
    /// Truth information
    int fPDG; // We set this to 0 for NC

  };

  //.......................................................................
  CVNValidation::CVNValidation(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNValidation::~CVNValidation()
  {  }

  //......................................................................
  void CVNValidation::reconfigure(const fhicl::ParameterSet& pset)
  {
    fResultLabel  = pset.get<std::string>("CVNResultLabel");
    fTruthLabel  = pset.get<std::string>("TruthLabel");
    fNumuEnergyLabel = pset.get<std::string> ("NumuEnergyLabel");
    fNueEnergyLabel = pset.get<std::string> ("NueEnergyLabel");  
  }

  //......................................................................
  void CVNValidation::beginJob()
  {

    art::ServiceHandle<art::TFileService> tfs;

    fValidTree = tfs->make<TTree>("CVNOutput", "CVN Output");
    fValidTree->Branch("pNumu",&fNumuProb,"fNumuProb/F");
    fValidTree->Branch("pNue",&fNueProb,"fNueProb/F");
    fValidTree->Branch("pNutau",&fNutauProb,"fNutauProb/F");
    fValidTree->Branch("pNC",&fNCProb,"fNCProb/F");
    fValidTree->Branch("numuEnergy",&fNumuEnergy,"fNumuEnergy/F");
    fValidTree->Branch("nueEnergy",&fNueEnergy,"fNueEnergy/F");
    fValidTree->Branch("pdg",&fPDG,"fPDG/I");
  }

  //......................................................................
  void CVNValidation::endJob()
  {

  }

  //......................................................................
  void CVNValidation::analyze(const art::Event& evt)
  {

    // Get the CVN results
    auto cvnResults = evt.getValidHandle<std::vector<cvn::Result>>(fResultLabel);
    if(!cvnResults.isValid()) return;
    if(cvnResults->size()==0) return;

    // Get the truth information
    auto truthInfo = evt.getValidHandle<std::vector<simb::MCTruth>>(fTruthLabel);    
    if(!truthInfo.isValid()) return;
    if(truthInfo->size()==0) return;
    if(!truthInfo->at(0).NeutrinoSet()) return;

    // Get the energy results
    auto numuEnergy = evt.getValidHandle<dune::EnergyRecoOutput>(fNumuEnergyLabel);     
    auto nueEnergy = evt.getValidHandle<dune::EnergyRecoOutput>(fNueEnergyLabel);     
    if(!numuEnergy.isValid() || !nueEnergy.isValid()) return;

    // Fill the CVN flavour branches 
    fNumuProb  = cvnResults->at(0).GetNumuProbability();    
    fNueProb   = cvnResults->at(0).GetNueProbability();    
    fNutauProb = cvnResults->at(0).GetNutauProbability();    
    fNCProb    = cvnResults->at(0).GetNCProbability();    
    // Reconstructed energy
    fNumuEnergy = numuEnergy->fNuLorentzVector.E();
    fNueEnergy  = nueEnergy->fNuLorentzVector.E();
    // PDG
    fPDG = truthInfo->at(0).GetNeutrino().Nu().PdgCode();
    if(truthInfo->at(0).GetNeutrino().CCNC() == simb::kNC){
      fPDG = 0;
    }

    fValidTree->Fill();

  }

DEFINE_ART_MODULE(cvn::CVNValidation)

} // end namespace cvn
////////////////////////////////////////////////////////////////////////







