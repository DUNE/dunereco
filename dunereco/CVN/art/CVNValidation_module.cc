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
#include "TVector3.h"

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

#include "dunereco/CVN/func/Result.h"

#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"


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

    // std::string fResultLabel;
    std::string fTruthLabel;
    // std::string fNumuEnergyLabel;
    // std::string fNueEnergyLabel;

    TTree*        fValidTree;

    /// Tree branch variables
    /// Neutrino flavour probabilities
    // float fNumuProb;
    // float fNueProb;
    // float fNutauProb;
    // float fNCProb;
    /// Reco energy for numu and nue probabilities
    // float fNumuEnergy;
    // float fNueEnergy;
    /// Truth information
    int fPDG; // We set this to 0 for NC
    float fNuEnergy;
    float fLepEnergy;
    float fLepAngle;
    float fVtxX;
    float fVtxY;
    float fVtxZ;
    int fMode;

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
    // fResultLabel  = pset.get<std::string>("CVNResultLabel");
    fTruthLabel  = pset.get<std::string>("TruthLabel");
    // fNumuEnergyLabel = pset.get<std::string> ("NumuEnergyLabel");
    // fNueEnergyLabel = pset.get<std::string> ("NueEnergyLabel");
  }

  //......................................................................
  void CVNValidation::beginJob()
  {

    art::ServiceHandle<art::TFileService> tfs;

    fValidTree = tfs->make<TTree>("CVNOutput", "CVN Output");
    // fValidTree->Branch("pNumu",&fNumuProb,"fNumuProb/F");
    // fValidTree->Branch("pNue",&fNueProb,"fNueProb/F");
    // fValidTree->Branch("pNutau",&fNutauProb,"fNutauProb/F");
    // fValidTree->Branch("pNC",&fNCProb,"fNCProb/F");
    // fValidTree->Branch("numuEnergy",&fNumuEnergy,"fNumuEnergy/F");
    // fValidTree->Branch("nueEnergy",&fNueEnergy,"fNueEnergy/F");
    fValidTree->Branch("pdg",&fPDG,"fPDG/I");
    fValidTree->Branch("mode",&fMode,"fMode/I");
    fValidTree->Branch("vtxx",&fVtxX,"fVtxX/F");
    fValidTree->Branch("vtxy",&fVtxY,"fVtxY/F");
    fValidTree->Branch("vtxz",&fVtxZ,"fVtxZ/F");
    fValidTree->Branch("nuenergy",&fNuEnergy,"fNuEnergy/F");
    fValidTree->Branch("lepenergy",&fLepEnergy,"fLepEnergy/F");
    fValidTree->Branch("lepangle",&fLepAngle,"fLepAngle/F");
  }

  //......................................................................
  void CVNValidation::endJob()
  {

  }

  //......................................................................
  void CVNValidation::analyze(const art::Event& evt)
  {

    // Get the CVN results
    // auto cvnResults = evt.getValidHandle<std::vector<cvn::Result>>(fResultLabel);
    // if(!cvnResults.isValid()) return;
    // if(cvnResults->size()==0) return;

    // Get the truth information
    auto truthInfo = evt.getValidHandle<std::vector<simb::MCTruth>>(fTruthLabel);
    if(!truthInfo.isValid()) return;
    if(truthInfo->size()==0) return;
    if(!truthInfo->at(0).NeutrinoSet()) return;
    simb::MCNeutrino true_neutrino = truthInfo->at(0).GetNeutrino();

    // Get the energy results
    // auto numuEnergy = evt.getValidHandle<dune::EnergyRecoOutput>(fNumuEnergyLabel);
    // auto nueEnergy = evt.getValidHandle<dune::EnergyRecoOutput>(fNueEnergyLabel);
    // if(!numuEnergy.isValid() || !nueEnergy.isValid()) return;

    // Fill the CVN flavour branches
    // fNumuProb  = cvnResults->at(0).GetNumuProbability();
    // fNueProb   = cvnResults->at(0).GetNueProbability();
    // fNutauProb = cvnResults->at(0).GetNutauProbability();
    // fNCProb    = cvnResults->at(0).GetNCProbability();
    // Reconstructed energy
    // fNumuEnergy = numuEnergy->fNuLorentzVector.E();
    // fNueEnergy  = nueEnergy->fNuLorentzVector.E();
    // PDG
    fPDG = true_neutrino.Nu().PdgCode();
    if(true_neutrino.CCNC() == simb::kNC){
      fPDG = 0;
    }
    fMode = true_neutrino.Mode();
    fNuEnergy = true_neutrino.Nu().E();
    fLepEnergy = true_neutrino.Lepton().E();
    TVector3 nu_mom = true_neutrino.Nu().Momentum().Vect();
    TVector3 lep_mom = true_neutrino.Lepton().Momentum().Vect();
    fLepAngle = nu_mom.Angle(lep_mom);

    TVector3 vtx = true_neutrino.Nu().EndPosition().Vect();
    fVtxX = vtx.X();
    fVtxY = vtx.Y();
    fVtxZ = vtx.Z();

    fValidTree->Fill();


  }

DEFINE_ART_MODULE(cvn::CVNValidation)

} // end namespace cvn
////////////////////////////////////////////////////////////////////////
