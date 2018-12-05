////////////////////////////////////////////////////////////////////////
// \file    CVNEventDump_module.cc
// \brief   Analyzer module for creating CVN PixelMap objects
// \author  Alexander Radovic - a.radovic@gmail.com
//          Leigh Whitehead - leigh.howard.whitehead@cern.ch
//           - Added in truth based fiducial volume cuts
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
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
//#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"

#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/TrainingData.h"
#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/FDSensOpt/MVAAlg/MVAAlg.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"



namespace cvn {
  class CVNEventDump : public art::EDAnalyzer {
  public:

    explicit CVNEventDump(fhicl::ParameterSet const& pset);
    ~CVNEventDump();

    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);
    void beginJob() override;
    void endJob() override;

  private:

    dunemva::MVAAlg fMVAAlg;

    std::string fPixelMapInput;
    std::string fGenieGenModuleLabel;
    std::string fEnergyNueLabel;
    std::string fEnergyNumuLabel;
    bool        fGetEnergyOutput;
    bool        fGetEventWeight;
    bool        fWriteMapTH2;
    bool        fApplyFidVol;
    bool        fUseTopology;
    unsigned int fTopologyHits; // Number of hits for a track to be considered detectable 
                                   // for topology definitions.

    TrainingData* fTrain;
    TTree*        fTrainTree;

    //art::ServiceHandle<cheat::BackTracker> fBT;

    /// Function to extract TH2 from PixelMap and write to TFile
    void WriteMapTH2(const art::Event& evt, int slice, const PixelMap& pm);

  };

  //.......................................................................
  CVNEventDump::CVNEventDump(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset),
    fMVAAlg(pset.get<fhicl::ParameterSet>("MVAAlg"))
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNEventDump::~CVNEventDump()
  {  }

  //......................................................................
  void CVNEventDump::reconfigure(const fhicl::ParameterSet& pset)
  {
    fPixelMapInput  = pset.get<std::string>("PixelMapInput");
    fGenieGenModuleLabel  = pset.get<std::string>("GenieGenModuleLabel");
    fWriteMapTH2    = pset.get<bool>       ("WriteMapTH2");
    fApplyFidVol    = pset.get<bool>("ApplyFidVol");
    fEnergyNueLabel = pset.get<std::string> ("EnergyNueLabel");  
    fEnergyNumuLabel = pset.get<std::string> ("EnergyNumuLabel");
    fGetEnergyOutput = pset.get<bool> ("GetEnergyOutput");
    fGetEventWeight = pset.get<bool> ("GetEventWeight");
    fUseTopology  = pset.get<bool>("UseTopology");
    fTopologyHits = pset.get<unsigned int>("TopologyHitsCut");    
  }

  //......................................................................
  void CVNEventDump::beginJob()
  {


    art::ServiceHandle<art::TFileService> tfs;

    fTrainTree = tfs->make<TTree>("CVNTrainTree", "Training records");
    fTrainTree->Branch("train", "cvn::TrainingData", &fTrain);


  }

  //......................................................................
  void CVNEventDump::endJob()
  {

  }

  //......................................................................
  void CVNEventDump::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    art::Handle< std::vector< cvn::PixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cvn::PixelMap > > pixelmaplist;
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmapListHandle))
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    // If we have no pixel map then return
    if(pixelmaplist.size() == 0) return;

    InteractionType interaction = kOther;

    // * monte carlo
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if(evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
      art::fill_ptr_vector(mclist, mctruthListHandle);

    //unsigned short nmclist=  mclist.size();

    //std::cout<<"mctruth: "<<nmclist<<std::endl;

    art::Ptr<simb::MCTruth> truth = mclist[0];
    simb::MCNeutrino truthN=truth->GetNeutrino();
    //truth = mclist[0];

    AssignLabels labels;

    interaction = labels.GetInteractionType(truthN);
    if(fUseTopology){
      labels.GetTopology(truth,fTopologyHits);
//      labels.PrintTopology();
    }
    float nuEnergy = 0;
    float lepEnergy = 0;
//    if(truth.NeutrinoSet()){
      nuEnergy = truthN.Nu().E();
      lepEnergy = truthN.Lepton().E();
//    }

    // If outside the fiducial volume don't waste any time filling other variables
    if(fApplyFidVol){
      // Get the interaction vertex from the end point of the neutrino. This is 
      // because the start point of the lepton doesn't make sense for taus as they
      // are decayed by the generator and not GEANT
      TVector3 vtx = truthN.Nu().EndPosition().Vect();
      bool isFid = (fabs(vtx.X())<310 && fabs(vtx.Y())<550 && vtx.Z()>50 && vtx.Z()<1244);
      if(!isFid) return;
    }

    float recoNueEnergy = 0.;
    float recoNumuEnergy = 0.;
    // Should we use the EnergyReco_module reconstructed energies?
    if(fGetEnergyOutput){
      // Get the nue info
      if(fEnergyNueLabel != ""){
        art::Handle<dune::EnergyRecoOutput> energyRecoNueHandle;
        evt.getByLabel(fEnergyNueLabel, energyRecoNueHandle);

        recoNueEnergy = energyRecoNueHandle->fNuLorentzVector.E();
      }
      // And the numu
      if(fEnergyNumuLabel != ""){
        art::Handle<dune::EnergyRecoOutput> energyRecoNumuHandle;
        evt.getByLabel(fEnergyNumuLabel, energyRecoNumuHandle);

        recoNumuEnergy = energyRecoNumuHandle->fNuLorentzVector.E();
      }
    }

    // If we don't want to get the event weight then leave it as 1.0.
    double eventWeight = 1.;
    if(fGetEventWeight){
      double mvaResult = 0.; // We don't care about this, but need it for the call
      fMVAAlg.Run(evt,mvaResult,eventWeight);
    }

    // Create the training data and add it to the tree
    TrainingData train(interaction, nuEnergy, lepEnergy, recoNueEnergy, recoNumuEnergy, eventWeight, *pixelmaplist[0]);
    // Set the topology information
    int topPDG     = labels.GetPDG();
    int nprot      = labels.GetNProtons();
    int npion      = labels.GetNPions();
    int npi0       = labels.GetNPizeros();
    int nneut      = labels.GetNNeutrons();
    int toptype    = labels.GetTopologyType();
    int toptypealt = labels.GetTopologyTypeAlt();
    if(fUseTopology){
      train.SetTopologyInformation(topPDG, nprot, npion, npi0, nneut, toptype, toptypealt);
    }
    fTrain = &train;
    fTrainTree->Fill();

    // Make a plot of the pixel map if required
    if (fWriteMapTH2) WriteMapTH2(evt, 0, train.fPMap);

  }

  //----------------------------------------------------------------------



  void CVNEventDump::WriteMapTH2(const art::Event& evt, int slice, const PixelMap& pm)
  {
      std::stringstream name;
      name << "PixelMap_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameL;
      nameL << "PixelTruthMap_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameX;
      nameX << "PixelMap_X_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameY;
      nameY << "PixelMap_Y_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameZ;
      nameZ << "PixelMap_Z_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      TH2F* hist  = pm.ToTH2();
      TH2F* histL = pm.ToLabTH2();
      TH2F* histX = pm.SingleViewToTH2(0);
      TH2F* histY = pm.SingleViewToTH2(1);
      TH2F* histZ = pm.SingleViewToTH2(2);
      hist->SetName(name.str().c_str());
      histL->SetName(nameL.str().c_str());
      histX->SetName(nameX.str().c_str());
      histY->SetName(nameY.str().c_str());
      histZ->SetName(nameZ.str().c_str());

      art::ServiceHandle<art::TFileService> tfs;

      TH2F* histWrite = tfs->make<TH2F>(*hist);
      histWrite->Write();
      TH2F* histWriteL = tfs->make<TH2F>(*histL);
      histWriteL->GetZaxis()->SetRangeUser(0,10);
      histWriteL->Write();
      TH2F* histWriteX = tfs->make<TH2F>(*histX);
      histWriteX->Write();
      TH2F* histWriteY = tfs->make<TH2F>(*histY);
      histWriteY->Write();
      TH2F* histWriteZ = tfs->make<TH2F>(*histZ);
      histWriteZ->Write();

      delete hist;
      delete histWrite;
      delete histL;
      delete histWriteL;
      delete histX;
      delete histWriteX;
      delete histY;
      delete histWriteY;
      delete histZ;
      delete histWriteZ;

  }

DEFINE_ART_MODULE(cvn::CVNEventDump)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







