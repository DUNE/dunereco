////////////////////////////////////////////////////////////////////////
// \file    CVNEventDumpProtoDUNE_module.cc
// \brief   Analyzer module for creating CVN PixelMap objects using protoDUNE particles
// \author  Leigh Whitehead - leigh.howard.whitehead@cern.ch
// 
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


namespace cvn {
  class CVNEventDumpProtoDUNE : public art::EDAnalyzer {
  public:

    explicit CVNEventDumpProtoDUNE(fhicl::ParameterSet const& pset);
    ~CVNEventDumpProtoDUNE();

    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);
    void beginJob() override;
    void endJob() override;

  private:

    std::string fPixelMapInput;
    bool        fWriteMapTH2;
    bool        fApplyFidVol;

    TrainingData* fTrain;
    TTree*        fTrainTree;

    /// Function to extract TH2 from PixelMap and write to TFile
    void WriteMapTH2(const art::Event& evt, int slice, const PixelMap& pm);

  };

  //.......................................................................
  CVNEventDumpProtoDUNE::CVNEventDumpProtoDUNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNEventDumpProtoDUNE::~CVNEventDumpProtoDUNE()
  {  }

  //......................................................................
  void CVNEventDumpProtoDUNE::reconfigure(const fhicl::ParameterSet& pset)
  {
    fPixelMapInput  = pset.get<std::string>("PixelMapInput");
    fWriteMapTH2    = pset.get<bool>       ("WriteMapTH2");
  }

  //......................................................................
  void CVNEventDumpProtoDUNE::beginJob()
  {


    art::ServiceHandle<art::TFileService> tfs;

    fTrainTree = tfs->make<TTree>("CVNTrainTree", "Training records");
    fTrainTree->Branch("train", "cvn::TrainingData", &fTrain);


  }

  //......................................................................
  void CVNEventDumpProtoDUNE::endJob()
  {

  }

  //......................................................................
  void CVNEventDumpProtoDUNE::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    art::Handle< std::vector< cvn::PixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cvn::PixelMap > > pixelmaplist;
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmapListHandle))
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    std::cout << "Found " << pixelmaplist.size() << " pixel maps in event" << std::endl;

    for(unsigned int p = 0; p < pixelmaplist.size(); ++p){

      // We will have to just fake the truth information we would usually have for the events
      InteractionType interaction = kOther;

      // Create the training data and add it to the tree
      TrainingData train(interaction, 0.0, 0.0, 0.0, 0.0, 1.0, *pixelmaplist[p]);
      fTrain = &train;
      fTrainTree->Fill();

      // Make a plot of the pixel map if required
      if (fWriteMapTH2) WriteMapTH2(evt, p, train.fPMap);
    
    }

  }

  //----------------------------------------------------------------------



  void CVNEventDumpProtoDUNE::WriteMapTH2(const art::Event& evt, int slice, const PixelMap& pm)
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

DEFINE_ART_MODULE(cvn::CVNEventDumpProtoDUNE)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







