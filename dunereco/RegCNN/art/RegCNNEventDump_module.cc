////////////////////////////////////////////////////////////////////////
// \file    RegCNNEventDump_module.cc
// \brief   Analyzer module for creating RegCNN PixelMap objects modified from CVNEventDump_module.cc
// \author  Ilsoo Seong - iseong@uci.edu
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
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larsim/MCCheater/BackTracker.h"

#include "dune/RegCNN/func/RegPixelMap.h"



namespace cnn {
  class RegCNNEventDump : public art::EDAnalyzer {
  public:

    explicit RegCNNEventDump(fhicl::ParameterSet const& pset);
    ~RegCNNEventDump();

    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);
    void beginJob() override;
    void endJob() override;

  private:

    std::string fPixelMapInput;
    std::string fGenieGenModuleLabel;
    bool        fWriteMapTH2;
    bool        fApplyFidVol;

    //art::ServiceHandle<cheat::BackTracker> fBT;

    /// Function to extract TH2 from PixelMap and write to TFile
    void WriteMapTH2(const art::Event& evt, int slice, const RegPixelMap& pm);

  };

  //.......................................................................
  RegCNNEventDump::RegCNNEventDump(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  RegCNNEventDump::~RegCNNEventDump()
  {  }

  //......................................................................
  void RegCNNEventDump::reconfigure(const fhicl::ParameterSet& pset)
  {
    fPixelMapInput  = pset.get<std::string>("PixelMapInput");
    fGenieGenModuleLabel  = pset.get<std::string>("GenieGenModuleLabel");
    fWriteMapTH2    = pset.get<bool>       ("WriteMapTH2");
    fApplyFidVol    = pset.get<bool>("ApplyFidVol");
    
  }

  //......................................................................
  void RegCNNEventDump::beginJob()
  {

  }

  //......................................................................
  void RegCNNEventDump::endJob()
  {

  }

  //......................................................................
  void RegCNNEventDump::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    art::Handle< std::vector< cnn::RegPixelMap > > pixelmapListHandle;
    std::vector< art::Ptr< cnn::RegPixelMap > > pixelmaplist;
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, pixelmapListHandle))
      art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    unsigned short nhits =  pixelmaplist.size();
    // If we have no hits then return
    if(nhits < 1) return;


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

    //float nuEnergy = 0;
    //float lepEnergy = 0;
    ////if(truth.NeutrinoSet()){
    //nuEnergy = truthN.Nu().E();
    //lepEnergy = truthN.Lepton().E();
    ////}

    // If outside the fiducial volume don't waste any time filling other variables
    if(fApplyFidVol){
      // Get the interaction vertex from the end point of the neutrino. This is 
      // because the start point of the lepton doesn't make sense for taus as they
      // are decayed by the generator and not GEANT
      TVector3 vtx = truthN.Nu().EndPosition().Vect();
      bool isFid = (fabs(vtx.X())<310 && fabs(vtx.Y())<550 && vtx.Z()>50 && vtx.Z()<1244);
      if(!isFid) return;
    }

    // Make a plot of the pixel map if required
    if (fWriteMapTH2) WriteMapTH2(evt, 0, *pixelmaplist[0]);

  }

  //----------------------------------------------------------------------



  void RegCNNEventDump::WriteMapTH2(const art::Event& evt, int slice, const RegPixelMap& pm)
  {
      std::stringstream name;
      name << "RegPixelMap_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameL;
      nameL << "PixelTruthMap_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameX;
      nameX << "RegPixelMap_X_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameY;
      nameY << "RegPixelMap_Y_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
      std::stringstream nameZ;
      nameZ << "RegPixelMap_Z_r" << evt.run() << "_s" << evt.subRun()<< "_e" << evt.event() << "_sl" << slice;
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

DEFINE_ART_MODULE(cnn::RegCNNEventDump)
} // end namespace cnn
////////////////////////////////////////////////////////////////////////







