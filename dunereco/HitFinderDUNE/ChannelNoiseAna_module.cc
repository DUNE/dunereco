////////////////////////////////////////////////////////////////////////
// Class:       ChannelNoiseAna
// Plugin Type: analyzer (art v2_04_00)
// File:        ChannelNoiseAna_module.cc
//
// Generated at Wed Oct  5 02:50:48 2016 by Matthew Thiesse using cetskelgen
// from cetlib version v1_20_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcore/Geometry/Geometry.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "dune/RunHistory/DetPedestalDUNE.h"

#include "RMSHitFinderAlg.h"

#include "TTree.h"

#include <vector>
#include <string>

namespace dune {
  class ChannelNoiseAna;
}


class dune::ChannelNoiseAna : public art::EDAnalyzer {
public:
  explicit ChannelNoiseAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ChannelNoiseAna(ChannelNoiseAna const &) = delete;
  ChannelNoiseAna(ChannelNoiseAna &&) = delete;
  ChannelNoiseAna & operator = (ChannelNoiseAna const &) = delete;
  ChannelNoiseAna & operator = (ChannelNoiseAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  std::string fWireModuleLabel;

  TTree * fTree;
  int run;
  int event;
  int channel;
  int wire;
  int tpc;
  float mean;
  float rms;
  float pedmean;
  float pedrms;

  dune::RMSHitFinderAlg fHitFinderAlg;

  art::ServiceHandle<geo::Geometry> fGeom;
  const lariov::DetPedestalProvider& fPedestalRetrievalAlg = *(lar::providerFrom<lariov::DetPedestalService>());

};


dune::ChannelNoiseAna::ChannelNoiseAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) ,
  fWireModuleLabel(p.get<std::string>("WireModuleLabel")) ,
  fHitFinderAlg(p.get<fhicl::ParameterSet>("RMSHitFinderAlg"))
{}

void dune::ChannelNoiseAna::analyze(art::Event const & e)
{
  run = e.run();
  event = e.event();

  art::Handle<std::vector<recob::Wire> > wireHandle;
  if (!e.getByLabel(fWireModuleLabel,wireHandle))
    {
      mf::LogError("ChannelNoiseAna") << "No recob::Wires found when expected.";
      return;
    }

  for (size_t iwire = 0; iwire < wireHandle->size(); ++iwire)
    {
      art::Ptr<recob::Wire> pwire(wireHandle,iwire);
      channel = pwire->Channel();
      wire = fGeom->ChannelToWire(pwire->Channel())[0].Wire;
      tpc = fGeom->ChannelToWire(pwire->Channel())[0].TPC;
      std::vector<float> signal = pwire->Signal();
      fHitFinderAlg.RobustRMSBase(signal,mean,rms);
      pedmean = fPedestalRetrievalAlg.PedMean(channel);
      pedrms = fPedestalRetrievalAlg.PedRms(channel);
      fTree->Fill();
    }
}

void dune::ChannelNoiseAna::beginJob()
{
  art::ServiceHandle<art::TFileService> fTfs;
  fTree = fTfs->make<TTree>("channoise","channoise");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("channel",&channel,"channel/I");
  fTree->Branch("wire",&wire,"wire/I");
  fTree->Branch("tpc",&tpc,"tpc/I");
  fTree->Branch("mean",&mean,"mean/F");
  fTree->Branch("rms",&rms,"rms/F");
  fTree->Branch("pedmean",&pedmean,"pedmean/F");
  fTree->Branch("pedrms",&pedrms,"pedrms/F");
}

DEFINE_ART_MODULE(dune::ChannelNoiseAna)
