////////////////////////////////////////////////////////////////////////
// Class:       CheckHits
// Plugin Type: analyzer (Unknown Unknown)
// File:        CheckHits_module.cc
//
// Generated Tue Sep 12 2023. Module from Alex Wilkinson to examine horizontal lines of hits observed in wiremod vs detvar variations
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"

#include <map>

class CheckFDHits : public art::EDAnalyzer {
public:
  explicit CheckFDHits(fhicl::ParameterSet const& p);

  CheckFDHits(CheckFDHits const&) = delete;
  CheckFDHits(CheckFDHits&&) = delete;
  CheckFDHits& operator=(CheckFDHits const&) = delete;
  CheckFDHits& operator=(CheckFDHits&&) = delete;

  void analyze(art::Event const& e) override;

  void beginJob() override;
  void endJob() override;

  void reset();

private:
  //const geo::GeometryCore* fGeom;
  const geo::Geometry* fGeom;
  unsigned int fTPCIndex;
  unsigned int fPlaneIndex;
  std::vector<unsigned int> fChsLocal;
  std::vector<unsigned int> fChsGlobal;

  std::string fSCLabel;
  std::string fWireLabel;
  std::string fHitLabel;

  TTree*                           fTreeWires;
  std::vector<std::vector<double>> fWires;
  std::vector<int>                 fWiresChs;
  std::vector<std::vector<double>> fSC;
  std::vector<int>                 fSCChs;
  std::vector<std::vector<double>> fHits;
  std::vector<int>                 fHitsChs;
  std::map<int, int>               fChTypes; // 0 induction, 1 collection
};


CheckFDHits::CheckFDHits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fTPCIndex   (p.get<unsigned int> ("TPCIndex")),
    fPlaneIndex (p.get<unsigned int> ("PlaneIndex")),
    fChsLocal   (p.get<std::vector<unsigned int>> ("ChannelsLocal")),
    fSCLabel    (p.get<std::string>  ("SCLabel")),
    fWireLabel  (p.get<std::string>  ("WireLabel")),
    fHitLabel   (p.get<std::string>  ("HitLabel"))
{
  consumes<std::vector<sim::SimChannel>>(fSCLabel);
  consumes<std::vector<recob::Wire>>(fWireLabel);
  consumes<std::vector<recob::Hit>>(fHitLabel);

  art::ServiceHandle<art::TFileService> tfs;

  fTreeWires = tfs->make<TTree>("sc_wire_hit", "sc_wire_hit");
  fTreeWires->Branch("simchannels", &fSC);
  fTreeWires->Branch("simchannels_chs", &fSCChs);
  fTreeWires->Branch("wires", &fWires);
  fTreeWires->Branch("wires_chs", &fWiresChs);
  fTreeWires->Branch("hits", &fHits);
  fTreeWires->Branch("hits_chs", &fHitsChs);
  fTreeWires->Branch("ch_types", &fChTypes);
}

void CheckFDHits::analyze(art::Event const& e)
{
  this->reset();

  art::Handle<std::vector<sim::SimChannel>> simChs;
  e.getByLabel(fSCLabel, simChs);
  art::Handle<std::vector<recob::Wire>> wires;
  e.getByLabel(fWireLabel, wires);
  art::Handle<std::vector<recob::Hit>> hits;
  e.getByLabel(fHitLabel, hits);

  for (sim::SimChannel simCh : *simChs) {
    raw::ChannelID_t ch = simCh.Channel();
    if (std::find(fChsGlobal.begin(), fChsGlobal.end(), ch) == fChsGlobal.end()) {
      continue;
    }
    std::vector<double> respVec(6000, 0.0); // Might be 4492 for VD, not sure
    for (sim::TDCIDE tickIde : simCh.TDCIDEMap()) {
      int tick = (int)tickIde.first;
      for (sim::IDE ide : tickIde.second) {
        respVec[tick] += (double)ide.energy;
      }
    }
    //fChTypes[(int)ch] = fGeom->SignalType(ch);
    bool isCollection = (fPlaneIndex == 2);
    fChTypes[(int)ch] = isCollection ? 1 : 0;
    fSC.push_back(respVec);
    fSCChs.push_back((int)ch);
  }
  for (recob::Wire wire : *wires) {
    raw::ChannelID_t ch = wire.Channel();
    if (std::find(fChsGlobal.begin(), fChsGlobal.end(), ch) == fChsGlobal.end()) {
      continue;
    }
    std::vector<double> respVec(6000, 0.0); // Might be 4492 for VD, not sure
    int tick = 0;
    for (float mag : wire.Signal()) {
      if (mag != 0) {
        respVec[tick] += (double)mag;
      }
      tick++;
    }
    //fChTypes[(int)ch] = fGeom->SignalType(ch);
    
    bool isCollection = (fPlaneIndex == 2);
    fChTypes[(int)ch] = isCollection ? 1 : 0;
    fWires.push_back(respVec);
    fWiresChs.push_back((int)ch);
  }

  std::map<unsigned int, std::vector<double>> hitResps;
  for (unsigned int ch : fChsGlobal) {
    hitResps[ch] = std::vector<double>(6000, 0.0);
  }
  for (recob::Hit hit : *hits) {
    raw::ChannelID_t ch = hit.Channel();
    if (std::find(fChsGlobal.begin(), fChsGlobal.end(), ch) == fChsGlobal.end()) {
      continue;
    }
    for (int tick = hit.StartTick(); tick <= hit.EndTick(); tick++) {
      hitResps[ch][tick] += hit.PeakAmplitude();
    }
  }
  for (std::pair<unsigned int, std::vector<double>> hitResp : hitResps) {
    bool isCollection = (fPlaneIndex == 2);
    fChTypes[(int)hitResp.first] = isCollection ? 1 : 0;
    //fChTypes[(int)hitResp.first] = fGeom->SignalType(hitResp.first);
    fHits.push_back(hitResp.second);
    fHitsChs.push_back((int)hitResp.first);
  }

  fTreeWires->Fill();
}

void CheckFDHits::beginJob()
{
  //fGeom = art::ServiceHandle<geo::Geometry>()->provider();
  fGeom = art::ServiceHandle<geo::Geometry>().get();
  const geo::CryostatID cID(0);
  const geo::TPCID tID(cID, fTPCIndex);
  const geo::PlaneID pID(tID, fPlaneIndex);
  //const readout::ROPID rID = fGeom->WirePlaneToROP(pID);
  std::cout << "TPC " << fTPCIndex
            << " Plane " << fPlaneIndex << " ("
            //<< (fGeom->SignalType(pID) == geo::SigType_t::kCollection ? "collection" : "induction")
            //<< ") First Channel in ROP " << fGeom->FirstChannelInROP(rID)
            << "\nChannels:\n";
  for (unsigned int chLocal : fChsLocal) { 
    unsigned int chGlobal = chLocal; //+ fGeom->FirstChannelInROP(rID);
    fChsGlobal.push_back(chGlobal);
    std::cout << chLocal << " -> " << chGlobal << "\n";
  }
}

void CheckFDHits::endJob()
{
}

void CheckFDHits::reset()
{
  fSC.clear();
  fSCChs.clear();
  fWires.clear();
  fWiresChs.clear();
  fHits.clear();
  fHitsChs.clear();
  fChTypes.clear();
}

DEFINE_ART_MODULE(CheckFDHits)
