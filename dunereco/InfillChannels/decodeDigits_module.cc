////////////////////////////////////////////////////////////////////////
// Class:       decodeDigits
// Plugin Type: analyzer (art v3_05_01)
// File:        decodeDigits_module.cc
//
// Generated at Wed Dec  9 11:34:56 2020 by Alexander Wilkinson using cetskelgen
// from cetlib version v3_10_00.
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

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/raw.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include <TH2F.h>
#include <TTree.h>
#include <TGeoVolume.h>
#include <TFile.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <utility>

namespace infill {
  class decodeDigits;
}


class infill::decodeDigits : public art::EDAnalyzer {
public:
  explicit decodeDigits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  decodeDigits(decodeDigits const&) = delete;
  decodeDigits(decodeDigits&&) = delete;
  decodeDigits& operator=(decodeDigits const&) = delete;
  decodeDigits& operator=(decodeDigits&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.
  const geo::GeometryCore* geom;

  art::ServiceHandle<art::TFileService> tfs;

  std::set<raw::ChannelID_t> badChannels;
  std::set<raw::ChannelID_t> noisyChannels;
  std::set<raw::ChannelID_t> deadChannels;

  std::set<readout::ROPID> activeRops;
};

infill::decodeDigits::decodeDigits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  consumes<std::vector<raw::RawDigit>>("daq");
}

void infill::decodeDigits::analyze(art::Event const& e)
{
  std::string sEvent = (
    "Ev" + std::to_string(e.id().event()) + "Run" + std::to_string(e.id().run()) + 
    "SRun" + std::to_string(e.id().subRun())
  );
  std::map<readout::ROPID, TH2F*> ropImages;

  // Prepare TH2s
  for (const readout::ROPID& rop : activeRops) {
    std::string sTitle = (
      sEvent + "_TPCset" + std::to_string(rop.TPCset) + "_ROP" + std::to_string(rop.ROP)
    );
    ropImages[rop] = tfs->make<TH2F>(
      sTitle.c_str(), sTitle.c_str(), geom->Nchannels(rop), 0, geom->Nchannels(rop), 6000, 0, 6000
    );
  }

  art::Handle<std::vector<raw::RawDigit>> digs;
  e.getByLabel("daq", digs);

  // Fill TH2s
  for(const raw::RawDigit& dig : *digs){ 
    readout::ROPID rop = geom->ChannelToROP(dig.Channel());
    if(ropImages.count(rop)){ // Ignores ROPs associated with only inactive TPCIDs
      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());

      const raw::ChannelID_t firstCh = geom->FirstChannelInROP(rop);
      for(unsigned int tick = 0; tick < adcs.size(); ++tick){
        const int adc = adcs[tick] ? int(adcs[tick]) - dig.GetPedestal() : 0;

        ropImages[rop]->Fill(dig.Channel() - firstCh, tick, adc);
      }
    }
  }

  tfs->file().Write();
  tfs->file().Flush();
  for(auto it: ropImages) delete it.second;
}

void infill::decodeDigits::beginJob()
{
  geom = art::ServiceHandle<geo::Geometry>()->provider();

  // Get active ROPs (not facing a wall)
  geo::ROP_id_iterator iRop, rBegin(geom, geo::ROP_id_iterator::begin_pos), 
    rEnd(geom, geo::ROP_id_iterator::end_pos);
  for (iRop = rBegin; iRop != rEnd; ++iRop) { // Iterate over ROPs in detector
    for (const geo::TPCID tpcId : geom->ROPtoTPCs(*iRop)) {
      const geo::TPCGeo tpc = geom->TPC(tpcId);
      const TGeoVolume* tpcVol = tpc.ActiveVolume();
      
      if (tpcVol->Capacity() > 1000000) { // At least one of the ROP's TPCIDs needs to be active
        activeRops.insert(*iRop);
        break;
      }
    }
  }
}

void infill::decodeDigits::endJob()
{
  std::cout << "Getting dead channels..." << std::endl;

  badChannels = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels();
  noisyChannels = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().NoisyChannels();

  std::merge(
    badChannels.begin(), badChannels.end(), noisyChannels.begin(), noisyChannels.end(), 
    std::inserter(deadChannels, deadChannels.begin())
  );

  // Print dead channels to terminal
  for(const readout::ROPID& rop : activeRops) {
    std::cout << "(TPCset " << rop.TPCset << ", ROP " << rop.ROP << "): ";
    const raw::ChannelID_t firstCh = geom->FirstChannelInROP(rop);
    for (const raw::ChannelID_t ch : deadChannels) {
      if (ch >= firstCh + geom->Nchannels(rop)) {
        break;
      }
      else if (ch >= firstCh) {
        std::cout << ch - firstCh << " ";
      }
    }
    std::cout << std::endl;
  }
}

DEFINE_ART_MODULE(infill::decodeDigits)
