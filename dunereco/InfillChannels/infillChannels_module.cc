////////////////////////////////////////////////////////////////////////
// Class:       infillChannels
// Plugin Type: analyzer (art v3_05_01)
// File:        infillChannels_module.cc
//
// Alex Wilkinson - 08/03/21
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/raw.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include <TGeoVolume.h>

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <array>

#include <torch/script.h>
#include <torch/torch.h>

namespace infill 
{
  class infillChannels;
}

class infill::infillChannels : public art::EDProducer 
{
public:
  explicit infillChannels(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  infillChannels(infillChannels const&) = delete;
  infillChannels(infillChannels&&) = delete;
  infillChannels& operator=(infillChannels const&) = delete;
  infillChannels& operator=(infillChannels&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.
  const geo::GeometryCore* geom;

  std::set<raw::ChannelID_t> badChannels; 
  std::set<raw::ChannelID_t> noisyChannels;
  std::set<raw::ChannelID_t> deadChannels;

  std::set<readout::ROPID> activeRops;

  std::string networkLocInduction;
  std::string networkLocCollection;
  torch::jit::script::Module inductionModule;
  torch::jit::script::Module collectionModule;
};

infill::infillChannels::infillChannels(fhicl::ParameterSet const& p)
  : EDProducer{p},
    networkLocInduction  (p.get<std::string> ("networkLocInduction")),
    networkLocCollection (p.get<std::string> ("networkLocCollection"))
{
  consumes<std::vector<raw::RawDigit>>("daq");

  produces<std::vector<raw::RawDigit>>("infilled");
  // produces<std::vector<raw::RawDigit>>("masked");
}

void infill::infillChannels::produce(art::Event& e)
{
  typedef std::array<short, 6000> vecAdc; // A way to get number of ticks? (avoid hardcoded 6000)
  std::map<raw::ChannelID_t, vecAdc> infilledAdcs;
  torch::Tensor maskedRopTensor;
  torch::Tensor infilledRopTensor; 

  art::Handle<std::vector<raw::RawDigit>> digs;
  e.getByLabel("daq", digs);

  // Get infilled adc ROP by ROP
  for (const readout::ROPID& currentRop : activeRops) {
    maskedRopTensor = torch::zeros(
      {1, 1 ,6000, geom->Nchannels(currentRop)}, torch::dtype(torch::kFloat32).device(torch::kCPU).requires_grad(false)
    );
    auto maskedRopTensorAccess = maskedRopTensor.accessor<float, 4>();

    const raw::ChannelID_t firstCh = geom->FirstChannelInROP(currentRop);

    // Fill ROP image
    for (const raw::RawDigit& dig : *digs) {
      if (deadChannels.count(dig.Channel())) continue;

      readout::ROPID rop = geom->ChannelToROP(dig.Channel());
      if (rop != currentRop) continue;

      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression()); // can put pedestal in here

      for (unsigned int tick = 0; tick < adcs.size(); ++tick) {
        const short adc = adcs[tick] ? short(adcs[tick]) - dig.GetPedestal() : 0;

        maskedRopTensorAccess[0][0][tick][dig.Channel() - firstCh] = adc;
      } 
    }

    // Do the Infill
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(maskedRopTensor);
    if (geom->SignalType(currentRop) == geo::kInduction) {
      torch::NoGradGuard no_grad_guard; 
      infilledRopTensor = inductionModule.forward(inputs).toTensor().detach();
    }
    else if (geom->SignalType(currentRop) == geo::kCollection) {
      torch::NoGradGuard no_grad_guard; 
      infilledRopTensor = collectionModule.forward(inputs).toTensor().detach();
    }

    // Store infilled ADC of dead channels
    auto infilledRopTensorAccess = infilledRopTensor.accessor<float, 4>();
    for (const raw::ChannelID_t ch : deadChannels) {
      if (ch > firstCh + geom->Nchannels(currentRop)) {// Using some ChannelToROP would be better (think this algo is in another module as well)
        break;
      }
      else if (ch > firstCh) {
        for (unsigned int tick = 0; tick < infilledRopTensor.sizes()[2]; ++tick) {
          infilledAdcs[ch][tick] = (short)std::round(infilledRopTensorAccess[0][0][tick][ch - firstCh]);
        }
      }
    }
  }

  // Encode infilled ADC into RawDigit and put back onto event 
  auto infilledDigs = std::make_unique<std::vector<raw::RawDigit>>();
  *infilledDigs = *digs;
  for (raw::RawDigit& dig : *infilledDigs) {
    if (infilledAdcs.count(dig.Channel())) {
      raw::RawDigit::ADCvector_t infilledAdc(
        infilledAdcs[dig.Channel()].begin(), infilledAdcs[dig.Channel()].end()
      );
      raw::Compress(infilledAdc, dig.Compression()); // need to consider compression parameters
      dig = raw::RawDigit(dig.Channel(), dig.Samples(), infilledAdc, dig.Compression());
    }
  }
  e.put(std::move(infilledDigs), "infilled"); 
}

void infill::infillChannels::beginJob()
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

  // Dead channels = bad channels + noisy channels
  badChannels = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels();
  noisyChannels = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().NoisyChannels();

  std::merge(
    badChannels.begin(), badChannels.end(), noisyChannels.begin(), noisyChannels.end(), 
    std::inserter(deadChannels, deadChannels.begin())
  );

  // Load torchscripts
  std::cout << "Loading modules..." << std::endl;
  try {
    inductionModule = torch::jit::load(networkLocInduction);
    std::cout << "Induction module loaded." << std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what();
  }
  try {
    collectionModule = torch::jit::load(networkLocCollection);
    std::cout << "Collection module loaded." << std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what();
  }
}

void infill::infillChannels::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(infill::infillChannels)
