////////////////////////////////////////////////////////////////////////
// Class:       InfillChannels
// Plugin Type: analyzer (art v3_05_01)
// File:        InfillChannels_module.cc
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
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include <TGeoVolume.h>

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <array>

#include <torch/script.h>
#include <torch/torch.h>

namespace Infill 
{
  class InfillChannels;
}

class Infill::InfillChannels : public art::EDProducer 
{
public:
  explicit InfillChannels(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  InfillChannels(InfillChannels const&) = delete;
  InfillChannels(InfillChannels&&) = delete;
  InfillChannels& operator=(InfillChannels const&) = delete;
  InfillChannels& operator=(InfillChannels&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  // Declare member data here.
  const geo::GeometryCore* fGeom;

  std::set<raw::ChannelID_t> fBadChannels; 
  std::set<raw::ChannelID_t> fNoisyChannels;
  std::set<raw::ChannelID_t> fDeadChannels;

  std::set<readout::ROPID> fActiveRops;
  
  const std::string fNetworkPath;
  const std::string fNetworkNameInduction;
  const std::string fNetworkNameCollection;
  torch::jit::script::Module fInductionModule;
  torch::jit::script::Module fCollectionModule;
  const std::string fInputLabel;
};

Infill::InfillChannels::InfillChannels(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fNetworkPath           (p.get<std::string> ("NetworkPath")),
    fNetworkNameInduction  (p.get<std::string> ("NetworkNameInduction")),
    fNetworkNameCollection (p.get<std::string> ("NetworkNameCollection")),
    fInputLabel            (p.get<std::string> ("InputLabel"))
{
  consumes<std::vector<raw::RawDigit>>(fInputLabel);

  produces<std::vector<raw::RawDigit>>();
}

void Infill::InfillChannels::produce(art::Event& e)
{
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
  // Networks expect a fixed image size
  if (detProp.NumberTimeSamples() > 6000) {
    std::cerr << "InfillChannels_module.cc: Networks cannot handle more than 6000 time ticks\n";
    std::abort(); 
  } 
      
  typedef std::array<short, 6000> vecAdc;
  std::map<raw::ChannelID_t, vecAdc> infilledAdcs;
  torch::Tensor maskedRopTensor;
  torch::Tensor infilledRopTensor; 

  auto digs = e.getHandle<std::vector<raw::RawDigit> >(fInputLabel);

  // Get infilled adc ROP by ROP
  for (const readout::ROPID& currentRop : fActiveRops) {
    maskedRopTensor = torch::zeros(
      {1, 1 ,6000, fGeom->Nchannels(currentRop)}, torch::dtype(torch::kFloat32).device(torch::kCPU).requires_grad(false)
    );
    auto maskedRopTensorAccess = maskedRopTensor.accessor<float, 4>();

    const raw::ChannelID_t firstCh = fGeom->FirstChannelInROP(currentRop);

    // Fill ROP image
    for (const raw::RawDigit& dig : *digs) {
      if (fDeadChannels.count(dig.Channel())) continue;

      readout::ROPID rop = fGeom->ChannelToROP(dig.Channel());
      if (rop != currentRop) continue;

      raw::RawDigit::ADCvector_t adcs(dig.Samples());
      raw::Uncompress(dig.ADCs(), adcs, dig.Compression());
      
      for (unsigned int tick = 0; tick < adcs.size(); ++tick) {
        const int adc = adcs[tick] ? int(adcs[tick]) - dig.GetPedestal() : 0;

        maskedRopTensorAccess[0][0][tick][dig.Channel() - firstCh] = adc;
      } 
    }

    // Do the Infill
    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(maskedRopTensor);
    if (fGeom->SignalType(currentRop) == geo::kInduction) {
      torch::NoGradGuard no_grad_guard; 
      infilledRopTensor = fInductionModule.forward(inputs).toTensor().detach();
    }
    else if (fGeom->SignalType(currentRop) == geo::kCollection) {
      torch::NoGradGuard no_grad_guard; 
      infilledRopTensor = fCollectionModule.forward(inputs).toTensor().detach();
    }

    // Store infilled ADC of dead channels
    auto infilledRopTensorAccess = infilledRopTensor.accessor<float, 4>();
    for (const raw::ChannelID_t ch : fDeadChannels) {
      if (fGeom->ChannelToROP(ch) == currentRop) {
        for (unsigned int tick = 0; tick < detProp.NumberTimeSamples(); ++tick) {
          infilledAdcs[ch][tick] = (short)std::round(infilledRopTensorAccess[0][0][tick][ch - firstCh]);
        }
      }
    } // Could break early if ch > last ch in currentRop to save a bit of looping?
  }

  // Encode infilled ADC into RawDigit and put back onto event 
  auto infilledDigs = std::make_unique<std::vector<raw::RawDigit>>();
  *infilledDigs = *digs;
  for (raw::RawDigit& dig : *infilledDigs) {
    if (infilledAdcs.count(dig.Channel())) {
      raw::RawDigit::ADCvector_t infilledAdc(
        infilledAdcs[dig.Channel()].begin(), (infilledAdcs[dig.Channel()].begin() + detProp.NumberTimeSamples())
      );

      // Get new pedestal
      auto infilledAdcMin = std::min_element(infilledAdc.begin(), infilledAdc.end());
      short ped = *infilledAdcMin < 0 ? std::abs(*infilledAdcMin) + 1 : 0;
      for (short& adc : infilledAdc) adc += ped;

      raw::Compress(infilledAdc, dig.Compression()); // need to consider compression parameters
      dig = raw::RawDigit(dig.Channel(), dig.Samples(), infilledAdc, dig.Compression());
      dig.SetPedestal(ped);
    }
  }
  e.put(std::move(infilledDigs)); 
}

void Infill::InfillChannels::beginJob()
{
  fGeom = art::ServiceHandle<geo::Geometry>()->provider();

  // Dead channels = bad channels + noisy channels
  fBadChannels = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().BadChannels();
  fNoisyChannels = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider().NoisyChannels();

  std::merge(
    fBadChannels.begin(), fBadChannels.end(), fNoisyChannels.begin(), fNoisyChannels.end(), 
    std::inserter(fDeadChannels, fDeadChannels.begin())
  );

  // Get active ROPs (not facing a wall and has dead channels)
  for (auto const& ropID : fGeom->Iterate<readout::ROPID>()) { // Iterate over ROPs in the detector
    bool hasDeadCh = false;
    for (raw::ChannelID_t ch : fDeadChannels) {
      if (fGeom->ChannelToROP(ch) == ropID) {
        hasDeadCh = true;
        break;
      }
    }
    if (!hasDeadCh) continue; // Don't need to infill ROPs without dead channels

    for (const geo::TPCID tpcId : fGeom->ROPtoTPCs(ropID)) {
      const geo::TPCGeo tpc = fGeom->TPC(tpcId);
      const TGeoVolume* tpcVol = tpc.ActiveVolume();
      
      if (tpcVol->Capacity() > 1000000) { // At least one of the ROP's TPCIDs needs to be active
        // Networks expect a fixed image size
        if(fGeom->SignalType(ropID) == geo::kInduction && fGeom->Nchannels(ropID) > 800) {
	  std::cerr << "InfillChannels_module.cc: Induction view network cannot handle more then 800 channels\n";
	  std::abort();
        }
        if(fGeom->SignalType(ropID) == geo::kCollection && fGeom->Nchannels(ropID) > 480) {
	  std::cerr << "InfillChannels_module.cc: Collection view network cannot handle more then 400 channels\n";
	  std::abort();
        }

        fActiveRops.insert(ropID);
        break;
      }
    }
  }
  
  // Check dead channels resemble the dead channels used for training
  raw::ChannelID_t chGap = 1;
  for (const raw::ChannelID_t ch : fDeadChannels) {
    if (fDeadChannels.count(ch + 1)) {
      ++chGap;
      continue;
    }
    if (fGeom->ChannelToROP(ch - chGap) == fGeom->ChannelToROP(ch + 1)) {
      if (fGeom->SignalType(ch) == geo::kCollection && chGap > 3) {
        std::cerr << "There are dead channel gap larger than what was seen in training --- ";
        std::cerr << "**Consider retraining collection plane infill network**" << std::endl;
      }
      else if (fGeom->SignalType(ch) == geo::kInduction && chGap > 2) {
        std::cerr << "There are dead channel gap larger than what was seen in training --- ";
        std::cerr << "**Consider retraining induction plane infill network**" << std::endl;
      }
    }
    chGap = 1;
  }

  // Load torchscripts
  std::cout << "Loading modules..." << std::endl;
  const char* networkPath = std::getenv(fNetworkPath.c_str());
  if (networkPath == nullptr) {
    std::cerr << "InfillChannels_module.cc: Environment variable " << fNetworkPath << " was not found";
    std::abort();
  }
  const std::string networkLocInduction = std::string(networkPath) + "/" + fNetworkNameInduction;
  const std::string networkLocCollection = std::string(networkPath) + "/" + fNetworkNameCollection;

  try {
    fInductionModule = torch::jit::load(networkLocInduction);
    std::cout << "Induction module loaded from " << networkLocInduction <<std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what();
  }
  try {
    fCollectionModule = torch::jit::load(networkLocCollection);
    std::cout << "Collection module loaded from " << networkLocCollection << std::endl;
  }
  catch (const c10::Error& err) {
    std::cerr << "error loading the model\n";
    std::cerr << err.what();
  }
}

void Infill::InfillChannels::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Infill::InfillChannels)
