////////////////////////////////////////////////////////////////////////
// Class:       ExtractWire
// Plugin Type: analyzer (Unknown Unknown)
// File:        ExtractWire_module.cc
//
// Generated at Tue Apr  1 08:40:55 2025 by Jacob Calcutt using cetskelgen
// from cetlib version 3.18.02.
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

#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/Ntuple.hpp"
#include "lardataobj/RecoBase/Wire.h"

namespace wire {
  class ExtractWire;
}


class wire::ExtractWire : public art::EDAnalyzer {
public:
  explicit ExtractWire(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExtractWire(ExtractWire const&) = delete;
  ExtractWire(ExtractWire&&) = delete;
  ExtractWire& operator=(ExtractWire const&) = delete;
  ExtractWire& operator=(ExtractWire&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  //Main --> Truth (i.e. 'perfect' APA)
  //Alt --> Input (i.e. 'bad' APA)
  art::InputTag fWireLabelMain, fWireLabelAlt;
};


wire::ExtractWire::ExtractWire(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fWireLabelMain(p.get<art::InputTag>("WireLabelMain")),
    fWireLabelAlt(p.get<art::InputTag>("WireLabelAlt")) {}

void wire::ExtractWire::analyze(art::Event const& e) {
  auto main_wires = e.getValidHandle<std::vector<recob::Wire>>(fWireLabelMain);
  // auto alt_wires = e.getValidHandle<std::vector<recob::Wire>>(fWireLabelAlt);


  for (const auto & wire : (*main_wires)) {
    std::cout << wire.Channel() << std::endl;
  }
}

DEFINE_ART_MODULE(wire::ExtractWire)
