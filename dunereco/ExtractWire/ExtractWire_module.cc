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
#include "hep_hpc/hdf5/make_ntuple.hpp"
#include "lardataobj/RecoBase/Wire.h"
#include <array>

namespace wire {
  class ExtractWire;
}


class wire::ExtractWire : public art::EDAnalyzer {

using wire_range_t = std::pair<size_t, size_t>;
using output_ranges_t = std::vector<std::pair<wire_range_t, size_t>>;

using ntuple_t
  = hep_hpc::hdf5::Ntuple<
    hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    hep_hpc::hdf5::Column<size_t, 1>,    // Wire Number and Tick
    hep_hpc::hdf5::Column<float, 1>  // Value of waveform
  >;

using event_ntuple_t
  = hep_hpc::hdf5::Ntuple<
    hep_hpc::hdf5::Column<int, 1>,    // event id (run, subrun, event)
    hep_hpc::hdf5::Column<size_t, 1>    // Wire Number and Tick
  >;

  
public:
  explicit ExtractWire(fhicl::ParameterSet const& p);
  ~ExtractWire() noexcept {};
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExtractWire(ExtractWire const&) = delete;
  ExtractWire(ExtractWire&&) = delete;
  ExtractWire& operator=(ExtractWire const&) = delete;
  ExtractWire& operator=(ExtractWire&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginSubRun(art::SubRun const& sr) override;
  void endSubRun(art::SubRun const& sr) override;

private:

  void WireLoop(
    const std::array<int, 3> & evtID,
    const art::ValidHandle<std::vector<recob::Wire>> & wires,
    std::vector<ntuple_t*> & ntuple_vec,
    std::vector<event_ntuple_t*> & event_ntuple_vec
  );

  void MakeNTuples(
    std::vector<ntuple_t*> & wires_ntuple,
    std::vector<event_ntuple_t*> & event_ntuple,
    std::string base_name);

  // Declare member data here.
  //Main --> Truth (i.e. 'perfect' APA)
  //Alt --> Input (i.e. 'bad' APA)
  art::InputTag fWireLabelForInput, fWireLabelForTruth;
  std::string fOutputName;
  output_ranges_t fOutputRanges;
  hep_hpc::hdf5::File fOutputFile;

  std::vector<ntuple_t*> fNTupleForInput, fNTupleForTruth; //
  std::vector<event_ntuple_t*> fEventNTupleForInput, fEventNTupleForTruth; // for event info
};


wire::ExtractWire::ExtractWire(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
    fWireLabelForInput(p.get<art::InputTag>("WireLabelForInput")),
    fWireLabelForTruth(p.get<art::InputTag>("WireLabelForTruth")),
    fOutputName(p.get<std::string>("OutputName", "extracted_wires")),
    fOutputRanges(p.get<output_ranges_t>("OutputRanges")) {

}

void wire::ExtractWire::MakeNTuples(
  std::vector<ntuple_t*> & wires_ntuple,
  std::vector<event_ntuple_t*> & event_ntuple,
  std::string base_name) {
  for (auto & range : fOutputRanges) {
    std::ostringstream ntuple_name;
    ntuple_name << base_name << "_" << range.second;
    wires_ntuple.push_back(
      new hep_hpc::hdf5::Ntuple(
        hep_hpc::hdf5::make_ntuple(
          {fOutputFile, ntuple_name.str(), 1000},
          hep_hpc::hdf5::make_column<int>("event_id", 3),
          hep_hpc::hdf5::make_column<size_t>("chan_tick", 2),
          hep_hpc::hdf5::make_scalar_column<float>("value", 1)
        )
      )
    );

    std::ostringstream event_ntuple_name;
    event_ntuple_name << base_name << "_event_" << range.second;
    event_ntuple.push_back(
      new hep_hpc::hdf5::Ntuple(
        hep_hpc::hdf5::make_ntuple(
          {fOutputFile, event_ntuple_name.str(), 1000},
          hep_hpc::hdf5::make_column<int>("event_id", 3),
          hep_hpc::hdf5::make_column<size_t>("roi_count", 1)
        )
      )
    );
  }
}

void wire::ExtractWire::beginSubRun(art::SubRun const& sr) {

  // struct timeval now;
  // gettimeofday(&now, NULL);

  // Open HDF5 output
  std::ostringstream fileName;
  fileName << fOutputName << "_r" << std::setfill('0') << std::setw(5) << sr.run()
    << "_s" << std::setfill('0') << std::setw(5) << sr.subRun() << ".h5";
  fOutputFile = hep_hpc::hdf5::File(fileName.str(), H5F_ACC_TRUNC);

  MakeNTuples(fNTupleForInput, fEventNTupleForInput, "input_wire_info");
  MakeNTuples(fNTupleForTruth, fEventNTupleForTruth, "truth_wire_info");
  
}

void wire::ExtractWire::endSubRun(art::SubRun const& sr) {
  fOutputFile.close();
  for (auto * ptr : fNTupleForInput) delete ptr;
  for (auto * ptr : fNTupleForTruth) delete ptr;
  for (auto * ptr : fEventNTupleForInput) delete ptr;
  for (auto * ptr : fEventNTupleForTruth) delete ptr;
}
void wire::ExtractWire::WireLoop(
    const std::array<int, 3> & evtID,
    const art::ValidHandle<std::vector<recob::Wire>> & wires,
    std::vector<ntuple_t*> & ntuple_vec,
    std::vector<event_ntuple_t*> & event_ntuple_vec) {


  std::map<size_t, size_t> counts;

  for (const auto & wire : (*wires)) {
    bool found_out = false;
    size_t out_index = 0;
    size_t start_wire = 0;
    for (const auto & [range, index] : fOutputRanges) {
      if (range.first <= wire.Channel() && wire.Channel() < range.second) {
        found_out = true;
        out_index = index;
        start_wire = range.first;
        break;
      }
    }
    if (!found_out) continue;

       
    // std::cout << wire.Channel() << " " << wire.View() << std::endl;
    // std::cout << "\t" << wire.SignalROI().size() << std::endl;
    // std::cout << out_index << std::endl;

    // std::cout << "Ranges:" << std::endl;
    for (const auto & range : wire.SignalROI().get_ranges()) {
      // std::cout << range.begin_index() << std::endl;
      size_t range_index = range.begin_index();
      for (const auto & v : range.data()) {
        // std::cout << "\t" << range_index << " " << v << std::endl;
        std::array<size_t, 2> chan_tick {(wire.Channel()-start_wire), range_index};
        ntuple_vec.at(out_index)->insert(
          evtID.data(),
          chan_tick.data(),
          v
        );
        ++range_index;

        if (counts.find(out_index) == counts.end()) {
          counts[out_index] = 0;
        }
        ++counts[out_index];
      }
    }
  }

  for (const auto & [range, index] : fOutputRanges) {
    size_t count = counts[index];
    event_ntuple_vec.at(index)->insert(
      evtID.data(),
      count
    );
  }

}

void wire::ExtractWire::analyze(art::Event const& e) {

  int run = e.id().run();
  int subrun = e.id().subRun();
  int event = e.id().event();

  std::array<int, 3> evtID { run, subrun, event };

  auto wires_for_input = e.getValidHandle<std::vector<recob::Wire>>(fWireLabelForInput);
  auto wires_for_truth = e.getValidHandle<std::vector<recob::Wire>>(fWireLabelForTruth);

  WireLoop(evtID, wires_for_input, fNTupleForInput, fEventNTupleForInput);
  WireLoop(evtID, wires_for_truth, fNTupleForTruth, fEventNTupleForTruth);

}

DEFINE_ART_MODULE(wire::ExtractWire)
