////////////////////////////////////////////////////////////////////////
// Class:       GCNH5
// Plugin Type: analyzer (art v3_01_02)
// File:        GCNH5_module.cc
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

// Data product includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNFeatureUtils.h"

#include "hep_hpc/hdf5/make_ntuple.hpp"

// Boost includes
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming

using std::vector;
using std::string;
using std::round;
using std::abs;
using std::endl;
using std::setfill;
using std::setw;

using hep_hpc::hdf5::Column;
using hep_hpc::hdf5::make_scalar_column;

namespace cvn {

  class GCNH5 : public art::EDAnalyzer {
  public:
    explicit GCNH5(fhicl::ParameterSet const& p);
    ~GCNH5() noexcept {};

    GCNH5(GCNH5 const&) = delete;
    GCNH5(GCNH5&&) = delete;
    GCNH5& operator=(GCNH5 const&) = delete;
    GCNH5& operator=(GCNH5&&) = delete;

    void reconfigure(fhicl::ParameterSet const& p);

    void beginSubRun(art::SubRun const& sr) override;
    void endSubRun(art::SubRun const& sr) override;
    void analyze(art::Event const& e) override;

  private:

    string fGraphModuleLabel;   ///< Name of graph producer module
    string fGraphInstanceLabel; ///< Name of graph instance
    string fTruthLabel;         ///< Name of truth producer module
    string fOutputName;         ///< Output filename
    bool fSaveEventTruth;       ///< Whether to save event-level truth information
    bool fSaveParticleTruth;    ///< Whether to include particle truth information

    hep_hpc::hdf5::File fFile;  ///< Output HDF5 file
    hep_hpc::hdf5::Ntuple<Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<int, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<int, 1>,
                          Column<int, 1>>* fGraphNtuple; ///< graph ntuple

    hep_hpc::hdf5::Ntuple<Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>>* fEventNtuple; ///< Event ntuple

    hep_hpc::hdf5::Ntuple<Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<int, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<float, 1>,
                          Column<string, 1>,
                          Column<string, 1>>* fParticleNtuple; ///< Particle ntuple

  };

  GCNH5::GCNH5(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {
    this->reconfigure(p);
  }

  void GCNH5::reconfigure(fhicl::ParameterSet const& p) {

    fGraphModuleLabel   = p.get<string>("GraphModuleLabel");
    fGraphInstanceLabel = p.get<string>("GraphInstanceLabel");
    fTruthLabel         = p.get<string>("TruthLabel"),
    fOutputName         = p.get<string>("OutputName");
    fSaveEventTruth     = p.get<bool>("SaveEventTruth");
    fSaveParticleTruth  = p.get<bool>("SaveParticleTruth");

  } // cvn::GCNH5::reconfigure

  void GCNH5::analyze(art::Event const& e) {

    // Get the graphVector
    std::vector<art::Ptr<GCNGraph>> graphVector;
    art::InputTag itag1(fGraphModuleLabel, fGraphInstanceLabel);
    auto graphHandle = e.getHandle<std::vector<GCNGraph>>(itag1);
    if (!graphHandle) {
      throw art::Exception(art::errors::ProductNotFound)
        << "Could not find GCNGraph vector with module label "
        << fGraphModuleLabel << " and instance label "
        << fGraphInstanceLabel << "!" << endl;
    }
    art::fill_ptr_vector(graphVector, graphHandle);

    if (graphVector.size() > 1) throw art::Exception(art::errors::LogicError)
      << "There shouldn't be more than one GCNGraph per producer per event,"
      << " but here there are " << graphVector.size() << "." << endl;

    if (graphVector.empty()) return;

    int run = e.id().run();
    int subrun = e.id().subRun();
    int event = e.id().event();

    for (size_t itNode = 0; itNode < graphVector[0]->GetNumberOfNodes(); ++itNode) {
      GCNGraphNode node = graphVector[0]->GetNode(itNode);

      vector<float> pos = node.GetPosition();
      vector<float> feat = node.GetFeatures();
      vector<float> truth = node.GetGroundTruth();

      fGraphNtuple->insert(run, subrun, event, (int)round(abs(pos[0])),
        pos[1], pos[2], feat[0], feat[1], feat[4], feat[2], feat[3],
        (int)feat[5], truth[0]); 
    }

    // Event truth
    if (fSaveEventTruth) {

      // Get MC truth
      auto truthHandle = e.getHandle<vector<simb::MCTruth>>(fTruthLabel);
      if (!truthHandle || truthHandle->size() != 1) {
        throw art::Exception(art::errors::LogicError)
          << "Expected to find exactly one MC truth object!";
      }
      simb::MCNeutrino nutruth = truthHandle->at(0).GetNeutrino();

      // Fill variables
      fEventNtuple->insert(run, subrun, event,
        nutruth.CCNC() == simb::kCC,
        nutruth.Nu().E(),
        nutruth.Lepton().E(),
        nutruth.Nu().Momentum().Vect().Unit().X(),
        nutruth.Nu().Momentum().Vect().Unit().Y(),
        nutruth.Nu().Momentum().Vect().Unit().Z());

    } // if saving event truth

    // // Particle tree
    if (fSaveParticleTruth) {
      std::vector<ptruth> ptree = GCNFeatureUtils::GetParticleTree(graphVector[0].get());
      for (auto p : ptree) {
        fParticleNtuple->insert(run, subrun, event,
          std::get<0>(p), std::get<1>(p), std::get<2>(p),
          std::get<3>(p), std::get<4>(p), std::get<5>(p),
          std::get<6>(p), std::get<7>(p), std::get<8>(p),
          std::get<9>(p), std::get<10>(p), std::get<11>(p));
      }
    }

  } // cvn::GCNH5::analyze

  void GCNH5::beginSubRun(art::SubRun const& sr) {

    // Open HDF5 output
    boost::uuids::random_generator generator;
    boost::uuids::uuid uuid = generator();
    std::ostringstream fileName;
    fileName << fOutputName << "_r" << setfill('0') << setw(5) << sr.run()
      << "_r" << setfill('0') << setw(5) << sr.subRun() << "_" << uuid
      << ".h5";

    fFile = hep_hpc::hdf5::File(fileName.str(), H5F_ACC_TRUNC);

    fGraphNtuple = new hep_hpc::hdf5::Ntuple(
      make_ntuple({fFile, "graph_table", 1000},
      make_scalar_column<int>("run"),
      make_scalar_column<int>("subrun"),
      make_scalar_column<int>("event"),
      make_scalar_column<int>("plane"),
      make_scalar_column<float>("wire"),
      make_scalar_column<float>("time"),
      make_scalar_column<float>("integral"),
      make_scalar_column<float>("rms"),
      make_scalar_column<int>("rawplane"),
      make_scalar_column<float>("rawwire"),
      make_scalar_column<float>("rawtime"),
      make_scalar_column<int>("tpc"),
      make_scalar_column<int>("true_id")));

    if (fSaveEventTruth) {
      fEventNtuple = new hep_hpc::hdf5::Ntuple(
        make_ntuple({fFile, "event_table", 1000},
        make_scalar_column<int>("run"),
        make_scalar_column<int>("subrun"),
        make_scalar_column<int>("event"),
        make_scalar_column<int>("is_cc"),
        make_scalar_column<float>("nu_energy"),
        make_scalar_column<float>("lep_energy"),
        make_scalar_column<float>("nu_dir_x"),
        make_scalar_column<float>("nu_dir_y"),
        make_scalar_column<float>("nu_dir_z")));
    }

    if (fSaveParticleTruth) {
      fParticleNtuple = new hep_hpc::hdf5::Ntuple(
        make_ntuple({fFile, "particle_table", 1000},
        make_scalar_column<int>("run"),
        make_scalar_column<int>("subrun"),
        make_scalar_column<int>("event"),
        make_scalar_column<int>("id"),
        make_scalar_column<int>("type"),
        make_scalar_column<int>("parent_id"),
        make_scalar_column<float>("momentum"),
        make_scalar_column<float>("start_x"),
        make_scalar_column<float>("start_y"),
        make_scalar_column<float>("start_z"),
        make_scalar_column<float>("end_x"),
        make_scalar_column<float>("end_y"),
        make_scalar_column<float>("end_z"),
        make_scalar_column<string>("start_process"),
        make_scalar_column<string>("end_process")));
    }
  }

  void GCNH5::endSubRun(art::SubRun const& sr) {
    delete fGraphNtuple;
    if (fSaveEventTruth) delete fEventNtuple;
    if (fSaveParticleTruth) delete fParticleNtuple;
    fFile.close();
  }

  DEFINE_ART_MODULE(GCNH5)

} // namespace cvn

