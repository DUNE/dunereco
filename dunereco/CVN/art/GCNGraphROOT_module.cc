////////////////////////////////////////////////////////////////////////
// Class:       GCNGraphROOT
// Plugin Type: analyzer (art v3_01_02)
// File:        GCNGraphROOT_module.cc
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
#include "dune/CVN/func/GCNParticleFlow.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// Boost includes
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming

namespace cvn {

  class GCNGraphROOT : public art::EDAnalyzer {
  public:
    explicit GCNGraphROOT(fhicl::ParameterSet const& p);

    GCNGraphROOT(GCNGraphROOT const&) = delete;
    GCNGraphROOT(GCNGraphROOT&&) = delete;
    GCNGraphROOT& operator=(GCNGraphROOT const&) = delete;
    GCNGraphROOT& operator=(GCNGraphROOT&&) = delete;

    void reconfigure(fhicl::ParameterSet const& p);

    void beginSubRun(art::SubRun const& sr) override;
    void endSubRun(art::SubRun const& sr) override;
    void analyze(art::Event const& e) override;

  private:

    std::string fGraphModuleLabel;   ///< Name of graph producer module
    std::string fGraphInstanceLabel; ///< Name of graph instance
    std::string fTruthLabel;         ///< Name of truth producer module
    std::string fOutputName;         ///< ROOT output filename
    std::string fTreeName;           ///< ROOT tree name
    bool fSaveEventTruth;            ///< Whether to save event-level truth information
    bool fSaveParticleFlow;          ///< Whether to include particle flow information

    std::vector<std::vector<float>> fPosition;    ///< Node positions
    std::vector<std::vector<float>> fFeatures;    ///< Node features
    std::vector<std::vector<float>> fGroundTruth; ///< Node ground truth

    std::vector<unsigned int> fEvent; ///< Event numbers

    bool fIsCC;       ///< Whether the neutrino interaction is charged current
    float fNuEnergy;  ///< True neutrino energy
    float fLepEnergy; ///< True lepton energy
    float fNuDirX;    ///< X component of true neutrino direction
    float fNuDirY;    ///< Y component of true neutrino direction
    float fNuDirZ;    ///< Z component of true neutrino direction

    std::map<unsigned int, unsigned int> fParticleFlow; ///< Particle flow map

    TFile* fFile;      ///< Output ROOT file
    TTree* fTree;      ///< ROOT tree for writing to file

  };


  GCNGraphROOT::GCNGraphROOT(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {
    this->reconfigure(p);
  }

  void GCNGraphROOT::reconfigure(fhicl::ParameterSet const& p) {

    fGraphModuleLabel   = p.get<std::string>("GraphModuleLabel");
    fGraphInstanceLabel = p.get<std::string>("GraphInstanceLabel");
    fTruthLabel         = p.get<std::string>("TruthLabel"),
    fOutputName         = p.get<std::string>("OutputName");
    fTreeName           = p.get<std::string>("TreeName");
    fSaveEventTruth     = p.get<bool>("SaveEventTruth");
    fSaveParticleFlow   = p.get<bool>("SaveParticleFlow");

  } // cvn::GCNGraphROOT::reconfigure

  void GCNGraphROOT::analyze(art::Event const& e) {

    // Get the graphVector
    art::Handle<std::vector<GCNGraph>> graphHandle;
    std::vector<art::Ptr<GCNGraph>> graphVector;
    if (!e.getByLabel(fGraphModuleLabel, fGraphInstanceLabel, graphHandle)) {
      throw art::Exception(art::errors::ProductNotFound)
        << "Could not find GCNGraph vector with module label "
        << fGraphModuleLabel << " and instance label "
        << fGraphInstanceLabel << "!" << std::endl;
    }
    art::fill_ptr_vector(graphVector, graphHandle);

    if (graphVector.size() > 1) throw art::Exception(art::errors::LogicError)
      << "There shouldn't be more than one GCNGraph per producer per event,"
      << " but here there are " << graphVector.size() << "." << std::endl;

    if (graphVector.empty()) return;

    // Empty vectors
    fPosition.clear();
    fFeatures.clear();
    fGroundTruth.clear();

    // Loop over nodes and refill them
    for (size_t itNode = 0; itNode < graphVector[0]->GetNumberOfNodes(); ++itNode) {
      GCNGraphNode node = graphVector[0]->GetNode(itNode);
      fPosition.push_back(node.GetPosition());
      fFeatures.push_back(node.GetFeatures());
      fGroundTruth.push_back(node.GetGroundTruth());
    }

    fEvent = std::vector<unsigned int>({e.id().run(), e.id().subRun(), e.id().event()});

    // Event truth
    if (fSaveEventTruth) {

      // Get MC truth
      art::Handle<std::vector<simb::MCTruth>> truthHandle;
      e.getByLabel(fTruthLabel, truthHandle);
      if (!truthHandle.isValid() || truthHandle->size() != 1) {
        throw art::Exception(art::errors::LogicError)
          << "Expected to find exactly one MC truth object!";
      }
      simb::MCNeutrino nutruth = truthHandle->at(0).GetNeutrino();

      // Fill variables
      fIsCC = (nutruth.CCNC() == simb::kCC);
      fNuEnergy = nutruth.Nu().E();
      fLepEnergy = nutruth.Lepton().E();
      fNuDirX = nutruth.Nu().Momentum().Vect().Unit().X();
      fNuDirX = nutruth.Nu().Momentum().Vect().Unit().Y();
      fNuDirX = nutruth.Nu().Momentum().Vect().Unit().Z();

    } // if saving event truth

    // Particle flow
    if (fSaveParticleFlow) {

      // Get particle flow
      art::Handle<std::vector<cvn::GCNParticleFlow>> pfHandle;
      e.getByLabel(fGraphModuleLabel, fGraphInstanceLabel, pfHandle);
      if (!pfHandle.isValid() || pfHandle->size() != 1) {
        throw art::Exception(art::errors::LogicError)
          << "Expected exactly one graph particle flow object.";
      }
      fParticleFlow = pfHandle->at(0).GetMap();

    } // if saving particle flow

    fTree->Fill();

  } // cvn::GCNGraphROOT::analyze

  /// Beginning of a subrun, make a new file
  void GCNGraphROOT::beginSubRun(art::SubRun const& sr) {

    // Open ROOT file
    boost::uuids::random_generator generator;
    boost::uuids::uuid uuid = generator();
    std::ostringstream fileName;
    fileName << fOutputName << "_" << uuid << ".root";
    fFile = TFile::Open(fileName.str().c_str(), "recreate");

    fTree = new TTree(fTreeName.c_str(), fTreeName.c_str());
    fTree->Branch("Position", &fPosition);
    fTree->Branch("Features", &fFeatures);
    fTree->Branch("GroundTruth", &fGroundTruth);
    fTree->Branch("Event", &fEvent);

    if (fSaveEventTruth) {
      fTree->Branch("IsCC", &fIsCC);
      fTree->Branch("NuEnergy", &fNuEnergy);
      fTree->Branch("LepEnergy", &fLepEnergy);
      fTree->Branch("NuDirX", &fNuDirX);
      fTree->Branch("NuDirY", &fNuDirY);
      fTree->Branch("NuDirZ", &fNuDirZ);
    } // If saving event truth

    if (fSaveParticleFlow) {
      fTree->Branch("ParticleFlow", &fParticleFlow);
    } // If saving particle flow

  } // function GCNGraphROOT::beginSubRun

  /// End of a subrun, write all events to a ROOT file
  void GCNGraphROOT::endSubRun(art::SubRun const& sr) {

    fFile->WriteTObject(fTree, fTreeName.c_str());
    delete fFile;

  } // cvn::GCNGraphROOT::endSubRun

  DEFINE_ART_MODULE(GCNGraphROOT)

} // namespace cvn

