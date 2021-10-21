////////////////////////////////////////////////////////////////////////
// Class:       CVNSparseROOT
// Plugin Type: analyzer (art v3_01_02)
// File:        CVNSparseROOT_module.cc
//
// Generated at Wed Apr 10 14:53:36 2019 by Jeremy Hewes using cetskelgen
// from cetlib version v3_05_01.
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

// CVN includes
#include "dune/CVN/func/SparsePixelMap.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// Boost includes
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming

namespace cvn {

  class CVNSparseROOT : public art::EDAnalyzer {
  public:
    explicit CVNSparseROOT(fhicl::ParameterSet const& p);

    CVNSparseROOT(CVNSparseROOT const&) = delete;
    CVNSparseROOT(CVNSparseROOT&&) = delete;
    CVNSparseROOT& operator=(CVNSparseROOT const&) = delete;
    CVNSparseROOT& operator=(CVNSparseROOT&&) = delete;

    void reconfigure(fhicl::ParameterSet const& p);

    void beginSubRun(art::SubRun const& sr) override;
    void endSubRun(art::SubRun const& sr) override;
    void analyze(art::Event const& e) override;

  private:

    std::string fMapModuleLabel; ///< Name of map producer module
    std::string fMapInstanceLabel; ///< Name of sparse map instance
    std::string fOutputName; ///< ROOT output filename
    std::string fTreeName; ///< ROOt tree name
    bool        fIncludeGroundTruth; ///< Whether to include per-pixel ground truth

    std::vector<std::vector<float>> fCoordinates; ///< Pixel coordinates
    std::vector<std::vector<float>> fFeatures; ///< Pixel features
    std::vector<std::vector<int>> fPixelPDG; ///< Pixel PDG truth
    std::vector<std::vector<int>> fPixelTrackID; ///< Pixel track ID
    std::vector<std::vector<float>> fPixelEnergies; ///< Pixel energy
    std::vector<std::vector<std::string>> fProcesses; // Physical process that creates the particle 

    std::vector<unsigned int> fEvent; ///< Event numbers
    unsigned int fView; ///< View numbers

    TFile* fFile; ///< Output ROOT file
    TTree* fTree; ///< ROOT tree for writing to file
    size_t fCacheSize; ///< Size of TTree cache

  };


  CVNSparseROOT::CVNSparseROOT(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {
    this->reconfigure(p);
  }

  void CVNSparseROOT::reconfigure(fhicl::ParameterSet const& p) {

    fMapModuleLabel     = p.get<std::string>("MapModuleLabel");
    fMapInstanceLabel   = p.get<std::string>("MapInstanceLabel");
    fOutputName         = p.get<std::string>("OutputName");
    fTreeName           = p.get<std::string>("TreeName");
    fIncludeGroundTruth = p.get<bool>("IncludeGroundTruth");
    fCacheSize          = p.get<size_t>("CacheSize");

  } // cvn::CVNSparseROOT::reconfigure

  void CVNSparseROOT::analyze(art::Event const& e) {

    // Get the sparse maps
    std::vector<art::Ptr<SparsePixelMap>> maps;
    art::InputTag itag1(fMapModuleLabel, fMapInstanceLabel);
    auto hMaps = e.getHandle<std::vector<SparsePixelMap>>(itag1);
    if (!hMaps) {
      throw art::Exception(art::errors::ProductNotFound)
        << "Could not find SparsePixelMap vector with module label "
        << fMapModuleLabel << " and instance label "
        << fMapInstanceLabel << "!" << std::endl;
    }
    art::fill_ptr_vector(maps, hMaps);

    if (maps.size() > 1) throw art::Exception(art::errors::LogicError)
      << "There shouldn't be more than one SparsePixelMap per producer per event,"
      << " but here there are " << maps.size() << "." << std::endl;

    if (maps.empty()) return;

    for (unsigned int it = 0; it < maps[0]->GetViews(); ++it) {
      fCoordinates = maps[0]->GetCoordinates(it);
      fFeatures = maps[0]->GetFeatures(it);
      if (fIncludeGroundTruth) {
        fPixelPDG      = maps[0]->GetPixelPDGs(it);
        fPixelTrackID  = maps[0]->GetPixelTrackIDs(it);
        fPixelEnergies = maps[0]->GetPixelEnergies(it);
        fProcesses     = maps[0]->GetProcesses(it);
      }
      fEvent = std::vector<unsigned int>({e.id().run(), e.id().subRun(), e.id().event()});
      fView = it;
      fTree->Fill();
    }

  } // cvn::CVNSparseROOT::analyze

  /// Beginning of a subrun, make a new file
  void CVNSparseROOT::beginSubRun(art::SubRun const& sr) {

    // Open ROOT file
    boost::uuids::random_generator generator;
    boost::uuids::uuid uuid = generator();
    std::ostringstream fileName;
    fileName << fOutputName << "_" << uuid << ".root";
    fFile = TFile::Open(fileName.str().c_str(), "recreate");

    fTree = new TTree(fTreeName.c_str(), fTreeName.c_str());
    fTree->SetCacheSize(fCacheSize);
    fTree->Branch("Coordinates", &fCoordinates);
    fTree->Branch("Features", &fFeatures);
    if (fIncludeGroundTruth) {
      fTree->Branch("PixelPDG", &fPixelPDG);
      fTree->Branch("PixelTrackID", &fPixelTrackID);
      fTree->Branch("PixelEnergy", &fPixelEnergies);
      fTree->Branch("Process", &fProcesses);
    }
    fTree->Branch("Event", &fEvent);
    fTree->Branch("View", &fView);

  } // function CVNSparseROOT::beginSubRun

  /// End of a subrun, write all events to a ROOT file
  void CVNSparseROOT::endSubRun(art::SubRun const& sr) {

    fFile->WriteTObject(fTree, fTreeName.c_str());
    delete fFile;

  } // cvn::CVNSparseROOT::endSubRun

  DEFINE_ART_MODULE(CVNSparseROOT)

} // namespace cvn

