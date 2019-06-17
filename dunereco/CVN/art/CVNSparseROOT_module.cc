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

// highfive includes
#include <dune/CVN/highfive/H5DataSet.hpp>
#include <dune/CVN/highfive/H5DataSpace.hpp>
#include <dune/CVN/highfive/H5File.hpp>

// ROOT includes
#include "TFile.h"
#include "TTree.h"

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

    std::vector<std::vector<unsigned int>> fCoordinatesU; ///< Pixel coordinates for U plane
    std::vector<std::vector<unsigned int>> fCoordinatesV; ///< Pixel coordinates for V plane
    std::vector<std::vector<unsigned int>> fCoordinatesY; ///< Pixel coordinates for Y plane
    std::vector<float> fValuesU; ///< Pixel values for U plane
    std::vector<float> fValuesV; ///< Pixel values for V plane
    std::vector<float> fValuesY; ///< Pixel values for Y plane
    std::vector<unsigned int> fEvent; ///< Event numbers

    TFile* fFile; ///< Output ROOT file
    TTree* fTree; ///< ROOT tree for writing to file

  };


  CVNSparseROOT::CVNSparseROOT(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {
    this->reconfigure(p);
  }

  void CVNSparseROOT::reconfigure(fhicl::ParameterSet const& p) {

    fMapModuleLabel =   p.get<std::string>("MapModuleLabel");
    fMapInstanceLabel = p.get<std::string>("MapInstanceLabel");
    fOutputName =       p.get<std::string>("OutputName");
    fTreeName =         p.get<std::string>("TreeName");

  } // cvn::CVNSparseROOT::reconfigure

  void CVNSparseROOT::analyze(art::Event const& e) {

    // Get the sparse maps
    art::Handle<std::vector<SparsePixelMap>> hMaps;
    std::vector<art::Ptr<SparsePixelMap>> maps;
    if (!e.getByLabel(fMapModuleLabel, fMapInstanceLabel, hMaps)) {
      throw art::Exception(art::errors::ProductNotFound)
        << "Could not find SparsePixelMap vector with module label "
        << fMapModuleLabel << " and instance label "
        << fMapInstanceLabel << "!" << std::endl;
    }
    art::fill_ptr_vector(maps, hMaps);

    if (maps.size() > 1) throw art::Exception(art::errors::LogicError)
      << "There shouldn't be more than one SparsePixelMap per producer per event,"
      << " but here there are " << maps.size() << "." << std::endl;

    if (maps.size() > 0) {
      fCoordinatesU = maps[0]->GetCoordinates(0);
      fCoordinatesV = maps[0]->GetCoordinates(1);
      fCoordinatesY = maps[0]->GetCoordinates(2);
      fValuesU = maps[0]->GetValues(0);
      fValuesV = maps[0]->GetValues(1);
      fValuesY = maps[0]->GetValues(2);
      fEvent = std::vector<unsigned int>({e.id().run(), e.id().subRun(), e.id().event()});
      fTree->Fill();
    }

  } // cvn::CVNSparseROOT::analyze

  /// Beginning of a subrun, make a new file
  void CVNSparseROOT::beginSubRun(art::SubRun const& sr) {

    // Open ROOT file
    std::ostringstream fileName;
    fileName << fOutputName << "_" << sr.id().run() << "_" << sr.id().subRun() << ".root";
    fFile = TFile::Open(fileName.str().c_str(), "recreate");

    fTree = new TTree(fTreeName.c_str(), fTreeName.c_str());
    fTree->Branch("CoordinatesU", &fCoordinatesU);
    fTree->Branch("CoordinatesV", &fCoordinatesV);
    fTree->Branch("CoordinatesY", &fCoordinatesY);
    fTree->Branch("ValuesU", &fValuesU);
    fTree->Branch("ValuesV", &fValuesV);
    fTree->Branch("ValuesY", &fValuesY);
    fTree->Branch("Event", &fEvent);

  } // function CVNSparseROOT::beginSubRun

  /// End of a subrun, write all events to a ROOT file
  void CVNSparseROOT::endSubRun(art::SubRun const& sr) {

    fFile->WriteTObject(fTree, fTreeName.c_str());
    delete fFile;

  } // cvn::CVNSparseROOT::endSubRun

  DEFINE_ART_MODULE(CVNSparseROOT)

} // namespace cvn
