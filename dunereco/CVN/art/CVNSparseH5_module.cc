////////////////////////////////////////////////////////////////////////
// Class:       CVNSparseH5
// Plugin Type: analyzer (art v3_01_02)
// File:        CVNSparseH5_module.cc
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

namespace cvn {

  class CVNSparseH5 : public art::EDAnalyzer {
  public:
    explicit CVNSparseH5(fhicl::ParameterSet const& p);

    CVNSparseH5(CVNSparseH5 const&) = delete;
    CVNSparseH5(CVNSparseH5&&) = delete;
    CVNSparseH5& operator=(CVNSparseH5 const&) = delete;
    CVNSparseH5& operator=(CVNSparseH5&&) = delete;

    void reconfigure(fhicl::ParameterSet const& p);

    void beginSubRun(art::SubRun const& sr) override;
    void endSubRun(art::SubRun const& sr) override;
    void analyze(art::Event const& e) override;

  private:

    std::string fMapModuleLabel; ///< Name of map producer module
    std::string fMapInstanceLabel; ///< Name of sparse map instance
    std::string fOutputName; ///< H5 output filename

    std::vector<std::vector<unsigned int>> fCoordinates; ///< Pixel coordinates
    std::vector<float> fValues; ///< Pixel values
    std::vector<unsigned int> fNPixels; ///< Pixel map boundaries
    std::vector<std::vector<unsigned int>> fEvents; ///< Event numbers

  };


  CVNSparseH5::CVNSparseH5(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {

    this->reconfigure(p);

  }

  void CVNSparseH5::reconfigure(fhicl::ParameterSet const& p) {

    fMapModuleLabel =   p.get<std::string>("MapModuleLabel");
    fMapInstanceLabel = p.get<std::string>("MapInstanceLabel");
    fOutputName =       p.get<std::string>("OutputName");

  } // cvn::CVNSparseH5::reconfigure

  void CVNSparseH5::analyze(art::Event const& e) {

    // Get the sparse maps
    std::cout << "Getting sparse pixel maps from art event!" << std::endl;
    art::Handle<std::vector<SparsePixelMap>> hMaps;
    std::vector<art::Ptr<SparsePixelMap>> maps;
    if (!e.getByLabel(fMapModuleLabel, fMapInstanceLabel, hMaps)) {
      throw art::Exception(art::errors::ProductNotFound)
        << "Could not find SparsePixelMap vector with module label "
        << fMapModuleLabel << " and instance label "
        << fMapInstanceLabel << "!" << std::endl;
    }
    art::fill_ptr_vector(maps, hMaps);

    std::cout << "Size of pixel map vector from art event is " << maps.size() << std::endl;

    std::cout << "Filling information vectors." << std::endl;
    for (auto map : maps) {
      fEvents.emplace_back(std::vector<unsigned int>({e.id().run(), e.id().subRun(), e.id().event()}));
      std::vector<unsigned int> npixels = map->GetNPixels();
      fNPixels.insert(fNPixels.end(), npixels.begin(), npixels.end());
      std::vector<std::vector<unsigned int>> coords = map->GetCoordinatesFlat();
      fCoordinates.insert(fCoordinates.end(), coords.begin(), coords.end());
      std::vector<float> values = map->GetValuesFlat();
      fValues.insert(fValues.end(), values.begin(), values.end());
    }

  } // cvn::CVNSparseH5::analyze

  /// Beginning of a subrun, clear everything out
  void CVNSparseH5::beginSubRun(art::SubRun const& sr) {

    std::cout << "Clearing event vectors." << std::endl;
    fEvents.clear();
    fNPixels.clear();
    fCoordinates.clear();
    fValues.clear();

  } // function CVNSparseH5::beginSubRun

  /// End of a subrun, write all graphs to a H5 file
  void CVNSparseH5::endSubRun(art::SubRun const& sr) {

    // Ignore empty subruns
    if (fEvents.empty()) return;

    try {
      // Open H5 file
      std::ostringstream file_name;
      file_name << fOutputName << "_" << sr.id().run() << "_" << sr.id().subRun() << ".h5";
      HighFive::File f(file_name.str(),
        HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

      // Write images
      f.createDataSet("events", fEvents);
      f.createDataSet("npixels", fNPixels);
      f.createDataSet("coordinates", fCoordinates);
      f.createDataSet("values", fValues);

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }

  } // cvn::CVNSparseH5::endSubRun

  DEFINE_ART_MODULE(CVNSparseH5)

} // namespace cvn
