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

// Data products
#include "nusimdata/SimulationBase/MCTruth.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"
#include "dune/CVN/func/SparsePixelMap.h"
#include "dune/CVN/func/AssignLabels.h"

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

    std::string fMapModuleLabel;      ///< Name of map producer module
    std::string fMapInstanceLabel;    ///< Name of sparse map instance
    std::string fOutputName;          ///< H5 output filename
    std::string fGenieGenModuleLabel; ///< MC truth producer module
    std::string fEnergyNueLabel;      ///< Nue hypothesis reco energy
    std::string fEnergyNumuLabel;     ///< Numu hypothesis reco energy
    std::string fEnergyNutauLabel;    ///< Calorimetric reco energy

    std::vector<std::vector<unsigned int>> fCoordinates; ///< Pixel coordinates
    std::vector<float> fValues;                          ///< Pixel values
    std::vector<unsigned int> fNPixels;                  ///< Pixel map boundaries
    std::vector<std::vector<unsigned int>> fEvents;      ///< Event numbers
    std::vector<unsigned int> fViews;                    ///< View of each pixel map

    unsigned int fTopologyHitsCut;
    std::vector<int> fPDG, fNProton, fNPion, fNPi0,
      fNNeutron, fTopType, fTopTypeAlt;
    std::vector<float> fNuEnergy, fLepEnergy,
      fNueEnergy, fNumuEnergy, fNutauEnergy;

  };


  CVNSparseH5::CVNSparseH5(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {

    this->reconfigure(p);

  }

  void CVNSparseH5::reconfigure(fhicl::ParameterSet const& p) {

    fMapModuleLabel =   p.get<std::string>("MapModuleLabel");
    fMapInstanceLabel = p.get<std::string>("MapInstanceLabel");
    fOutputName =       p.get<std::string>("OutputName");

    fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");
    fEnergyNueLabel      = p.get<std::string>("EnergyNueLabel");
    fEnergyNumuLabel     = p.get<std::string>("EnergyNumuLabel");
    fEnergyNutauLabel    = p.get<std::string>("EnergyNutauLabel");

    fTopologyHitsCut = p.get<unsigned int>("TopologyHitsCut");

  } // cvn::CVNSparseH5::reconfigure

  void CVNSparseH5::analyze(art::Event const& e) {

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

    // MC information
    art::Handle<std::vector<simb::MCTruth>> hMCTruth;
    std::vector<art::Ptr<simb::MCTruth>> MCTruthList;
    if (e.getByLabel(fGenieGenModuleLabel, hMCTruth))
      art::fill_ptr_vector(MCTruthList, hMCTruth);

    art::Ptr<simb::MCTruth> MCTruth = MCTruthList[0];
    simb::MCNeutrino trueNeutrino = MCTruth->GetNeutrino();

    AssignLabels labels;
    labels.GetTopology(MCTruth, fTopologyHitsCut);

    float recoNueEnergy = 0;
    float recoNumuEnergy = 0;
    float recoNutauEnergy = 0;

    // Get nue info
    if (fEnergyNueLabel != "") {
      art::Handle<dune::EnergyRecoOutput> hEReco;
      e.getByLabel(fEnergyNueLabel, hEReco);
      recoNueEnergy = hEReco->fNuLorentzVector.E();
    }

    // Get numu info
    if (fEnergyNueLabel != "") {
      art::Handle<dune::EnergyRecoOutput> hEReco;
      e.getByLabel(fEnergyNumuLabel, hEReco);
      recoNumuEnergy = hEReco->fNuLorentzVector.E();
    }

    // Get nutau info
    if (fEnergyNutauLabel != "") {
      art::Handle<dune::EnergyRecoOutput> hEReco;
      e.getByLabel(fEnergyNutauLabel, hEReco);
      recoNutauEnergy = hEReco->fNuLorentzVector.E();
    }


    for (auto map : maps) {
      for (size_t view = 0; view < map->GetViews(); ++view) {

        // Sparse pixel map information
        fEvents.emplace_back(std::vector<unsigned int>({e.id().run(), e.id().subRun(), e.id().event()}));
        fViews.emplace_back(view);
        fNPixels.emplace_back(map->GetNPixels(view));
        std::vector<std::vector<unsigned int>> coords = map->GetCoordinates(view);
        fCoordinates.insert(fCoordinates.end(), coords.begin(), coords.end());
        std::vector<float> values = map->GetValues(view);
        fValues.insert(fValues.end(), values.begin(), values.end());

        // True lepton and neutrino energies
        fNuEnergy.push_back(trueNeutrino.Nu().E());
        fLepEnergy.push_back(trueNeutrino.Lepton().E());

        // Reco energies
        fNueEnergy.push_back(recoNueEnergy);
        fNumuEnergy.push_back(recoNumuEnergy);
        fNutauEnergy.push_back(recoNutauEnergy);

        // Other truth information
        fPDG.push_back(labels.GetPDG());
        fNProton.push_back(labels.GetNProtons());
        fNPion.push_back(labels.GetNPions());
        fNPi0.push_back(labels.GetNPizeros());
        fNNeutron.push_back(labels.GetNNeutrons());
        fTopType.push_back(labels.GetTopologyType());
        fTopTypeAlt.push_back(labels.GetTopologyTypeAlt());

      } // for view
    } // for map

  } // cvn::CVNSparseH5::analyze

  /// Beginning of a subrun, clear everything out
  void CVNSparseH5::beginSubRun(art::SubRun const& sr) {

    fEvents.clear();
    fViews.clear();
    fNPixels.clear();
    fCoordinates.clear();
    fValues.clear();

    fPDG.clear();
    fNProton.clear();
    fNPion.clear();
    fNPi0.clear();
    fNNeutron.clear();
    fTopType.clear();
    fTopTypeAlt.clear();

    fNuEnergy.clear();
    fLepEnergy.clear();
    fNueEnergy.clear();
    fNumuEnergy.clear();
    fNutauEnergy.clear();

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
      f.createDataSet("views", fViews);
      f.createDataSet("npixels", fNPixels);
      f.createDataSet("coordinates", fCoordinates);
      f.createDataSet("values", fValues);

      // Write topology information
      f.createDataSet("pdg", fPDG);
      f.createDataSet("n_protons", fNProton);
      f.createDataSet("n_pions", fNPion);
      f.createDataSet("n_pi0s", fNPi0);
      f.createDataSet("n_neutrons", fNNeutron);
      f.createDataSet("toptype", fTopType);
      f.createDataSet("toptypealt", fTopTypeAlt);

      f.createDataSet("nu_energy", fNuEnergy);
      f.createDataSet("lep_energy", fLepEnergy);
      f.createDataSet("nue_energy", fNueEnergy);
      f.createDataSet("numu_energy", fNumuEnergy);
      f.createDataSet("nutau_energy", fNutauEnergy);

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }

  } // cvn::CVNSparseH5::endSubRun

  DEFINE_ART_MODULE(CVNSparseH5)

} // namespace cvn
