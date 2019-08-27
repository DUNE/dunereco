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

    bool fIncludeMCTruth;             ///< Whether to write MC truth information to HDF5
    std::string fGenieGenModuleLabel; ///< MC truth producer module
    std::vector<int> fPDG, fNProton, fNPion, fNPi0,
      fNNeutron, fTopType, fTopTypeAlt; ///< True topology information
    std::vector<float> fNuEnergy, fLepEnergy; ///< True energy information

    bool fIncludeRecoEnergy;          ///< Whether to write reco energy information to HDF5
    std::string fEnergyNueLabel;      ///< Nue hypothesis reco energy
    std::string fEnergyNumuLabel;     ///< Numu hypothesis reco energy
    std::string fEnergyNutauLabel;    ///< Calorimetric reco energy
    std::vector<float> fNueEnergy, fNumuEnergy, fNutauEnergy; ///< Reconstructed energies

    bool fIncludePixelTruth;          ///< Whether to write per-pixel ground truth to HDF5
    std::vector<int> fPixelTrackID;   ///< Pixel true track ID
    std::vector<int> fPixelPDG;       ///< Pixel true particle PDG

    std::vector<std::vector<unsigned int>> fCoordinates; ///< Pixel coordinates
    std::vector<float> fValues;                          ///< Pixel values
    std::vector<unsigned int> fNPixels;                  ///< Pixel map boundaries
    std::vector<std::vector<unsigned int>> fEvents;      ///< Event numbers
    std::vector<unsigned int> fViews;                    ///< View of each pixel map

    unsigned int fTopologyHitsCut;    ///< Minimum number of hits to create image

  };


  CVNSparseH5::CVNSparseH5(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {

    this->reconfigure(p);

  }

  void CVNSparseH5::reconfigure(fhicl::ParameterSet const& p) {

    fMapModuleLabel =   p.get<std::string>("MapModuleLabel");
    fMapInstanceLabel = p.get<std::string>("MapInstanceLabel");
    fOutputName =       p.get<std::string>("OutputName");

    fIncludeMCTruth      = p.get<bool>("IncludeMCTruth");
    fGenieGenModuleLabel = p.get<std::string>("GenieGenModuleLabel");

    fIncludeRecoEnergy   = p.get<bool>("IncludeRecoEnergy");
    fEnergyNueLabel      = p.get<std::string>("EnergyNueLabel");
    fEnergyNumuLabel     = p.get<std::string>("EnergyNumuLabel");
    fEnergyNutauLabel    = p.get<std::string>("EnergyNutauLabel");

    fIncludePixelTruth   = p.get<bool>("IncludePixelTruth");

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
    int pdg(0), nProton(0), nPion(0), nPi0(0), nNeutron(0),
      topType(0), topTypeAlt(0);
    float trueEnergy(0.), lepEnergy(0.);
    if (fIncludeMCTruth) {
      art::Handle<std::vector<simb::MCTruth>> hMCTruth;
      std::vector<art::Ptr<simb::MCTruth>> MCTruthList;
      if (e.getByLabel(fGenieGenModuleLabel, hMCTruth))
        art::fill_ptr_vector(MCTruthList, hMCTruth);
      art::Ptr<simb::MCTruth> MCTruth = MCTruthList[0];
      simb::MCNeutrino trueNeutrino = MCTruth->GetNeutrino();
      AssignLabels labels;
      labels.GetTopology(MCTruth, fTopologyHitsCut);
      trueEnergy = trueNeutrino.Nu().E();
      lepEnergy  = trueNeutrino.Lepton().E();
      pdg        = labels.GetPDG();
      nProton    = labels.GetNProtons();
      nPion      = labels.GetNPions();
      nPi0       = labels.GetNPizeros();
      nNeutron   = labels.GetNNeutrons();
      topType    = labels.GetTopologyType();
      topTypeAlt = labels.GetTopologyTypeAlt();
    } // if including MC truth

    float recoNueEnergy(0.), recoNumuEnergy(0.), recoNutauEnergy(0.);
    if (fIncludeRecoEnergy) {
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
    } // if including reco energy

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

        // Per pixel ground truth
        if (fIncludePixelTruth) {
          std::vector<int> trackID = map->GetPixelTrackID(view);
          fPixelTrackID.insert(fPixelTrackID.end(), trackID.begin(), trackID.end());
          std::vector<int> pdg = map->GetPixelPDG(view);
          fPixelPDG.insert(fPixelPDG.end(), pdg.begin(), pdg.end());
        }

        // MC truth information
        if (fIncludeMCTruth) {
          fNuEnergy.push_back(trueEnergy);
          fLepEnergy.push_back(lepEnergy);
          fPDG.push_back(pdg);
          fNProton.push_back(nProton);
          fNPion.push_back(nPion);
          fNPi0.push_back(nPi0);
          fNNeutron.push_back(nNeutron);
          fTopType.push_back(topType);
          fTopTypeAlt.push_back(topTypeAlt);
        } // if including MC truth

        // Reco energy information
        if (fIncludeRecoEnergy) {
          fNueEnergy.push_back(recoNueEnergy);
          fNumuEnergy.push_back(recoNumuEnergy);
          fNutauEnergy.push_back(recoNutauEnergy);
        }

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

    fPixelTrackID.clear();
    fPixelPDG.clear();

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

      // MC truth information
      if (fIncludeMCTruth) {
        f.createDataSet("nu_energy", fNuEnergy);
        f.createDataSet("lep_energy", fLepEnergy);
        f.createDataSet("pdg", fPDG);
        f.createDataSet("n_protons", fNProton);
        f.createDataSet("n_pions", fNPion);
        f.createDataSet("n_pi0s", fNPi0);
        f.createDataSet("n_neutrons", fNNeutron);
        f.createDataSet("toptype", fTopType);
        f.createDataSet("toptypealt", fTopTypeAlt);
      }

      // Reconstructed information
      if (fIncludeRecoEnergy) {
        f.createDataSet("nue_energy", fNueEnergy);
        f.createDataSet("numu_energy", fNumuEnergy);
        f.createDataSet("nutau_energy", fNutauEnergy);
      }

      // Pixel ground truth
      if (fIncludePixelTruth) {
        f.createDataSet("pixel_track_id", fPixelTrackID);
        f.createDataSet("pixel_pdg", fPixelPDG);
      }

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }

  } // cvn::CVNSparseH5::endSubRun

  DEFINE_ART_MODULE(CVNSparseH5)

} // namespace cvn
