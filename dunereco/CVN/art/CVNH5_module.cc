////////////////////////////////////////////////////////////////////////
// \file    CVNH5_module.cc
// \brief   Analyzer module for creating CVN HDF5 file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <experimental/filesystem>

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/Exception.h"

// Data products
#include "nusimdata/SimulationBase/MCTruth.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/TrainingData.h"
#include "dune/CVN/func/InteractionType.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/CVNImageUtils.h"

// Compression
#include <dune/CVN/highfive/H5DataSet.hpp>
#include <dune/CVN/highfive/H5DataSpace.hpp>
#include <dune/CVN/highfive/H5File.hpp>

namespace cvn {

  class CVNH5 : public art::EDAnalyzer {
  public:

    explicit CVNH5(fhicl::ParameterSet const& pset);
    ~CVNH5();

    void beginSubRun(const art::SubRun& sr) override;
    void endSubRun(const art::SubRun& sr) override;
    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);

  private:

    std::string fOutputDir;
    std::string fPixelMapInput;
    bool fSetLog;
    std::vector<bool> fReverseViews;
    unsigned int fTopologyHitsCut;
    unsigned int fNPixWire, fNPixTime;

    std::string fGenieGenModuleLabel;
    std::string fEnergyNueLabel;
    std::string fEnergyNumuLabel;
    std::string fEnergyNutauLabel;

    std::string out_dir;

    std::vector<ImageVector> fImageVectors;
    std::vector<int> fPDG, fNProton, fNPion, fNPi0,
      fNNeutron, fTopType, fTopTypeAlt;
    std::vector<float> fNuEnergy, fLepEnergy,
      fNueEnergy, fNumuEnergy, fNutauEnergy;

    void write_files(TrainingData td, unsigned int n);

  };

  //......................................................................
  CVNH5::CVNH5(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNH5::~CVNH5()
  {}

  //......................................................................
  void CVNH5::reconfigure(const fhicl::ParameterSet& pset)
  {
    fOutputDir = pset.get<std::string>("OutputDir", "");
    fPixelMapInput = pset.get<std::string>("PixelMapInput");
    fSetLog = pset.get<bool>("SetLog");
    fReverseViews = pset.get<std::vector<bool>>("ReverseViews");
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");
    fNPixWire = pset.get<unsigned int>("NPixWire");
    fNPixTime = pset.get<unsigned int>("NPixTime");

    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
    fEnergyNueLabel = pset.get<std::string>("EnergyNueLabel");
    fEnergyNumuLabel = pset.get<std::string>("EnergyNumuLabel");
    fEnergyNutauLabel = pset.get<std::string>("EnergyNutauLabel");
  }

  //......................................................................
  void CVNH5::beginSubRun(const art::SubRun& sr)
  {
    fImageVectors.clear();

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
  }

  //......................................................................
  void CVNH5::endSubRun(const art::SubRun& sr) {

    try {
      std::ostringstream fileName;
      fileName << "CVN_" << sr.run() << "_" << sr.subRun() << ".h5";

      // Create a new file using the default property lists.
      HighFive::File f(fileName.str(),
        HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

      // Create the dataset
      f.createDataSet("imagevector", fImageVectors);

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
      // catch and print any HDF5 error
      std::cerr << err.what() << std::endl;
    }
  }

  //......................................................................
  void CVNH5::analyze(const art::Event& evt)
  {
    // Get the pixel maps
    art::Handle<std::vector<cvn::PixelMap>> h_pixelmaps;
    std::vector<art::Ptr<cvn::PixelMap>> pixelmaps;
    if (evt.getByLabel(fPixelMapInput, fPixelMapInput, h_pixelmaps))
      art::fill_ptr_vector(pixelmaps, h_pixelmaps);

    // If no pixel maps, quit
    if (pixelmaps.size() == 0) return;

    // MC information
    art::Handle<std::vector<simb::MCTruth>> h_mctruth;
    std::vector<art::Ptr<simb::MCTruth>> mctruth_list;
    if (evt.getByLabel(fGenieGenModuleLabel, h_mctruth))
      art::fill_ptr_vector(mctruth_list, h_mctruth);

    art::Ptr<simb::MCTruth> mctruth = mctruth_list[0];
    simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();

    AssignLabels labels;
    labels.GetTopology(mctruth, fTopologyHitsCut);

    // True lepton and neutrino energies
    fNuEnergy.push_back(true_neutrino.Nu().E());
    fLepEnergy.push_back(true_neutrino.Lepton().E());

    // Put a containment cut here

    float reco_nue_energy = 0;
    float reco_numu_energy = 0;
    float reco_nutau_energy = 0;

    // Get nue info
    if (fEnergyNueLabel != "") {
      art::Handle<dune::EnergyRecoOutput> h_ereco;
      evt.getByLabel(fEnergyNueLabel, h_ereco);
      reco_nue_energy = h_ereco->fNuLorentzVector.E();
    }

    // Get numu info
    if (fEnergyNueLabel != "") {
      art::Handle<dune::EnergyRecoOutput> h_ereco;
      evt.getByLabel(fEnergyNumuLabel, h_ereco);
      reco_numu_energy = h_ereco->fNuLorentzVector.E();
    }

    // Get nutau info
    if (fEnergyNutauLabel != "") {
      art::Handle<dune::EnergyRecoOutput> h_ereco;
      evt.getByLabel(fEnergyNutauLabel, h_ereco);
      reco_nutau_energy = h_ereco->fNuLorentzVector.E();
    }

    fNueEnergy.push_back(reco_nue_energy);
    fNumuEnergy.push_back(reco_numu_energy);
    fNutauEnergy.push_back(reco_nutau_energy);

    fPDG.push_back(labels.GetPDG());
    fNProton.push_back(labels.GetNProtons());
    fNPion.push_back(labels.GetNPions());
    fNPi0.push_back(labels.GetNPizeros());
    fNNeutron.push_back(labels.GetNNeutrons());
    fTopType.push_back(labels.GetTopologyType());
    fTopTypeAlt.push_back(labels.GetTopologyTypeAlt());

    CVNImageUtils image_utils(fNPixWire, fNPixTime, 3);
    ImageVector vec(3, std::vector<std::vector<unsigned char>>(fNPixWire,
     std::vector<unsigned char>(fNPixTime, 0)));

    image_utils.SetLogScale(fSetLog);
    image_utils.SetViewReversal(fReverseViews);
    image_utils.SetPixelMapSize(pixelmaps[0]->fNWire, pixelmaps[0]->fNTdc);
    image_utils.ConvertPixelMapToImageVector(*pixelmaps[0], vec);
    fImageVectors.push_back(vec);
  }

DEFINE_ART_MODULE(cvn::CVNH5)
} // namespace cvn
