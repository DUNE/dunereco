////////////////////////////////////////////////////////////////////////
// \file    CVNZlibMaker_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
//           - wrote the zlib code used in this module
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>

#include "boost/filesystem.hpp"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
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
#include "zlib.h"

namespace fs = boost::filesystem;

namespace cvn {

  class CVNZlibMaker : public art::EDAnalyzer {
  public:

    explicit CVNZlibMaker(fhicl::ParameterSet const& pset);
    ~CVNZlibMaker();

    void beginJob() override;
    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);

  private:

    std::string fOutputDir;
    std::string fPixelMapInput;
    bool fSetLog;
    std::vector<bool> fReverseViews;
    unsigned int fTopologyHitsCut;

    std::string fGenieGenModuleLabel;
    std::string fEnergyNueLabel;
    std::string fEnergyNumuLabel;
    std::string fEnergyNutauLabel;

    unsigned int fPlaneLimit;
    unsigned int fTDCLimit;

    std::string out_dir;

    void write_files(TrainingData td, unsigned int n);

  };

  //......................................................................
  CVNZlibMaker::CVNZlibMaker(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNZlibMaker::~CVNZlibMaker()
  {  }

  //......................................................................
  void CVNZlibMaker::reconfigure(const fhicl::ParameterSet& pset)
  {
    fOutputDir = pset.get<std::string>("OutputDir", "");
    fPixelMapInput = pset.get<std::string>("PixelMapInput");
    fSetLog = pset.get<bool>("SetLog");
    fReverseViews = pset.get<std::vector<bool>>("ReverseViews");
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");

    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
    fEnergyNueLabel = pset.get<std::string>("EnergyNueLabel");
    fEnergyNumuLabel = pset.get<std::string>("EnergyNumuLabel");
    fEnergyNutauLabel = pset.get<std::string>("EnergyNutauLabel");

    fPlaneLimit = pset.get<unsigned int>("PlaneLimit");
    fTDCLimit = pset.get<unsigned int>("TDCLimit");
  }

  //......................................................................
  void CVNZlibMaker::beginJob()
  {
    // Set the output directory.
    // First we look for CONDOR_DIR_INPUT and set it to the grid location.
    // If it isn't there, we look for a FHICL parameter.
    // Otherwise it just goes to the current directory.
/*
    char const* grid_dir = getenv("CONDOR_DIR_INPUT");
    if (grid_dir != NULL) {
      char const* tmp_grid_dir = getenv("TMP");
      if (tmp_grid_dir == NULL)
        throw art::Exception(art::errors::NotFound)
          << "Could not find environment variable \"TMP\" which "
          << "is usually set for condor_lar environment." << std::endl;
      out_dir = std::string(tmp_grid_dir) + "/out";
    }
*/
    if (fOutputDir != "")
      out_dir = fOutputDir;

    else
      out_dir = ".";

    // Throw an error if the specified output directory doesn't exist
    if (!fs::exists(out_dir))
      throw art::Exception(art::errors::FileOpenError)
        << "Output directory " << out_dir << " does not exist!" << std::endl;

    // std::cout << "Writing files to output directory " << out_dir << std::endl;
  }

  //......................................................................
  void CVNZlibMaker::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    std::vector<art::Ptr<cvn::PixelMap>> pixelmaps;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto h_pixelmaps = evt.getHandle<std::vector<cvn::PixelMap>>(itag1);
    if (h_pixelmaps)
      art::fill_ptr_vector(pixelmaps, h_pixelmaps);

    // If no pixel maps, quit
    if (pixelmaps.size() == 0) return;

    InteractionType interaction = kOther;

    // MC information
    std::vector<art::Ptr<simb::MCTruth>> mctruth_list;
    auto h_mctruth = evt.getHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
    if (h_mctruth)
      art::fill_ptr_vector(mctruth_list, h_mctruth);

    art::Ptr<simb::MCTruth> mctruth = mctruth_list[0];
    simb::MCNeutrino true_neutrino = mctruth->GetNeutrino();

    AssignLabels labels;

    interaction = labels.GetInteractionType(true_neutrino);
    labels.GetTopology(mctruth, fTopologyHitsCut);

    // True lepton and neutrino energies
    float nu_energy = true_neutrino.Nu().E();
    float lep_energy = true_neutrino.Lepton().E();

    // Put a containment cut here

    float reco_nue_energy = 0;
    float reco_numu_energy = 0;
    float reco_nutau_energy = 0;

    // Get nue info
    if (fEnergyNueLabel != "") {
      auto h_ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyNueLabel);
      reco_nue_energy = h_ereco->fNuLorentzVector.E();
    }

    // Get numu info
    if (fEnergyNueLabel != "") {
      auto h_ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyNumuLabel);
      reco_numu_energy = h_ereco->fNuLorentzVector.E();
    }

    // Get nutau info
    if (fEnergyNutauLabel != "") {
      auto h_ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyNutauLabel);
      reco_nutau_energy = h_ereco->fNuLorentzVector.E();
    }

    // Hard-coding event weight for now
    // Should probably fix this at some point
    int event_weight = 1;

    TrainingData train(interaction, nu_energy, lep_energy,
      reco_nue_energy, reco_numu_energy, reco_nutau_energy,
      event_weight, *pixelmaps[0]);

    int pdg        = labels.GetPDG();
    int n_proton   = labels.GetNProtons();
    int n_pion     = labels.GetNPions();
    int n_pi0      = labels.GetNPizeros();
    int n_neutron  = labels.GetNNeutrons();
    int toptype    = labels.GetTopologyType();
    int toptypealt = labels.GetTopologyTypeAlt();

    train.SetTopologyInformation(pdg, n_proton, n_pion,
      n_pi0, n_neutron, toptype, toptypealt);

    this->write_files(train, evt.event());
  }

  //......................................................................
  void CVNZlibMaker::write_files(TrainingData td, unsigned int n)
  {
    // cropped from 2880 x 500 to 500 x 500 here 
    std::vector<unsigned char> pixel_array(3 * fPlaneLimit * fTDCLimit);
   
    CVNImageUtils image_utils(fPlaneLimit, fTDCLimit, 3);
    image_utils.SetPixelMapSize(td.fPMap.NWire(), td.fPMap.NTdc());
    image_utils.SetLogScale(fSetLog);
    image_utils.SetViewReversal(fReverseViews);
    image_utils.ConvertPixelMapToPixelArray(td.fPMap, pixel_array);

    ulong src_len = 3 * fPlaneLimit * fTDCLimit; // pixelArray length
    ulong dest_len = compressBound(src_len);     // calculate size of the compressed data               
    char* ostream = (char *) malloc(dest_len);  // allocate memory for the compressed data

    int res = compress((Bytef *) ostream, &dest_len, (Bytef *) &pixel_array[0], src_len);

    // Buffer error

    if (res == Z_BUF_ERROR)
      std::cout << "Buffer too small!" << std::endl;

    // Memory error
    else if (res ==  Z_MEM_ERROR)
      std::cout << "Not enough memory for compression!" << std::endl;

    // Compression ok 
    else {

      // Create output files 
      std::string image_file_name = out_dir + "/event_" + std::to_string(n) + ".gz";
      std::string info_file_name = out_dir + "/event_" + std::to_string(n) + ".info";

      std::ofstream image_file (image_file_name, std::ofstream::binary);
      std::ofstream info_file  (info_file_name);

      if(image_file.is_open() && info_file.is_open()) {

        // Write compressed data to file

        image_file.write(ostream, dest_len);

        image_file.close(); // close file

        // Write records to file

        // Category

        info_file << td.fInt << std::endl;

        // Energy

        info_file << td.fNuEnergy << std::endl;
        info_file << td.fLepEnergy << std::endl;
        info_file << td.fRecoNueEnergy << std::endl;
        info_file << td.fRecoNumuEnergy << std::endl;
        info_file << td.fRecoNutauEnergy << std::endl;
        info_file << td.fEventWeight << std::endl;

        // Topology

        info_file << td.fNuPDG << std::endl;
        info_file << td.fNProton << std::endl;
        info_file << td.fNPion << std::endl;
        info_file << td.fNPizero << std::endl;         
        info_file << td.fNNeutron << std::endl;

        info_file << td.fTopologyType << std::endl;
        info_file << td.fTopologyTypeAlt;        

        info_file.close(); // close file
      }
      else {

        if (image_file.is_open())
          image_file.close();
        else 
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << image_file_name << "!" << std::endl;

        if (info_file.is_open())
          info_file.close();
        else
          throw art::Exception(art::errors::FileOpenError)
            << "Unable to open file " << info_file_name << "!" << std::endl;
      }
    }
    
    free(ostream);  // free allocated memory

  } // cvn::CVNZlibMaker::write_files

DEFINE_ART_MODULE(cvn::CVNZlibMaker)
} // namespace cvn
