////////////////////////////////////////////////////////////////////////
// \file    CVNZlibMakerProtoDUNE_module.cc
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
#include "canvas/Utilities/Exception.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
// Data products
#include "nusimdata/SimulationBase/MCParticle.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/CVNImageUtils.h"

// Compression
#include "zlib.h"

namespace fs = boost::filesystem;

namespace cvn {

  class CVNZlibMakerProtoDUNE : public art::EDAnalyzer {
  public:

    explicit CVNZlibMakerProtoDUNE(fhicl::ParameterSet const& pset);
    ~CVNZlibMakerProtoDUNE();

    void beginJob() override;
    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);

    struct PrimaryTrainingInfo {

      PrimaryTrainingInfo(TVector3 vtx, unsigned short inter, int pdg, float en){
        vertex = vtx;
        interaction = inter;
        pdgCode = pdg;
        energy = en;
      };

      TVector3 vertex;
      unsigned short interaction;
      int pdgCode;
      float energy;

    };

  private:

    std::string fOutputDir;
    std::string fPixelMapInput;
    bool fSetLog;
    std::vector<bool> fReverseViews;

    std::string fLArG4ModuleLabel;

    std::string out_dir;

    void write_files(const PrimaryTrainingInfo &primary, const art::Ptr<cvn::PixelMap> pm, unsigned int n) const;

  };

  //......................................................................
  CVNZlibMakerProtoDUNE::CVNZlibMakerProtoDUNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  CVNZlibMakerProtoDUNE::~CVNZlibMakerProtoDUNE()
  {  }

  //......................................................................
  void CVNZlibMakerProtoDUNE::reconfigure(const fhicl::ParameterSet& pset)
  {
    fOutputDir = pset.get<std::string>("OutputDir", "");
    fPixelMapInput = pset.get<std::string>("PixelMapInput");
    fSetLog = pset.get<bool>("SetLog");
    fReverseViews = pset.get<std::vector<bool>>("ReverseViews");
    fLArG4ModuleLabel = pset.get<std::string>("LArG4ModuleLabel");
  }

  //......................................................................
  void CVNZlibMakerProtoDUNE::beginJob()
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
  void CVNZlibMakerProtoDUNE::analyze(const art::Event& evt)
  {

    // Get the pixel maps
    std::vector<art::Ptr<cvn::PixelMap>> pixelmaps;
    art::InputTag itag1(fPixelMapInput, fPixelMapInput);
    auto h_pixelmaps = evt.getHandle<std::vector<cvn::PixelMap>>(itag1);
    if (h_pixelmaps)
      art::fill_ptr_vector(pixelmaps, h_pixelmaps);

    // If no pixel maps, quit
    if (pixelmaps.size() == 0) return;

    // MC information
    art::ServiceHandle< cheat::ParticleInventoryService > piService;

    // The truth information we need is the interaction vertex and interaction type
    AssignLabels labels;
    unsigned short beamParticleInteraction;

    // We only have one beam particle called "primary", so find it. It should be the
    // first particle, but best not to assume
    float beamParticleEnergy = 0;
    TVector3 beamParticleVtx; // This is the interaction vertex
    int beamParticlePDG = 0;

    bool gotPrimary = false;
    for(auto const m : piService->ParticleList()){
      const simb::MCParticle* particle = m.second;
      if(particle->Process().compare("primary")==0){
        beamParticleEnergy = particle->E();
        beamParticleVtx.SetXYZ(particle->EndX(),particle->EndY(),particle->EndZ());
        beamParticlePDG = particle->PdgCode();
        beamParticleInteraction = labels.GetProtoDUNEBeamInteractionType(*particle);
        gotPrimary = true;
        break;
      }
    }

    if(!gotPrimary) return;

    const PrimaryTrainingInfo beamPrimary(beamParticleVtx,beamParticleInteraction,beamParticlePDG,beamParticleEnergy);

    this->write_files(beamPrimary, pixelmaps.at(0), evt.event());
  }

  //......................................................................
  void CVNZlibMakerProtoDUNE::write_files(const PrimaryTrainingInfo &primary, const art::Ptr<cvn::PixelMap> pm, unsigned int n) const
  {
    CVNImageUtils image_utils;
    std::vector<unsigned char> pixel_array(3 * pm->NWire() * pm->NTdc());

    image_utils.DisableRegionSelection();
    image_utils.SetLogScale(fSetLog);
    image_utils.SetViewReversal(fReverseViews);
    image_utils.ConvertPixelMapToPixelArray(*(pm.get()),pixel_array);

    ulong src_len = 3 * pm->NWire() * pm->NTdc(); // pixelArray length
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
      std::string image_file_name = out_dir + "/cvn_event_" + std::to_string(n) + ".gz";
      std::string info_file_name = out_dir + "/cvn_event_" + std::to_string(n) + ".info";

      std::ofstream image_file (image_file_name, std::ofstream::binary);
      std::ofstream info_file  (info_file_name);

      if(image_file.is_open() && info_file.is_open()) {

        // Write compressed data to file

        image_file.write(ostream, dest_len);

        image_file.close(); // close file

        // Write truth information
        
        info_file << primary.vertex.X() << std::endl;
        info_file << primary.vertex.Y() << std::endl;
        info_file << primary.vertex.Z() << std::endl;
        info_file << primary.energy << std::endl;
        info_file << primary.interaction << std::endl;
        info_file << primary.pdgCode << std::endl;
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

  } // cvn::CVNZlibMakerProtoDUNE::write_files

DEFINE_ART_MODULE(cvn::CVNZlibMakerProtoDUNE)
} // namespace cvn
