////////////////////////////////////////////////////////////////////////
// \file    GCNZlibMakerProtoDUNE_module.cc
// \brief   Analyzer module for creating CVN gzip file objects
// \author  Jeremy Hewes - jhewes15@fnal.gov
//          Saul Alonso-Monsalve - saul.alonso.monsalve@cern.ch
//           - wrote the zlib code used in this module
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>
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
#include "dune/CVN/func/GCNGraph.h"

// Compression
#include "zlib.h"

namespace fs = boost::filesystem;

namespace cvn {

  class GCNZlibMakerProtoDUNE : public art::EDAnalyzer {
  public:

    explicit GCNZlibMakerProtoDUNE(fhicl::ParameterSet const& pset);
    ~GCNZlibMakerProtoDUNE();

    void beginJob() override;
    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);

  private:

    std::string fOutputDir;
    std::string fGraphLabel;
    unsigned int fTopologyHitsCut;

    std::string fLArG4ModuleLabel;

    std::string out_dir;

  };

  //......................................................................
  GCNZlibMakerProtoDUNE::GCNZlibMakerProtoDUNE(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  GCNZlibMakerProtoDUNE::~GCNZlibMakerProtoDUNE()
  {  }

  //......................................................................
  void GCNZlibMakerProtoDUNE::reconfigure(const fhicl::ParameterSet& pset)
  {
    fOutputDir = pset.get<std::string>("OutputDir", "");
    fGraphLabel = pset.get<std::string>("GraphLabel");
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");

    fLArG4ModuleLabel = pset.get<std::string>("LArG4ModuleLabel");
  }

  //......................................................................
  void GCNZlibMakerProtoDUNE::beginJob()
  {
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
  void GCNZlibMakerProtoDUNE::analyze(const art::Event& evt)
  {

    std::cout << "GCNZlibMakerProtoDUNE: looking for graphs with label " << fGraphLabel << std::endl;

    // Get the graphs
    std::vector<art::Ptr<cvn::GCNGraph>> graphs;
    auto h_graphs = evt.getHandle<std::vector<cvn::GCNGraph>>(fGraphLabel);
    if (h_graphs)
      art::fill_ptr_vector(graphs, h_graphs);

    // If no graphs, quit
    if (graphs.size() == 0) return;

    // MC information
    art::ServiceHandle< cheat::ParticleInventoryService > piService;

    AssignLabels labels;
    unsigned short beamParticleInteraction;

    // We only have one beam particle->called "primary", so find it. It should be the
    // first particle-> but best not to assume
    float beamParticleEnergy = 0;
    TVector3 beamParticleVtx; // This is the interaction vertex
    int beamParticlePDG = 0;

//    bool gotPrimary = false;
    for(auto const m: piService->ParticleList()){
      const simb::MCParticle* particle = m.second;
      if(particle->Process().compare("primary")==0){
        beamParticleEnergy = particle->E();
        beamParticleVtx.SetXYZ(particle->EndX(),particle->EndY(),particle->EndZ());
        beamParticlePDG = particle->PdgCode();
        beamParticleInteraction = labels.GetProtoDUNEBeamInteractionType(*particle);
//        gotPrimary = true;
        break;
      }
    }

//    if(!gotPrimary) return;
    unsigned int counter = 0;
    for(auto const graph : graphs){
      // If the graph has no nodes then give up
      if(graph->GetNumberOfNodes() == 0) return;

//      std::cout << "GCNZlibMakerProtoDUNE: found graph with " << graph->GetNumberOfNodes() << " nodes" << std::endl;


      // Now write the zlib file using this information
      // We need to extract all of the information into a single vector to write
      // into the compressed file format
      const std::vector<float> vectorToWrite = graph->ConvertGraphToVector();
 
      ulong src_len = vectorToWrite.size() *sizeof(float);
      ulong dest_len = compressBound(src_len);     // calculate size of the compressed data               
      char* ostream = (char *) malloc(dest_len);  // allocate memory for the compressed data

      int res = compress((Bytef *) ostream, &dest_len, (Bytef *) &vectorToWrite[0], src_len);

      // Buffer error
      if (res == Z_BUF_ERROR)
        std::cout << "Buffer too small!" << std::endl;
      // Memory error
      else if (res ==  Z_MEM_ERROR)
        std::cout << "Not enough memory for compression!" << std::endl;
      // Compression ok 
      else {

        // Create output files 
        std::stringstream image_file_name; 
        image_file_name << out_dir << "/gcn_event_" << evt.event() << "_" << counter << ".gz";
        std::stringstream info_file_name;
        info_file_name  << out_dir << "/gcn_event_" << evt.event() << "_" << counter << ".info";

        std::ofstream image_file (image_file_name.str(), std::ofstream::binary);
        std::ofstream info_file  (info_file_name.str());

        if(image_file.is_open() && info_file.is_open()) {

          // Write the graph to the file and close it
          image_file.write(ostream, dest_len);
          image_file.close(); // close file

          // Write the auxillary information to the text file
          info_file << beamParticleVtx.X() << std::endl;
          info_file << beamParticleVtx.Y() << std::endl;
          info_file << beamParticleVtx.Z() << std::endl;
          info_file << beamParticleEnergy << std::endl;
          info_file << beamParticleInteraction << std::endl; // Interaction type first
          info_file << beamParticlePDG << std::endl;

          // Number of nodes and node features is needed for unpacking
          info_file << graph->GetNumberOfNodes() << std::endl;
          info_file << graph->GetNode(0).GetNumberOfFeatures() << std::endl;

          info_file.close(); // close file
        }
        else {

          if (image_file.is_open())
            image_file.close();
          else 
            throw art::Exception(art::errors::FileOpenError)
              << "Unable to open file " << image_file_name.str() << "!" << std::endl;

          if (info_file.is_open())
            info_file.close();
          else
            throw art::Exception(art::errors::FileOpenError)
              << "Unable to open file " << info_file_name.str() << "!" << std::endl;
        }
      }
      ++counter; 
      free(ostream);  // free allocated memory
    }
  }
    
DEFINE_ART_MODULE(cvn::GCNZlibMakerProtoDUNE)
} // namespace cvn
