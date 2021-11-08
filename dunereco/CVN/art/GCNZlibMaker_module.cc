////////////////////////////////////////////////////////////////////////
// \file    GCNZlibMaker_module.cc
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
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/Exception.h"

// Data products
#include "nusimdata/SimulationBase/MCTruth.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "dune/CVN/func/AssignLabels.h"
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/InteractionType.h"

// Compression
#include "zlib.h"

namespace fs = boost::filesystem;

namespace cvn {

  class GCNZlibMaker : public art::EDAnalyzer {
  public:

    explicit GCNZlibMaker(fhicl::ParameterSet const& pset);
    ~GCNZlibMaker();

    void beginJob() override;
    void analyze(const art::Event& evt) override;
    void reconfigure(const fhicl::ParameterSet& pset);

  private:

    std::string fOutputDir;
    std::string fGraphLabel;
    unsigned int fTopologyHitsCut;

    std::string fGenieGenModuleLabel;
    std::string fEnergyNueLabel;
    std::string fEnergyNumuLabel;
    std::string fEnergyNutauLabel;

    std::string out_dir;

  };

  //......................................................................
  GCNZlibMaker::GCNZlibMaker(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  GCNZlibMaker::~GCNZlibMaker()
  {  }

  //......................................................................
  void GCNZlibMaker::reconfigure(const fhicl::ParameterSet& pset)
  {
    fOutputDir = pset.get<std::string>("OutputDir", "");
    fGraphLabel = pset.get<std::string>("GraphLabel");
    fTopologyHitsCut = pset.get<unsigned int>("TopologyHitsCut");

    fGenieGenModuleLabel = pset.get<std::string>("GenieGenModuleLabel");
    fEnergyNueLabel = pset.get<std::string>("EnergyNueLabel");
    fEnergyNumuLabel = pset.get<std::string>("EnergyNumuLabel");
    fEnergyNutauLabel = pset.get<std::string>("EnergyNutauLabel");
  }

  //......................................................................
  void GCNZlibMaker::beginJob()
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
  void GCNZlibMaker::analyze(const art::Event& evt)
  {

    std::cout << "GCNZlibMaker: looking for graphs with label " << fGraphLabel << std::endl;

    // Get the graphs
    std::vector<art::Ptr<cvn::GCNGraph>> graphs;
    auto h_graphs = evt.getHandle<std::vector<cvn::GCNGraph>>(fGraphLabel);
    if (h_graphs)
      art::fill_ptr_vector(graphs, h_graphs);

    // If no graphs, quit
    std::cout << "Found " << graphs.size() << " graphs" << std::endl;
    if (graphs.size() == 0) return;

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

    // We only want interactions in the fiducial volume for training
    // Get the interaction vertex from the end point of the neutrino. This is 
    // because the start point of the lepton doesn't make sense for taus as they
    // are decayed by the generator and not GEANT
    TVector3 vtx = true_neutrino.Nu().EndPosition().Vect();
    bool isFid = (fabs(vtx.X())<310. && fabs(vtx.Y())<550. && vtx.Z()>50. && vtx.Z()<1244.);
    if(!isFid) return;

    float reco_nue_energy = 0.;
    float reco_numu_energy = 0.;
    float reco_nutau_energy = 0.;

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

    float event_weight = -1.0; // We currently just set this to -1

    // If the graph has no nodes then give up
    for(unsigned int i = 0; i < graphs.size(); ++i){
      const art::Ptr<cvn::GCNGraph> g = graphs.at(i);
      
      std::cout << "GCNZlibMaker: found graph with " << g->GetNumberOfNodes() << " nodes" << std::endl;
      if(g->GetNumberOfNodes() == 0) continue;
  
      // Now write the zlib file using this information
      // We need to extract all of the information into a single vector to write
      // into the compressed file format
      const std::vector<float> vectorToWrite = g->ConvertGraphToVector();
   
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

        std::stringstream modifier;
        if(graphs.size() > 1){
          modifier << "_" << i;
        }
  
        // Create output files 
        std::stringstream image_file_name; 
        image_file_name << out_dir << "/event_" << evt.event() << modifier.str() << ".gz";
        std::stringstream info_file_name;
        info_file_name << out_dir << "/event_" << evt.event() << modifier.str() << ".info";
  
        std::ofstream image_file (image_file_name.str(), std::ofstream::binary);
        std::ofstream info_file  (info_file_name.str());
  
        if(image_file.is_open() && info_file.is_open()) {
  
          // Write the graph to the file and close it
          image_file.write(ostream, dest_len);
          image_file.close(); // close file
  
          // Write the auxillary information to the text file
          info_file << interaction << std::endl; // Interaction type first
  
          // True and reconstructed energy variables
          info_file << nu_energy << std::endl;
          info_file << lep_energy << std::endl;
          info_file << reco_nue_energy << std::endl;
          info_file << reco_numu_energy << std::endl;
          info_file << reco_nutau_energy << std::endl;
          info_file << event_weight << std::endl;
  
          info_file << labels.GetPDG() << std::endl;
          info_file << labels.GetNProtons() << std::endl;
          info_file << labels.GetNPions() << std::endl;
          info_file << labels.GetNPizeros() << std::endl;
          info_file << labels.GetNNeutrons() << std::endl;
          info_file << labels.GetTopologyType() << std::endl;
          info_file << labels.GetTopologyTypeAlt() << std::endl;
  
          // Number of nodes and node features is needed for unpacking
          info_file << g->GetNumberOfNodes() << std::endl;
          info_file << g->GetNumberOfNodeCoordinates() << std::endl;
          info_file << g->GetNumberOfNodeFeatures() << std::endl;
  
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
      free(ostream);  // free allocated memory
      ostream = nullptr;
    }
    
    return;
  }
    
DEFINE_ART_MODULE(cvn::GCNZlibMaker)
} // namespace cvn
