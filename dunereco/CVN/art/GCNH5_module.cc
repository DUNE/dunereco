////////////////////////////////////////////////////////////////////////
// Class:       GCNH5
// Plugin Type: analyzer (art v3_01_02)
// File:        GCNH5_module.cc
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
#include "dune/CVN/func/GCNGraph.h"

// highfive includes
#include <dune/CVN/highfive/H5DataSet.hpp>
#include <dune/CVN/highfive/H5DataSpace.hpp>
#include <dune/CVN/highfive/H5File.hpp>

namespace cvn {

  class GCNH5 : public art::EDAnalyzer {
  public:
    explicit GCNH5(fhicl::ParameterSet const& p);

    GCNH5(GCNH5 const&) = delete;
    GCNH5(GCNH5&&) = delete;
    GCNH5& operator=(GCNH5 const&) = delete;
    GCNH5& operator=(GCNH5&&) = delete;

    void reconfigure(fhicl::ParameterSet const& p);

    void endSubRun(art::SubRun const& sr) override;
    void analyze(art::Event const& e) override;

  private:

    std::string fGraphModuleLabel; /// Name of graph producer module
    std::string fOutputName; /// H5 output filename

    std::vector<std::vector<unsigned int>> fEvents; /// Event numbers
    std::vector<float> fGraphs; /// Output graph vector
    unsigned int fNFeatures; /// Number of features in each node
    std::vector<unsigned int> fNNodes; /// Number of nodes in each graph

  };


  GCNH5::GCNH5(fhicl::ParameterSet const& p)
    : EDAnalyzer{p} {

    this->reconfigure(p);

  }

  void GCNH5::reconfigure(fhicl::ParameterSet const& p) {

    fGraphModuleLabel = p.get<std::string>("GraphModuleLabel");
    fOutputName = p.get<std::string>("OutputName");

  } // cvn::GCNH5::reconfigure

  void GCNH5::analyze(art::Event const& e) {

    // Get the graphs
    art::Handle<std::vector<GCNGraph>> h_graphs;
    std::vector<art::Ptr<GCNGraph>> graphs;
    if (e.getByLabel(fGraphModuleLabel, h_graphs))
      art::fill_ptr_vector(graphs, h_graphs);

    for (auto graph : graphs) {
      fEvents.push_back(std::vector<unsigned int>({e.id().run(), e.id().subRun(), e.id().event()}));
      std::vector<float> graph_vec = graph->ConvertGraphToVector();
      fGraphs.insert(fGraphs.end(), graph_vec.begin(), graph_vec.end());
      fNFeatures = graph->GetNode(0).GetNumberOfFeatures() + 3; // Assume same number of features in each node and graph... is that a safe assumption?
      fNNodes.push_back(graph->GetNumberOfNodes());
    }
    

  } // cvn::GCNH5::analyze

  /// End of a subrun, write all graphs to a H5 file
  void GCNH5::endSubRun(art::SubRun const& sr) {

    // Ignore empty subruns
    if (fGraphs.empty()) return;

    try {

      // Open H5 file
      std::ostringstream file_name;
      file_name << fOutputName << "_" << sr.id().run() << "_" << sr.id().subRun() << ".h5";
      HighFive::File f(file_name.str(),
        HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

      // Write images
      f.createDataSet("graphvector", fGraphs);

      // Write metadata
      f.createDataSet("events", fEvents);
      f.createDataSet("n_nodes", fNNodes);
      f.createDataSet("n_features", fNFeatures);

    } catch (HighFive::Exception& err) {
      std::cerr << err.what() << std::endl;
    }

    fEvents.clear();
    fGraphs.clear();
    fNNodes.clear();

  } // cvn::GCNH5::endSubRun

  DEFINE_ART_MODULE(GCNH5)

} // namespace cvn
