////////////////////////////////////////////////////////////////////////
// \file    GCNGraphMaker2D_module.cc
// \brief   Producer module for creating three 2D CVN Graph objects
// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch 
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNFeatureUtils.h"

namespace cvn {

  class GCNGraphMaker2D : public art::EDProducer {
  public:
    explicit GCNGraphMaker2D(fhicl::ParameterSet const& pset);
    ~GCNGraphMaker2D();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    /// Module label for input space points
    std::string fPixelMapLabel;

    /// Minimum charge to created a graph node from a pixel
    float fChargeThreshold;

    /// Radius for calculating number of neighbours
    float fNeighbourPixels;

  };



  //.......................................................................
  GCNGraphMaker2D::GCNGraphMaker2D(fhicl::ParameterSet const& pset): art::EDProducer(pset),
  fPixelMapLabel (pset.get<std::string>    ("PixelMapLabel")),
  fChargeThreshold  (pset.get<float> ("ChargeThreshold")),
  fNeighbourPixels (pset.get<float>("NeighbourPixels"))
  {
    produces< std::vector<cvn::GCNGraph>   >();
  }

  //......................................................................
  GCNGraphMaker2D::~GCNGraphMaker2D()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void GCNGraphMaker2D::beginJob()
  {  }

  //......................................................................
  void GCNGraphMaker2D::endJob()
  {
  }

  //......................................................................
  void GCNGraphMaker2D::produce(art::Event& evt)
  {
    auto pixelMaps = evt.getValidHandle<std::vector<cvn::PixelMap>>(fPixelMapLabel);

    // Create the Graph vector and fill it if we have enough hits
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>);

    if(pixelMaps->size() > 0){

      cvn::GCNFeatureUtils graphUtil;

      // Get the three graphs from the pixel map
      // These graphs contain N nodes with (wire,time) position and (charge) features
      std::vector<cvn::GCNGraph> graphs2D = graphUtil.ExtractGraphsFromPixelMap(pixelMaps->at(0),fChargeThreshold);

      // For each of the graphs we want to add a number of neighbours feature
      for(cvn::GCNGraph g : graphs2D){
        std::map<unsigned int, unsigned int> neighbourMap = graphUtil.Get2DGraphNeighbourMap(g,fNeighbourPixels);
        std::cout << "Built graph with " << g.GetNumberOfNodes() << " nodes" << std::endl;
        for(unsigned int n = 0; n < g.GetNumberOfNodes(); ++n){
          cvn::GCNGraphNode &node = g.GetNodeEditable(n);
          node.AddFeature(neighbourMap.at(n));
        }
        // Add the graph to the output vector
        graphs->push_back(g);
      } 
      
    }
    // Write our graph to the event
    evt.put(std::move(graphs));
  }

  DEFINE_ART_MODULE(cvn::GCNGraphMaker2D)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







