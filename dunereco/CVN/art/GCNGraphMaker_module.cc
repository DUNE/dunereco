////////////////////////////////////////////////////////////////////////
// \file    GCNGraphMaker_module.cc
// \brief   Producer module for creating CVN Graph objects
// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch 
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"

// LArSoft includes
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNFeatureUtils.h"

namespace cvn {

  class GCNGraphMaker : public art::EDProducer {
  public:
    explicit GCNGraphMaker(fhicl::ParameterSet const& pset);
    ~GCNGraphMaker();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    /// Module label for input space points
    std::string fSpacePointLabel;

    /// Minimum number of space points to produce a graph
    unsigned short fMinClusterHits;

    /// Radius for calculating number of neighbours
    float fNeighbourRadius;

  };



  //.......................................................................
  GCNGraphMaker::GCNGraphMaker(fhicl::ParameterSet const& pset):
  fSpacePointLabel (pset.get<std::string>    ("SpacePointLabel")),
  fMinClusterHits  (pset.get<unsigned short> ("MinClusterHits")),
  fNeighbourRadius (pset.get<float>("NeighbourRadius"))
  {

    produces< std::vector<cvn::GCNGraph>   >();

  }

  //......................................................................
  GCNGraphMaker::~GCNGraphMaker()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void GCNGraphMaker::beginJob()
  {  }

  //......................................................................
  void GCNGraphMaker::endJob()
  {
  }

  //......................................................................
  void GCNGraphMaker::produce(art::Event& evt)
  {

    // Get the space points from the event
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> pointList;
    if(evt.getByLabel(fSpacePointLabel,spacePointHandle)){
      art::fill_ptr_vector(pointList, spacePointHandle);
    }

    // Create the Graph vector and fill it if we have enough hits
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>);
//    std::cout << "GCNGraphMaker: checking if we have enough space points (" << pointList.size() << " / " << fMinClusterHits << ")" << std::endl;
    if(pointList.size() >= fMinClusterHits){

      cvn::GCNGraph newGraph;

      // Get the utility to help us calculate features
      cvn::GCNFeatureUtils graphUtil;

      // We can calculate the number of neighbours for each space point with some radius
      // Store the number of neighbours for each spacepoint ID
      const std::map<int,unsigned int> neighbourMap = graphUtil.GetAllNeighbours(evt,fNeighbourRadius,fSpacePointLabel);

      for(art::Ptr<recob::SpacePoint> sp : pointList){

        // Get the position
        std::vector<float> position;
        // Why does this use an array... we want a vector in any case
        const double *pos = sp->XYZ();
        for(unsigned int p = 0; p < 3; ++p) position.push_back(pos[p]);
        
        // Calculate some features
        std::vector<float> features;
        // The neighbour map gives us our first feature
        features.push_back(neighbourMap.at(sp->ID()));
        // How about charge?
        features.push_back(graphUtil.GetSpacePointCharge(*sp,evt,fSpacePointLabel));

        newGraph.AddNode(position,features);
      }

      std::cout << "GCNGraphMaker: produced GCNGraph object with " << newGraph.GetNumberOfNodes() << " nodes" << std::endl;

      // Add out graph to the vector
      graphs->push_back(newGraph);
    }

    // Write our graph to the event
    evt.put(std::move(graphs));
  }


DEFINE_ART_MODULE(cvn::GCNGraphMaker)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







