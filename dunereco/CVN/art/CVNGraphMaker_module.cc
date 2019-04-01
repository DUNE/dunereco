////////////////////////////////////////////////////////////////////////
// \file    CVNGraphMaker_module.cc
// \brief   Producer module for creating CVN Graph objects
// \author  Leigh H. Whitehead - leigh.howard.whitehead@cern.ch 
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// ROOT includes
#include "TTree.h"
#include "TH2F.h"

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
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/TrainingData.h"

namespace cvn {

  class CVNGraphMaker : public art::EDProducer {
  public:
    explicit CVNGraphMaker(fhicl::ParameterSet const& pset);
    ~CVNGraphMaker();

    void produce(art::Event& evt);
    void beginJob();
    void endJob();



  private:
    /// Module label for input space points
    std::string fSpacePointLabel;

    /// Module label for output graphs
    std::string fGraphLabel;

    /// Minimum number of space points to produce a graph
    unsigned short fMinClusterHits;

    /// Radius for calculating number of neighbours
    float fNeighbourRadius;

    float GetDistanceBetweenPoints(const double *p1, const double *p2) const;

    void GetNeighbours(std::map<int,unsigned int> &neighbourMap, const std::vector<art::Ptr<recob::SpacePoint>> &points) const;
  };



  //.......................................................................
  CVNGraphMaker::CVNGraphMaker(fhicl::ParameterSet const& pset):
  fSpacePointLabel (pset.get<std::string>    ("SpacePointLabel")),
  fGraphLabel      (pset.get<std::string>    ("GraphLabel")),
  fMinClusterHits  (pset.get<unsigned short> ("MinClusterHits")),
  fNeighbourRadius (pset.get<float>("NeighbourRadius"))
  {

    produces< std::vector<cvn::PixelMap>   >(fGraphLabel);

  }

  //......................................................................
  CVNGraphMaker::~CVNGraphMaker()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void CVNGraphMaker::beginJob()
  {  }

  //......................................................................
  void CVNGraphMaker::endJob()
  {
  }

  //......................................................................
  void CVNGraphMaker::produce(art::Event& evt)
  {

    // Get the space points from the event
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> pointList;
    if(evt.getByLabel(fSpacePointLabel,spacePointHandle)){
      art::fill_ptr_vector(pointList, spacePointHandle);
    }

    // Create the Graph vector and fill it if we have enough hits
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>);
    if(pointList.size() < fMinClusterHits){

      cvn::GCNGraph newGraph;

      // We can calculate the number of neighbours for each space point with some radius
      // Store the number of neighbours for each spacepoint ID
      std::map<int,unsigned int> neighbourMap;
      this->GetNeighbours(neighbourMap,pointList);

      for(art::Ptr<recob::SpacePoint> sp : pointList){

        // Get the position
        std::vector<float> position;
        // Why does this use an array... we want a vector in any case
        const double *pos = sp->XYZ();
        for(unsigned int p = 0; p < 3; ++p) position.push_back(pos[p]);
        
        // Calculate some features
        std::vector<float> features;
        features.push_back(neighbourMap.at(sp->ID()));

        // Get the features... what will this be?!
        newGraph.AddNode(position,features);
      }

      // Add out graph to the vector
      graphs->push_back(newGraph);
    }

    // Write our graph to the event
    evt.put(std::move(graphs), fGraphLabel);
  }

  //----------------------------------------------------------------------
  void CVNGraphMaker::GetNeighbours(std::map<int,unsigned int> &neighbourMap, const std::vector<art::Ptr<recob::SpacePoint>> &points) const{

    for(art::Ptr<recob::SpacePoint> sp0 : points){
      // We want an entry even if it ends up being zero
      neighbourMap[sp0->ID()] = 0;

      for(art::Ptr<recob::SpacePoint> sp1 : points){

        if(sp0->ID() == sp1->ID()) continue;

        float dist = this->GetDistanceBetweenPoints(sp0->XYZ(),sp1->XYZ());
        if(dist < fNeighbourRadius){
          ++neighbourMap[sp0->ID()];
        }
      }
    }

  }

  float CVNGraphMaker::GetDistanceBetweenPoints(const double *p1, const double *p2) const{
    float dx = p2[0] - p1[0];
    float dy = p2[1] - p1[1];
    float dz = p2[2] - p1[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
  }


DEFINE_ART_MODULE(cvn::CVNGraphMaker)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







