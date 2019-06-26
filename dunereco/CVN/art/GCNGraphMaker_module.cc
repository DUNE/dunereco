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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

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

    /// Do we want collection plane hits only?
    bool fCollectionPlaneOnly;

    /// Include true GEANT ID of primary associated particle as graph feature
    bool fTrueParticleFeature;

  };



  //.......................................................................
  GCNGraphMaker::GCNGraphMaker(fhicl::ParameterSet const& pset): art::EDProducer(pset),
  fSpacePointLabel (pset.get<std::string>    ("SpacePointLabel")),
  fMinClusterHits  (pset.get<unsigned short> ("MinClusterHits")),
  fNeighbourRadius (pset.get<float>("NeighbourRadius")),
  fCollectionPlaneOnly(pset.get<bool>("CollectionPlaneOnly")),
  fTrueParticleFeature(pset.get<bool>("TrueParticleFeature"))

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
    auto allSP = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
    const art::FindManyP<recob::Hit> sp2Hit(allSP, evt, fSpacePointLabel);

    // Create the Graph vector and fill it if we have enough hits
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>);

    if(allSP->size() >= fMinClusterHits){

      cvn::GCNGraph newGraph;

      // Get the utility to help us calculate features
      cvn::GCNFeatureUtils graphUtil;

      // We can calculate the number of neighbours for each space point with some radius
      // Store the number of neighbours for each spacepoint ID
      const std::map<int, unsigned int> neighbourMap = graphUtil.GetAllNeighbours(evt, fNeighbourRadius, fSpacePointLabel);

      // Get the charge and true ID for each spacepoint
      auto chargeMap = graphUtil.GetSpacePointChargeMap(evt, fSpacePointLabel);
      const std::map<unsigned int, unsigned int> *trueIDMap = nullptr;
      if (fTrueParticleFeature) trueIDMap = new const std::map<unsigned int, unsigned int>(graphUtil.GetTrueG4ID(evt, fSpacePointLabel));

      for (const recob::SpacePoint &sp : *allSP) {

        // Do we only want collection plane spacepoints?
        if (fCollectionPlaneOnly) {
          // Are there any associated hits from the collection plane?
          bool collectionHit = false;
          for (auto hit : sp2Hit.at(sp.ID())) {
            if (hit->View() == 2) {
              collectionHit = true;
              break;
            }
          }
          // If not, skip this spacepoint
          if (!collectionHit) continue;
        }

        // Get the position
        std::vector<float> position;
        // Why does this use an array... we want a vector in any case
        const double *pos = sp.XYZ();
        for(unsigned int p = 0; p < 3; ++p) position.push_back(pos[p]);
        
        // Calculate some features
        std::vector<float> features;
        // The neighbour map gives us our first feature
        features.push_back(neighbourMap.at(sp.ID()));
        // Now charge and true ID
        features.push_back(chargeMap.at(sp.ID()));
        if (fTrueParticleFeature) features.push_back(trueIDMap->at(sp.ID()));

        newGraph.AddNode(position, features);
      }

      std::cout << "GCNGraphMaker: produced GCNGraph object with " << newGraph.GetNumberOfNodes()
      << " nodes from " << allSP->size() << " spacepoints." << std::endl;

      // Add out graph to the vector
      graphs->push_back(newGraph);
    }

    // Write our graph to the event
    evt.put(std::move(graphs));
  }


DEFINE_ART_MODULE(cvn::GCNGraphMaker)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







