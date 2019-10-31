////////////////////////////////////////////////////////////////////////
// \file    GCNGraphMakerProtoDUNE_module.cc
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

// LArSoft includes
#include "lardataobj/RecoBase/SpacePoint.h"

#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNFeatureUtils.h"

#include "dune/Protodune/Analysis/ProtoDUNEPFParticleUtils.h"
#include "dune/Protodune/Analysis/ProtoDUNESliceUtils.h"

namespace cvn {

  class GCNGraphMakerProtoDUNE : public art::EDProducer {
    public:
      explicit GCNGraphMakerProtoDUNE(fhicl::ParameterSet const& pset);
      ~GCNGraphMakerProtoDUNE();

      void produce(art::Event& evt);
      void beginJob();
      void endJob();



    private:
      /// Module label for input space points
      std::string fSpacePointLabel;

      /// Minimum number of space points to produce a graph
      unsigned short fMinClusterHits;

      /// Radii for calculating number of neighbours for any number of cut values
      std::vector<float> fNeighbourRadii;

      /// Use only the beam slice as determined by Pandora
      bool fUseBeamSliceOnly;

      /// Slicing module (typically Pandora)
      std::string fSliceLabel;

      /// PFParticle module (typically Pandora)
      std::string fParticleLabel;
  };



  //.......................................................................
  GCNGraphMakerProtoDUNE::GCNGraphMakerProtoDUNE(fhicl::ParameterSet const& pset): art::EDProducer(pset),
  fSpacePointLabel (pset.get<std::string>    ("SpacePointLabel")),
  fMinClusterHits  (pset.get<unsigned short> ("MinClusterHits")),
  fNeighbourRadii (pset.get<std::vector<float>>("NeighbourRadii")),
  fUseBeamSliceOnly(pset.get<bool>("UseBeamSliceOnly")),
  fSliceLabel      (pset.get<std::string>("SliceModuleLabel")),
  fParticleLabel   (pset.get<std::string>("ParticleModuleLabel"))
  {

    produces< std::vector<cvn::GCNGraph> >();

  }

  //......................................................................
  GCNGraphMakerProtoDUNE::~GCNGraphMakerProtoDUNE()
  {
    //======================================================================
    // Clean up any memory allocated by your module
    //======================================================================
  }

  //......................................................................
  void GCNGraphMakerProtoDUNE::beginJob()
  {  }

  //......................................................................
  void GCNGraphMakerProtoDUNE::endJob()
  {
  }

  //......................................................................
  void GCNGraphMakerProtoDUNE::produce(art::Event& evt)
  {


    // Create the Graph vector and fill it if we have enough hits
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>);
    //    std::cout << "GCNGraphMakerProtoDUNE: checking if we have enough space points (" << pointList.size() << " / " << fMinClusterHits << ")" << std::endl;
    // Vector for all of the space points
    // Get the space points from the event
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    if(evt.getByLabel(fSpacePointLabel,spacePointHandle)){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }
    std::vector<art::Ptr<recob::SpacePoint>> graphSpacePoints;

    cvn::GCNGraph newGraph;

    // Get the utility to help us calculate features
    cvn::GCNFeatureUtils graphUtil;

    // We can calculate the number of neighbours for each space point with some radius
    // Store the number of neighbours for each spacepoint ID
    std::vector<std::map<int,unsigned int>> neighbourMap;
    // Get the charge map too
    std::map<unsigned int,float> chargeMap;

    if(fUseBeamSliceOnly){
      // We need to get a vector of space points from the beam slice
      protoana::ProtoDUNEPFParticleUtils pfpUtil;
      const unsigned short beamSliceIndex = pfpUtil.GetBeamSlice(evt,fParticleLabel);

      if(beamSliceIndex < 500){

        // Get a map of slice index to all of the PFParticles within it
        const std::map<unsigned int,std::vector<const recob::PFParticle*>> particleSliceMap = pfpUtil.GetAllPFParticleSliceMap(evt, fParticleLabel);
        for(const recob::PFParticle* p : particleSliceMap.at(beamSliceIndex)){
          // Get the SpacePoints associated to the PFParticle
          const std::vector<const recob::SpacePoint*> particlePoints = pfpUtil.GetPFParticleSpacePoints(*p, evt, fParticleLabel);
          std::cout << "Got beam slice particle with " << particlePoints.size() << " space points" << std::endl;
          for(const recob::SpacePoint* s : particlePoints){
            graphSpacePoints.push_back(allSpacePoints.at(s->ID()));
          }
        }
        neighbourMap  = graphUtil.GetNeighboursForRadii(evt,fNeighbourRadii,graphSpacePoints);
      }
    }
    else{
      if(evt.getByLabel(fSpacePointLabel,spacePointHandle)){
        art::fill_ptr_vector(graphSpacePoints, spacePointHandle);
      }
      neighbourMap = graphUtil.GetNeighboursForRadii(evt,fNeighbourRadii,graphSpacePoints);
    }

    if(graphSpacePoints.size() != 0 && graphSpacePoints.size() >= fMinClusterHits){
      
      chargeMap = graphUtil.GetSpacePointChargeMap(evt,fParticleLabel);
      std::map<int,std::pair<int,int>> twoNearest = graphUtil.GetTwoNearestNeighbours(evt,fSpacePointLabel);

      for(art::Ptr<recob::SpacePoint> sp : graphSpacePoints){

        // Get the position
        std::vector<float> position;
        // Why does this use an array... we want a vector in any case
        const double *pos = sp->XYZ();
        for(unsigned int p = 0; p < 3; ++p) position.push_back(pos[p]);

        // Calculate some features
        std::vector<float> features;

        // The neighbour map gives us our first feature(s)
        for(unsigned int m = 0; m < fNeighbourRadii.size(); ++m){
          features.push_back(neighbourMap[m].at(sp->ID()));
        }

        // How about charge?
        features.push_back(chargeMap.at(sp->ID()));

        // Angle and dot product between node and its two nearest neighbours
        float angle = -999.;
        float dotProduct = -999.;
        int n1ID = twoNearest[sp->ID()].first;
        int n2ID = twoNearest[sp->ID()].second;
        const recob::SpacePoint &n1 = *(graphSpacePoints[n1ID].get());
        const recob::SpacePoint &n2 = *(graphSpacePoints[n2ID].get());
        graphUtil.GetAngleAndDotProduct(*(sp.get()),n1,n2,dotProduct,angle);
        features.push_back(dotProduct);
        features.push_back(angle);

        newGraph.AddNode(position,features);
      }

      std::cout << "GCNGraphMakerProtoDUNE: produced GCNGraph object with " << newGraph.GetNumberOfNodes() << " nodes" << std::endl;

      // Add out graph to the vector
      graphs->push_back(newGraph);
    }

    // Write our graph to the event
    evt.put(std::move(graphs));
  }


  DEFINE_ART_MODULE(cvn::GCNGraphMakerProtoDUNE)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////







