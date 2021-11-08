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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNFeatureUtils.h"

#include "dune/CVN/func/CVNProtoDUNEUtils.h"

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
      bool fUseAllSlices;

      /// Slicing module (typically Pandora)
      std::string fSliceLabel;

      /// PFParticle module (typically Pandora)
      std::string fParticleLabel;

      // Consider EM activity as different to the parent
      bool fUseEM;

      // Use hits not energy for truth matching
      bool fUseHitsForTruthMatching;
  };



  //.......................................................................
  GCNGraphMakerProtoDUNE::GCNGraphMakerProtoDUNE(fhicl::ParameterSet const& pset): art::EDProducer(pset),
  fSpacePointLabel (pset.get<std::string>    ("SpacePointLabel")),
  fMinClusterHits  (pset.get<unsigned short> ("MinClusterHits")),
  fNeighbourRadii (pset.get<std::vector<float>>("NeighbourRadii")),
  fUseBeamSliceOnly(pset.get<bool>("UseBeamSliceOnly")),
  fUseAllSlices(pset.get<bool>("UseAllSlices")),
  fSliceLabel      (pset.get<std::string>("SliceModuleLabel")),
  fParticleLabel   (pset.get<std::string>("ParticleModuleLabel")),
  fUseEM           (pset.get<bool>("UseEM",true)),
  fUseHitsForTruthMatching (pset.get<bool>("UseHitsForTruthMatching",true))
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
    std::vector<art::Ptr<recob::SpacePoint>> allSpacePoints;
    auto spacePointHandle = evt.getHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
    if (spacePointHandle){
      art::fill_ptr_vector(allSpacePoints, spacePointHandle);
    }

    // Graph space points for each slice we want to consider
    std::map<unsigned int,std::map<unsigned int,art::Ptr<recob::SpacePoint>>> allGraphSpacePoints;

    // Get the utility to help us calculate features
    cvn::GCNFeatureUtils graphUtil;

    // We can calculate the number of neighbours for each space point with some radius
    // Store the number of neighbours for each spacepoint ID, for each slice. If we
    // only want the beam slice, or all of the slices together, this will have one
    // element
    std::map<unsigned int,std::vector<std::map<int,unsigned int>>> neighbourMap;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    if(fUseAllSlices || fUseBeamSliceOnly){
      // We need to get a vector of space points from the beam slice
      cvn::CVNProtoDUNEUtils pfpUtil;

      const unsigned short beamSliceIndex = pfpUtil.GetBeamSlice(evt,fParticleLabel);
      if(beamSliceIndex > 500 && fUseBeamSliceOnly){
        std::cerr << "Requested beam slice only yet one didn't exist, returning no graphs" << std::endl; 
        evt.put(std::move(graphs));
        return;
      }

      // Get a map of slice index to all of the PFParticles within it
      const std::map<unsigned int,std::vector<const recob::PFParticle*>> particleSliceMap = pfpUtil.GetAllPFParticleSliceMap(evt, fParticleLabel);
      for(const std::pair<unsigned int,std::vector<const recob::PFParticle*>> m : particleSliceMap){
        unsigned int sliceID = m.first;
        std::map<unsigned int,art::Ptr<recob::SpacePoint>> sliceSpacePoints;
        // If we want the beam slice then make sure we have it
        if(fUseBeamSliceOnly && (sliceID != beamSliceIndex)) continue;

        for(const recob::PFParticle* p : m.second){
          // Get the SpacePoints associated to the PFParticle
          const std::vector<const recob::SpacePoint*> particlePoints = pfpUtil.GetPFParticleSpacePoints(*p, evt, fParticleLabel);
          //std::cout << "Got beam slice particle with " << particlePoints.size() << " space points" << std::endl;
          for(const recob::SpacePoint* s : particlePoints){
            sliceSpacePoints.insert(std::make_pair(s->ID(),allSpacePoints.at(s->ID())));
          }
        }
        allGraphSpacePoints[sliceID] = sliceSpacePoints;
        neighbourMap.insert(std::make_pair(sliceID,graphUtil.GetNeighboursForRadii(evt,fNeighbourRadii,sliceSpacePoints)));
      }
    }
    else{
      std::vector<art::Ptr<recob::SpacePoint>> eventSpacePoints;
      spacePointHandle = evt.getHandle<std::vector<recob::SpacePoint>>(fSpacePointLabel);
      if (spacePointHandle){
        art::fill_ptr_vector(eventSpacePoints, spacePointHandle);
      }
      //Convert this vector to a map
      std::map<unsigned int,art::Ptr<recob::SpacePoint>> mapVec;
      for(art::Ptr<recob::SpacePoint> p : eventSpacePoints){
        mapVec.insert(std::make_pair(p->ID(),p));
      }
      neighbourMap.insert(std::make_pair(0,graphUtil.GetNeighboursForRadii(evt,fNeighbourRadii,eventSpacePoints)));
      allGraphSpacePoints.insert(std::make_pair(0,mapVec));
    }

    // Get some maps we'll need to fill our features
    std::cout << "Found all neighbours for " << neighbourMap.size() << " slices, building graphs..." << std::endl;

    // Function is linear in number of points so just do it once
    std::map<unsigned int,float> chargeMap = graphUtil.GetSpacePointChargeMap(evt,fSpacePointLabel);

    // Mean hit RMS for each spacepoint
    std::map<unsigned int, float> hitRMSMap = graphUtil.GetSpacePointMeanHitRMSMap(evt,fSpacePointLabel);

    // The true particle PDG code is needed for training node classifiers
    std::map<unsigned int,int> trueIDMap = graphUtil.GetTruePDG(clockData, evt, fSpacePointLabel, !fUseEM, fUseHitsForTruthMatching);

    // Now we want to produce a graph for each one of the slices
    for(const std::pair<unsigned int,std::map<unsigned int,art::Ptr<recob::SpacePoint>>> &sps : allGraphSpacePoints){
      if(sps.second.size() >= fMinClusterHits){
       
        cvn::GCNGraph newGraph;

        std::map<int,std::pair<int,int>> twoNearest = graphUtil.GetTwoNearestNeighbours(evt,sps.second);

        std::cout << "Constructing graph for slice " << sps.first << " with " << sps.second.size() << " nodes." << std::endl;
        for(const std::pair<unsigned int,art::Ptr<recob::SpacePoint>> &sp : sps.second){
  
          // Get the position
          std::vector<float> position;
          // Why does this use an array... we want a vector in any case
          const double *pos = sp.second->XYZ();
          for(unsigned int p = 0; p < 3; ++p) position.push_back(pos[p]);
  
          // Calculate some features
          std::vector<float> features;
  
          // The neighbour map gives us our first feature(s)
          for(unsigned int m = 0; m < fNeighbourRadii.size(); ++m){
            features.push_back(neighbourMap.at(sps.first)[m].at(sp.second->ID()));
          }
  
          // How about charge?
          features.push_back(chargeMap.at(sp.second->ID()));
  
          // Now the hit width
//          features.push_back(hitRMSMap.at(sp.second->ID()));

          // Angle and dot product between node and its two nearest neighbours
          float angle = -999.;
          float dotProduct = -999.;
          int n1ID = twoNearest[sp.second->ID()].first;
          int n2ID = twoNearest[sp.second->ID()].second;
          const recob::SpacePoint n1 = *(sps.second.at(n1ID).get());
          const recob::SpacePoint n2 = *(sps.second.at(n2ID).get());
          graphUtil.GetAngleAndDotProduct(*(sp.second.get()),n1,n2,dotProduct,angle);
          features.push_back(dotProduct);
          features.push_back(angle);
  
          // We set the "ground truth" as the particle PDG code in this case
          std::vector<float> truePDG;
          truePDG.push_back(static_cast<float>(trueIDMap.at(sp.second->ID())));
          newGraph.AddNode(position,features,truePDG);
//          if(abs(trueIDMap.at(sp.second->ID())) != 13) std::cout << "Adding node " << sp.second->ID() << " with neighbours " << n1ID << " and " << n2ID << " and PDG = " << truePDG[0] << std::endl;
        }
  
        std::cout << "GCNGraphMakerProtoDUNE: produced GCNGraph object with " << newGraph.GetNumberOfNodes() << " nodes" << std::endl;
  
        // Add out graph to the vector
        graphs->push_back(newGraph);
      }
    }

    // Write our graph to the event
    evt.put(std::move(graphs));
  }

  DEFINE_ART_MODULE(cvn::GCNGraphMakerProtoDUNE)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
