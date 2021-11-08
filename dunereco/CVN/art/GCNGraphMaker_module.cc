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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNParticleFlow.h"
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
    std::string fSpacePointModuleLabel;
    std::string fSpacePointInstanceLabel;

    /// Minimum number of space points to produce a graph
    unsigned short fMinClusterHits;

    /// Number of neighbours as a graph feature - enable & define radius
    bool fUseNeighbourRadius;
    float fNeighbourRadius;

    /// 2D node features
    bool fInclude2DFeatures;

    /// Do we want collection plane hits only?
    bool fCollectionPlaneOnly;

    /// Include true GEANT ID of primary associated particle as graph feature
    bool fSaveTrueParticle;

    /// Include ground truth for node - and if so, define proximity
    bool fUseNodeDeghostingGroundTruth;
    float fTruthRadius;
    bool fUseNodeDirectionGroundTruth;

    /// Whether to save particle hierarchy for particle flow ground truth
    bool fSaveParticleFlow;

  };

  //.......................................................................
  GCNGraphMaker::GCNGraphMaker(fhicl::ParameterSet const& pset): art::EDProducer(pset),
  fSpacePointModuleLabel (pset.get<std::string>    ("SpacePointModuleLabel")),
  fSpacePointInstanceLabel (pset.get<std::string>  ("SpacePointInstanceLabel", "")),
  fMinClusterHits  (pset.get<unsigned short> ("MinClusterHits")),
  fUseNeighbourRadius(pset.get<bool>("UseNeighbourRadius")),
  fNeighbourRadius (pset.get<float>("NeighbourRadius")),
  fInclude2DFeatures (pset.get<bool>("Include2DFeatures")),
  fCollectionPlaneOnly(pset.get<bool>("CollectionPlaneOnly")),
  fSaveTrueParticle(pset.get<bool>("SaveTrueParticle")),
  fUseNodeDeghostingGroundTruth(pset.get<bool>("UseNodeDeghostingGroundTruth")),
  fTruthRadius(pset.get<float>("TruthRadius")),
  fUseNodeDirectionGroundTruth(pset.get<bool>("UseNodeDirectionGroundTruth")),
  fSaveParticleFlow(pset.get<bool>("SaveParticleFlow"))

  {
    produces< std::vector<cvn::GCNGraph>   >();
    if (fSaveParticleFlow) produces< std::vector<cvn::GCNParticleFlow> >();
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
  {}

  //......................................................................
  void GCNGraphMaker::endJob()
  {}

  //......................................................................
  void GCNGraphMaker::produce(art::Event& evt)
  {
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints;
    art::InputTag itag1(fSpacePointModuleLabel,fSpacePointInstanceLabel);
    auto spacePointHandle = evt.getHandle<std::vector<recob::SpacePoint>>(itag1);
    if (!spacePointHandle) {

      throw art::Exception(art::errors::LogicError)
        << "Could not find spacepoints with module label "
        << fSpacePointModuleLabel << " and instance label "
        << fSpacePointInstanceLabel << "!";
    }
    art::fill_ptr_vector(spacePoints, spacePointHandle);
    art::FindManyP<recob::Hit> fmp(spacePointHandle, evt, fSpacePointModuleLabel);
    std::vector<std::vector<art::Ptr<recob::Hit>>> sp2Hit(spacePoints.size());
    for (size_t spIdx = 0; spIdx < sp2Hit.size(); ++spIdx) {
      sp2Hit[spIdx] = fmp.at(spIdx);
    } // for spacepoint

    // Create the Graph vector and fill it if we have enough hits
    std::unique_ptr<std::vector<cvn::GCNGraph>> graphs(new std::vector<cvn::GCNGraph>);
    std::unique_ptr<std::vector<cvn::GCNParticleFlow>> gpf(nullptr);
    if (fSaveParticleFlow) {
      gpf = std::make_unique<std::vector<cvn::GCNParticleFlow> >();
    }

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    if(spacePoints.size() >= fMinClusterHits) {

      cvn::GCNGraph newGraph;

      // Get the utility to help us calculate features
      cvn::GCNFeatureUtils graphUtil;

      // We can calculate the number of neighbours for each space point with some radius
      // Store the number of neighbours for each spacepoint ID
      // const std::map<int, unsigned int> neighbourMap = graphUtil.GetAllNeighbours(evt, fNeighbourRadius, fSpacePointModuleLabel);

      // Get the charge and true ID for each spacepoint
      auto chargeMap = graphUtil.GetSpacePointChargeMap(spacePoints, sp2Hit);
      std::unique_ptr<std::map<unsigned int, int> const> trueIDMap{nullptr};
      if (fSaveTrueParticle) {
        trueIDMap = std::make_unique<std::map<unsigned int, int>>(graphUtil.GetTrueG4ID(clockData, spacePoints, sp2Hit));
      } 

      // Get 2D hit features if requested
      std::map<unsigned int, std::vector<float>> hitMap;
      if (fInclude2DFeatures) {
        hitMap = graphUtil.Get2DFeatures(spacePoints, sp2Hit);
      }

      // Get ground truth if requested
      if (fUseNodeDirectionGroundTruth && !fUseNodeDeghostingGroundTruth) {
        throw art::Exception(art::errors::LogicError)
          << "You must enable deghosting ground truth if using direction ground truth!";
      }
      std::vector<float> nodeDeghostingGroundTruth;
      std::vector<std::vector<float>>* nodeDirectionGroundTruth = nullptr;
      if (fUseNodeDirectionGroundTruth) nodeDirectionGroundTruth = new std::vector<std::vector<float>>();
      if (fUseNodeDeghostingGroundTruth) {
        nodeDeghostingGroundTruth = graphUtil.GetNodeGroundTruth(clockData, spacePoints,
          sp2Hit, fTruthRadius, nodeDirectionGroundTruth);
      }

      std::set<unsigned int> trueParticles;
      for (size_t spIdx = 0; spIdx < spacePoints.size(); ++spIdx) {
        const art::Ptr<recob::SpacePoint> sp = spacePoints[spIdx];
        // Do we only want collection plane spacepoints?

        if (fCollectionPlaneOnly) {
          // Are there any associated hits from the collection plane?
          bool collectionHit = false;
          for (art::Ptr<recob::Hit> hit : sp2Hit[spIdx]) {
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
        const double *pos = sp->XYZ();
        for (size_t p = 0; p < 3; ++p) position.push_back(pos[p]);
        
        // Calculate some features
        std::vector<float> features, truth;
        // The neighbour map gives us our first feature
        // features.push_back(neighbourMap.at(sp->ID()));
        // Now charge and true ID
        features.push_back(chargeMap.at(sp->ID()));

        if (fInclude2DFeatures) {
          features.insert(features.end(), hitMap.at(sp->ID()).begin(), hitMap.at(sp->ID()).end());
        }

        // Now ground truth info
        if (fSaveTrueParticle) {
          truth.push_back(trueIDMap->at(sp->ID()));
          trueParticles.insert(abs(trueIDMap->at(sp->ID())));
        }

        // Add deghosting ground truth if requested
        if (fUseNodeDeghostingGroundTruth) {
          truth.push_back(nodeDeghostingGroundTruth[spIdx]);
          // Also add direction ground truth
          if (fUseNodeDirectionGroundTruth) {
            truth.insert(truth.end(), nodeDirectionGroundTruth->at(spIdx).begin(),
              nodeDirectionGroundTruth->at(spIdx).end());
          }
        }

        // Add a node with the requested features & ground truth
        newGraph.AddNode(position, features, truth);
      }

      if (fSaveParticleFlow) {
        gpf->push_back(graphUtil.GetParticleFlowMap(trueParticles));
      }

      mf::LogInfo("GCNGraphMaker") << "Produced GCNGraph object with "
        << newGraph.GetNumberOfNodes() << " nodes from " << spacePoints.size() << " spacepoints.";

      // Add out graph to the vector
      graphs->push_back(newGraph);
    }

    // Write our graph to the event
    evt.put(std::move(graphs));
    if (fSaveParticleFlow) evt.put(std::move(gpf));
  }

DEFINE_ART_MODULE(cvn::GCNGraphMaker)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////
