/////////////////////////////////////////////////////////////////////////////////////
// Class:       LowEEvent                                                          //
// Module Type: producer                                                           //
// Module Label: solarevent                                                        //
// File:        LowEEvent_module.cc                                                //
//                                                                                 //
// Event producer for LowE analysis. It collects information about                 //
// particles, blips/clusters/pandora and flashes.                                  //
//                                                                                 //
// Can be used as template for other LowEEvent producers (radiological, supernova) //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef LowEEventProducer_H
#define LowEEventProducer_H 1

#include "dunereco/LowEUtils/LowEUtils.h"
#include "dunereco/LowEUtils/LowEEvent.h"
#include "dunereco/LowEUtils/LowECluster.h"
#include "dunecore/ProducerUtils/ProducerUtils.h"

using namespace solar;
using namespace producer;

namespace lowe
{
    class SolarEvent : public art::EDProducer
    {
    public:
        explicit SolarEvent(const fhicl::ParameterSet &);
        virtual ~SolarEvent();
    
        void beginJob();
        void endJob();
        void reconfigure(fhicl::ParameterSet const &pset);
    
        // The producer routine, called once per event.
        void produce(art::Event &);
    
    private:
        // The parameters we'll read from the .fcl file.
        std::string fSignalLabel; // Input tag for MCTruth label
        std::string fGEANTLabel; // Input tag for GEANT label
        std::string fHitLabel; // Input tag for Hit collection
        std::string fClusterLabel; // Input tag for LowECluster collection
        std::string fOpFlashLabel; // Input tag for OpFlash collection
        bool fDebug;
        lowe::LowEUtils *lowe;
        producer::ProducerUtils *producer;
    };
    } // namespace lowe

//------------------------------------------------------------------------------

namespace lowe
{
    DEFINE_ART_MODULE(SolarEvent)
}

#endif // LowEEventProducer_H

//------------------------------------------------------------------------------

namespace lowe
{
    //--------------------------------------------------------------------------
    // Constructor
    SolarEvent::SolarEvent(const fhicl::ParameterSet &p)
        : EDProducer{p},
          fSignalLabel(p.get<std::string>("SignalLabel")),
        //   fGEANTLabel(p.get<std::string>("GEANT4Label")),
          fClusterLabel(p.get<std::string>("ClusterLabel")),
          fOpFlashLabel(p.get<std::string>("OpFlashLabel")),
          fDebug(p.get<bool>("Debug")),
          producer(new producer::ProducerUtils(p))
    {
        reconfigure(p);
        produces<std::vector<solar::LowECluster>>();
    }

    //--------------------------------------------------------------------------
    void SolarEvent::reconfigure(fhicl::ParameterSet const &p)
    {
        // Reconfigure parameters if needed
    }

    //--------------------------------------------------------------------------
    // Destructor
    SolarEvent::~SolarEvent()
    {
        delete producer;
    }
    
    //--------------------------------------------------------------------------
    void SolarEvent::beginJob()
    {
        // Initialization code if needed
    }
    
    //--------------------------------------------------------------------------
    void SolarEvent::endJob()
    {
        // Cleanup code if needed
    }
    
    //--------------------------------------------------------------------------
    void SolarEvent::produce(art::Event &evt)
    {
        // Initialize the objects for LowEUtils and ProducerUtils instances
        std::vector<bool> EventCandidateFound = {}; // Vector to store if a cluster is a candidate for an event
        std::vector<std::vector<int>> EventCandidateIdx = {}; // Indices of Main and Adjacent clusters
        std::vector<std::vector<art::Ptr<solar::LowECluster>>> EventCandidateVector = {}; // Main and Adjacent clusters
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
        
        // Prepare the output data: a vector of LowEEvents and associations to LowEClusters and OpFlashes
        std::unique_ptr<std::vector<lowe::LowEEvent>> SolarEventPtr(new std::vector<lowe::LowEEvent>);
        std::unique_ptr<art::Assns<lowe::LowEEvent, solar::LowECluster>> ClusterAssnPtr(new art::Assns<lowe::LowEEvent, solar::LowECluster>);
        std::unique_ptr<art::Assns<lowe::LowEEvent, recob::OpFlash>> FlashAssnPtr(new art::Assns<lowe::LowEEvent, recob::OpFlash>);

        // Get input data: LowEClusters and OpFlashes (need to implement blips/clusters/pandora handling)
        std::vector<art::Ptr<solar::LowECluster>> ClusterPtr;
        art::Handle<std::vector<solar::LowECluster>> ClusterHandle;
        std::vector<art::Ptr<recob::OpFlash>> FlashPtr;
        art::Handle<std::vector<recob::OpFlash>> FlashHandle;

        if (evt.getByLabel(fClusterLabel, ClusterHandle))
        {
            art::fill_ptr_vector(ClusterPtr, ClusterHandle);
            producer->PrintInColor("Found valid handle for LowEClusters!", ProducerUtils::GetColor("green"));
        }
        else
        {
            producer->PrintInColor("Valid handle not found for LowEClusters!", ProducerUtils::GetColor("red"), "Debug");
        }
        if (evt.getByLabel(fHitLabel, FlashHandle))
        {
            art::fill_ptr_vector(FlashPtr, FlashHandle);
            producer->PrintInColor("Found valid handle for OpFlashes!", ProducerUtils::GetColor("green"));
        }
        else
        {
            producer->PrintInColor("Valid handle not found for OpFlashes!", ProducerUtils::GetColor("red"), "Debug");
        }

        // Process the clusters to generate groups of Main and Adjacent clusters
        lowe->FindPrimaryClusters(ClusterPtr, EventCandidateFound, EventCandidateVector, EventCandidateIdx, clockData, evt);

        // Fill the LowEEvent object with the found clusters
        for (size_t i = 0; i < EventCandidateVector.size(); ++i)
        {
            if (EventCandidateFound[i])
            {
                lowe::LowEEvent event;
                event.setClustersIdx(EventCandidateIdx[i]); // Set the indices of the clusters in the event
                event.setPtrClusters(EventCandidateVector[i]);
                event.setPtrFlashes(FlashPtr);
                SolarEventPtr->push_back(event);
            }
        }
        
        // Create associations between LowEEvents and LowEClusters
        for (size_t i = 0; i < SolarEventPtr->size(); ++i)
        {
            const auto &event = SolarEventPtr->at(i);
            // Create associations with LowEClusters by retrieving the cluster indices and accessing the clusters
            art::PtrVector<solar::LowECluster> clusterPtr;
            for (const auto &clusterIdx : event.getClustersIdx())
            {
                if (clusterIdx < int(ClusterPtr.size()))
                {
                    clusterPtr.push_back(ClusterPtr[clusterIdx]);
                }
                else
                {
                    producer->PrintInColor("Cluster index out of bounds: " + std::to_string(clusterIdx), ProducerUtils::GetColor("red"), "Debug");

                }
            }
            util::CreateAssn(*this, evt, *SolarEventPtr, clusterPtr, *(ClusterAssnPtr.get()), i);

            // Create associations with OpFlashes
            // util::CreateAssn(*this, evt, *SolarEventPtr, event.getPtrFlashes(), *(FlashAssnPtr.get()), i);
        }

        // Store the LowEEvents and associations in the event
        evt.put(std::move(SolarEventPtr));
        evt.put(std::move(ClusterAssnPtr), "LowEClusterAssns");
        // evt.put(std::move(FlashAssnPtr), "LowEFlashAssns");
    }
} // namespace lowe

