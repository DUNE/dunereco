/////////////////////////////////////////////////////////////////////////////////////
// Class:        LowEEvent                                                         //
// Module Type:  producer                                                          //
// Module Label: solarevent                                                        //
// File:         LowEEvent_module.cc                                               //
//                                                                                 //
// Event producer for LowE analysis. It collects information about                 //
// particles, blips/clusters/pandora and flashes.                                  //
//                                                                                 //
// Can be used as template for other LowEEvent producers (radiological, supernova) //
/////////////////////////////////////////////////////////////////////////////////////

#ifndef SolarEvent_H
#define SolarEvent_H 1

#include "dunereco/LowEUtils/LowEEvent.h"

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
        void reconfigure(fhicl::ParameterSet const &p);
    
        // The producer routine, called once per event.
        void produce(art::Event &);
    
    private:
        // The parameters we'll read from the .fcl file.
        std::string fSignalLabel;   // Input tag for MCTruth label
        std::string fGEANTLabel;    // Input tag for GEANT label
        std::string fClusterLabel;  // Input tag for LowECluster collection
        std::string fOpFlashLabel;  // Input tag for OpFlash collection
        std::string fFlashAlgoType; //Choose the type of algorithm for the Flashmatch
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

#endif // SolarEvent_H

//------------------------------------------------------------------------------

namespace lowe
{
    //--------------------------------------------------------------------------
    // Constructor
    SolarEvent::SolarEvent(fhicl::ParameterSet const &p)
        : EDProducer{p},
        fSignalLabel(p.get<std::string>("SignalLabel", "marley")),
        fGEANTLabel(p.get<std::string>("GEANTLabel", "largeant")),
        fClusterLabel(p.get<std::string>("ClusterLabel", "solarcluster")),
        fOpFlashLabel(p.get<std::string>("OpFlashLabel", "solarflash")),
        fFlashAlgoType(p.get<std::string>("FlashAlgoType")),  
        fDebug(p.get<bool>("Debug", false)),      
        lowe(new lowe::LowEUtils(p)),
        producer(new producer::ProducerUtils(p))
    {
        reconfigure(p);
        produces<std::vector<lowe::LowEEvent>>();
        produces<art::Assns<lowe::LowEEvent, solar::LowECluster>>("LowEClusterAssns");
        produces<art::Assns<lowe::LowEEvent, recob::OpFlash>>("LowEFlashAssns");
    }

    //--------------------------------------------------------------------------
    void SolarEvent::reconfigure(fhicl::ParameterSet const &p)
    {
        // Reconfigure the parameters from the FHiCL file
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
        std::string sSolarEvent = "SolarEvent::produce called\n";
        std::string sSolarEventFinding = "";
        std::string sSolarEventErrors = "";
        // Initialize the backtracking variables
        std::vector<std::set<int>> trackids = {};
        std::map<int, simb::MCParticle> ThisGeneratorParts;
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
        std::vector<art::Ptr<recob::OpFlash>> FlashPtr;

        // Get MCTruth information
        auto SignalHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(fSignalLabel); // Get generator handle
        std::string SolarTruthInfo = "SolarEvent::produce: MCTruth information for label " + fSignalLabel;
        art::FindManyP<simb::MCParticle> Assn(SignalHandle, evt, fGEANTLabel);            // Assign labels to MCPArticles
        producer->FillMyMaps(ThisGeneratorParts, Assn, SignalHandle);                  // Fill empty list with previously assigned particles
        SolarTruthInfo = SolarTruthInfo + "\n# of particles " + ProducerUtils::str(int(ThisGeneratorParts.size())) + " for " + fSignalLabel;

        const simb::MCNeutrino &nue = SignalHandle->at(0).GetNeutrino();
        SolarTruthInfo = SolarTruthInfo + "\n\tNeutrino Energy: " + ProducerUtils::str(1e3*nue.Nu().E()) + " MeV with vertex (x,y,z) " + ProducerUtils::str(nue.Nu().Vx()) + ", " + ProducerUtils::str(nue.Nu().Vy()) + ", " + ProducerUtils::str(nue.Nu().Vz());

        float eKE = 0;
        float eVertex[3] = {-1e6, -1e6, -1e6}; // Initialize electron vertex to invalid values
        if (ThisGeneratorParts.size() > 0) {
            std::string ElectronTruthInfo;
            for (std::map<int, simb::MCParticle>::iterator iter = ThisGeneratorParts.begin(); iter != ThisGeneratorParts.end(); iter++)
            {
                std::set<int> ThisGeneratorIDs = {};
                trackids.push_back(ThisGeneratorIDs);
                trackids[0].insert(iter->first);
                int pdg = iter->second.PdgCode();

                if (pdg == 11 && 1e3*iter->second.E()-1e3*iter->second.Mass() > eKE) {
                    eKE = 1e3*iter->second.E()-1e3*iter->second.Mass();
                    eVertex[0] = iter->second.EndX();
                    eVertex[1] = iter->second.EndY();
                    eVertex[2] = iter->second.EndZ();
                    ElectronTruthInfo = "\n\tMain e⁻ Energy: " + ProducerUtils::str(eKE) + "MeV with vertex (x,y,z) " + ProducerUtils::str(eVertex[0]) + ", " + ProducerUtils::str(eVertex[1]) + ", " + ProducerUtils::str(eVertex[2]);
                }
            }

            if (eKE > 0) {
                SolarTruthInfo = SolarTruthInfo + ElectronTruthInfo;
            }
        }

        else {
            std::set<int> ThisGeneratorIDs = {};
            trackids.push_back(ThisGeneratorIDs);
        }

        producer->PrintInColor(SolarTruthInfo, producer::ProducerUtils::GetColor("magenta"));

        auto ClusterHandle = evt.getValidHandle<std::vector<solar::LowECluster>>(fClusterLabel); // Get LowECluster handle
        auto FlashHandle = evt.getValidHandle<std::vector<recob::OpFlash>>(fOpFlashLabel); // Get OpFlash handle
        for (const auto &cluster : *ClusterHandle) {
            ClusterPtr.push_back(art::Ptr<solar::LowECluster>(ClusterHandle, &cluster - &(*ClusterHandle)[0]));
        }
        for (const auto &flash : *FlashHandle) {
            FlashPtr.push_back(art::Ptr<recob::OpFlash>(FlashHandle, &flash - &(*FlashHandle)[0]));
        }

        // Process the clusters to generate groups of Main and Adjacent clusters
        lowe->FindPrimaryClusters(ClusterPtr, EventCandidateFound, EventCandidateVector, EventCandidateIdx, clockData, evt);
        
        if (EventCandidateVector.empty()) {
            sSolarEventErrors += "No event candidates found.\n";
        }

        sSolarEventFinding += "Found " + ProducerUtils::str(EventCandidateVector.size()) + " event candidates.\n";

        // Fill the LowEEvent object with the found clusters
        for (size_t i = 0; i < EventCandidateVector.size(); ++i)
        {
            sSolarEventFinding += "Processing event candidate " + ProducerUtils::str(i) + " with " + ProducerUtils::str(EventCandidateVector[i].size()) + " clusters.\n";
            lowe::LowEEvent event;                 
            event.setClustersIdx(EventCandidateIdx[i]); // Set the indices of the clusters in the event
            event.setPtrClusters(EventCandidateVector[i]); // Set the clusters from the vector of Ptr<solar::LowECluster>
            event.setPtrFlashes(FlashPtr); // Set the flashes from the vector of Ptr<recob::OpFlash>

            if (fFlashAlgoType=="SolarFlashMatch") {
               
               int matchedFlashIndex = -1;
            
               if (!FlashPtr.empty()) {
                   matchedFlashIndex = lowe->MatchPDSFlash(EventCandidateVector[i], FlashPtr, clockData, evt, fDebug);
                
                   if (matchedFlashIndex >= 0) {
                       art::Ptr<recob::OpFlash> matchedFlashPtr = FlashPtr[matchedFlashIndex];
                       event.setPtrMatchedFlash(matchedFlashPtr); // Set the matched flash from a Ptr<recob::OpFlash>
                       sSolarEventFinding += "Matched flash index: " + ProducerUtils::str(matchedFlashIndex) + "\n";
                    }
                
                   else {
                       sSolarEventFinding += "No flash candidate matched to the cluster.\n";
                    }
                }

               else {
                   sSolarEventFinding += "No OpFlashes available for matching.\n";
               }

               // Check that method isMatchedFlashValid() returns true
               if (event.isMatchedFlashValid()) {
                   sSolarEventFinding += "Matched flash is valid.\n";
                }

               else {
                   sSolarEventErrors += "Matched flash is not valid.\n";
                }
            }

            else if (fFlashAlgoType == "likelihoodFlashMatch") {
                std::cout << "likelihoodFlashmatch!!! "<< std ::endl;
            }
            
            // Add the event to the SolarEventPtr vector
            SolarEventPtr->push_back(event);
        }
        
        // Create associations between LowEEvents and LowEClusters
        for (size_t i = 0; i < SolarEventPtr->size(); ++i)
        {
            const auto &event = SolarEventPtr->at(i);
            std::vector<art::Ptr<solar::LowECluster>> clusterPtr = {};
            for (const auto &clusterIdx : event.getClustersIdx())
            {
                if (clusterIdx < int(ClusterPtr.size())) {
                    clusterPtr.push_back(ClusterPtr[clusterIdx]);
                }

                else {
                    // producer->PrintInColor("Cluster index out of bounds: " + ProducerUtils::str(clusterIdx), ProducerUtils::GetColor("red"), "Debug");
                    sSolarEventErrors += "Cluster index out of bounds: " + ProducerUtils::str(clusterIdx) + "\n";
                }
            }

            util::CreateAssn(*this, evt, *SolarEventPtr, clusterPtr, *(ClusterAssnPtr.get()), i);
            util::CreateAssn(*this, evt, *SolarEventPtr, FlashPtr, *(FlashAssnPtr.get()), i);
        }

        // Store the LowEEvents and associations in the event
        evt.put(std::move(SolarEventPtr));
        evt.put(std::move(ClusterAssnPtr), "LowEClusterAssns");
        evt.put(std::move(FlashAssnPtr), "LowEFlashAssns");

        // Print the debug information
        if (fDebug) {
            producer->PrintInColor(sSolarEvent, producer::ProducerUtils::GetColor("yellow"));
            producer->PrintInColor(sSolarEventFinding, producer::ProducerUtils::GetColor("green"), "Debug");
            producer->PrintInColor(sSolarEventErrors, producer::ProducerUtils::GetColor("red"), "Debug");
        }
    }
} // namespace lowe
