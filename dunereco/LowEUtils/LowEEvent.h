/**
 *  @file LowEEvent.h
 *
 *  @brief LowEEvent class for low energy neutrino studies in DUNE.
 *
 *  @author sergio.manthey@ciemat.es
 *
 */

#ifndef LOWEEVENT_H
#define LOWEEVENT_H

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"

// DUNE includes
#include "dunereco/LowEUtils/LowEUtils.h"

// std includes
#include <vector>
#include <list>
#include <set>
#include <map>

using namespace solar;

namespace lowe {
    /**
    *  @brief LowEEvent is de1fined as a collection of 1 reco object (blip, cluster or pandora) and a list of PDS flash candidates.
    */
    class LowEEvent {
    public:
        LowEEvent(); // Default constructor: initializes an empty LowEEvent with no clusters or flashes.

        // Radiological Constructor: used for radiological studies, where we have a list of Blips and OpFlashes.
        // Missing implementation... Maybe better to just constrcut with 1 blip and a list of flashes?
        // LowEEvent(const std::vector<blipobj::Blip>& Blips, const std::vector<recob::OpFlash>& flashes)
        //     : fBlips(), fFlashes(flashes), fMatchedFlash() {}
        
        // Solar Constructor: used for solar neutrino studies, where we have a collection of LowEClusters and OpFlashes.
        LowEEvent(const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes);
        
        // Supernova Constructor: used for supernova studies, where we have a pandora particle and OpFlashes.
        // Missing implementation...

        LowEEvent(const LowEEvent&);
        LowEEvent& operator=(LowEEvent const&);

        // Destructor
        ~LowEEvent() = default; // Default destructor

        // Initialize the LowEEvent with clusters and flashes
        void initialize(const std::vector<solar::LowECluster>& clusters, const std::vector<recob::OpFlash>& flashes);

        // Getters
        const std::vector<solar::LowECluster>& getClusters() const { return fClusters; } // returns the clusters. 1st is main cluster, rest are adjacent.
        const solar::LowECluster& getMainCluster() const
        { 
            static solar::LowECluster emptyCluster; // static instance to return when fClusters is empty
            return fClusters.empty() ? emptyCluster : fClusters.front(); 
        } // returns the main cluster, if any.        
        const solar::LowECluster* getAdjacentCluster(size_t i) const { return (i < fClusters.size()) ? &fClusters[i] : nullptr; } // returns the i-th adjacent cluster, if any.
        size_t getNumClusters() const { return fClusters.size(); } // returns the number of clusters in the event.
        const std::vector<recob::OpFlash>& getFlashes() const { return fFlashes; }
        const recob::OpFlash* getFlash(size_t i) const { return (i < fFlashes.size()) ? &fFlashes[i] : nullptr; } // returns the i-th flash, if any.
        const recob::OpFlash& getMatchedFlash() const { return fMatchedFlash; } // returns the matched flash, if any.
        size_t getNumFlashes() const { return fFlashes.size(); } // returns the number of flashes in the event.
        bool hasMatchedFlash() const { return fMatchedFlash.InBeamFrame(); } // returns true if there is a matched flash.

        // Setters
        void setClusters(const std::vector<solar::LowECluster>& clusters) { fClusters = clusters; }
        void setFlashes(const std::vector<recob::OpFlash>& flashes) { fFlashes = flashes; }
        void setMatchedFlash(const recob::OpFlash& flash) { fMatchedFlash = flash; } // sets the matched flash, if any.
        void addCluster(const solar::LowECluster& cluster) { fClusters.push_back(cluster); }
        void addFlash(const recob::OpFlash& flash) { fFlashes.push_back(flash); }
        void addMatchedFlash(const recob::OpFlash& flash) { fMatchedFlash = flash; } // sets the matched flash, if any.
        void clear() { fClusters.clear(); fFlashes.clear(); fMatchedFlash = recob::OpFlash(); } // clears the event data.
        // Function to match the clusters with the flashes using LowEMatchingUtils
        int matchHighestPDSFlash(
            const std::vector<art::Ptr<solar::LowECluster>>& SolarClusterVector,
            const std::vector<art::Ptr<recob::OpFlash>>& PDSFlashes,
            art::Ptr<recob::OpFlash>& MatchedFlash,
            const detinfo::DetectorClocksData& clockData,
            const art::Event& evt,
            bool debug = false)
        {
            // Call the matching utility function to find the best match
            LowEUtils Utils{fhicl::ParameterSet()}; // Changed parentheses to braces
            return Utils.MatchPDSFlash(SolarClusterVector, PDSFlashes, MatchedFlash, clockData, evt, debug);
        }
                
            
    private:
        std::vector<solar::LowECluster> fClusters; // collection of clusters in the event
        std::vector<recob::OpFlash> fFlashes; // collection of PDS flash candidates in the event
        recob::OpFlash fMatchedFlash; // matched flash, if any
        struct LowEInteraction;

    }; // class LowEEvent

    // Typedefs for convenience
    typedef std::vector<lowe::LowEEvent> LowEEventVector; // vector of LowEEvent objects
    typedef std::list<lowe::LowEEvent> LowEEventList; // list of LowEEvent objects
    typedef std::set<lowe::LowEEvent> LowEEventSet; // set of LowEEvent objects
    typedef std::map<size_t, lowe::LowEEvent> LowEEventMap; // map of LowEEvent objects indexed by size_t (e.g., event ID)

    // Define a struct for LowEEvent to hold truth information if needed in the future
    struct LowEInteraction
    {
        int InteractionType; // Type of interaction (e.g., charged current, elastic scattering, etc.)
        float NeutrinoEnergy; // Neutrino energy in MeV
        float NeutrinoVertex[3]; // Interaction vertex as an array of 3 floats

        // Add fields for truth information here, e.g., neutrino energy, interaction type, etc.
        // For now, this is just a placeholder.
    }; // struct LowEInteraction

} // namespace lowe
#endif // LOWEEVENT_H
