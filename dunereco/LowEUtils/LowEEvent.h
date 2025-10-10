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

// DUNE includes
#include "dunereco/LowEUtils/LowEUtils.h"

// std includes
#include <vector>

using namespace solar;
using namespace producer;

namespace lowe {
    /**
    *  @brief LowEEvent is defined as a collection of 1 reco object (blip, cluster or pandora) and a list of PDS flash candidates.
    */
    class LowEEvent {
    public:
        LowEEvent();

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
        const std::vector<int>& getClustersIdx() const { return fClustersIdx; } // returns the indices of the clusters in the event as recovered from the original LowECluster collection
        const std::vector<solar::LowECluster>& getClusters() const { return fClusters; } // returns the clusters. 1st is main cluster, rest are adjacent.
        const std::vector<art::Ptr<solar::LowECluster>> getPtrClusters() const 
        { 
            std::vector<art::Ptr<solar::LowECluster>> ptrs;
            for (size_t i = 0; i < fClusters.size(); ++i) {
                ptrs.emplace_back(); // Use default constructor for empty Ptr
            }
            return ptrs; 
        } // returns the clusters as Ptrs. 1st is main cluster, rest are adjacent.
        const solar::LowECluster& getMainCluster() const
        { 
            static solar::LowECluster emptyCluster; // static instance to return when fClusters is empty
            return fClusters.empty() ? emptyCluster : fClusters.front(); 
        } // returns the main cluster, if any.        
        const solar::LowECluster* getAdjacentCluster(size_t i) const { return (i < fClusters.size()) ? &fClusters[i] : nullptr; } // returns the i-th adjacent cluster, if any.
        size_t getNumClusters() const { return fClusters.size(); } // returns the number of clusters in the event.
        const std::vector<recob::OpFlash>& getFlashes() const { return fFlashes; }
        const recob::OpFlash* getFlash(size_t i) const { return (i < fFlashes.size()) ? &fFlashes[i] : nullptr; } // returns the i-th flash, if any.
        const std::vector<art::Ptr<recob::OpFlash>> getPtrFlashes() const 
        { 
            std::vector<art::Ptr<recob::OpFlash>> ptrs;
            for (size_t i = 0; i < fFlashes.size(); ++i) {
                ptrs.emplace_back(); // Use default constructor for empty Ptr
            }
            return ptrs; 
        } // returns the flashes as Ptrs.
        const recob::OpFlash& getMatchedFlash() const { return fMatchedFlash; } // returns the matched flash, if any.
        size_t getNumFlashes() const { return fFlashes.size(); } // returns the number of flashes in the event.
        bool isFlashMatched() const { return fIsFlashMatched; } // returns true if there is a matched flash.
        bool isMatchedFlashValid() const { return fMatchedFlash.InBeamFrame(); } // returns true if there is a matched flash.

        // Setters
        void setClustersIdx(const std::vector<int>& clustersIdx) { fClustersIdx = clustersIdx; } // sets the indices of the clusters in the event as recovered from the original LowECluster collection
        void setClusters(const std::vector<solar::LowECluster>& clusters) { fClusters = clusters; }
        void setPtrClusters(const std::vector<art::Ptr<solar::LowECluster>>& clusters) 
        {
            fClusters.clear();
            for (const auto& cluster : clusters) {
                fClusters.push_back(*cluster);
            }
        } // sets the clusters from a vector of Ptr<solar::LowECluster>
        void setFlashes(const std::vector<recob::OpFlash>& flashes) { fFlashes = flashes; }
        void setPtrFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes) 
        {
            fFlashes.clear();
            for (const auto& flash : flashes) {
                fFlashes.push_back(*flash);
            }
        } // sets the flashes from a vector of Ptr<recob::OpFlash>
        void setMatchedFlash(const recob::OpFlash& flash) 
        { 
            fMatchedFlash = flash;
            fIsFlashMatched = true; // Set the flag indicating a matched flash exists
        } // sets the matched flash, if any.
        void setPtrMatchedFlash(const art::Ptr<recob::OpFlash>& flash) 
        { 
            if (flash) {
                fMatchedFlash = *flash;
                fIsFlashMatched = true; // Set the flag indicating a matched flash exists 
            } else {
                fMatchedFlash = recob::OpFlash(); // Reset to default if null Ptr
            }
        } // sets the matched flash from a Ptr<recob::OpFlash>
        void addCluster(const solar::LowECluster& cluster) { fClusters.push_back(cluster); }
        void addFlash(const recob::OpFlash& flash) { fFlashes.push_back(flash); }
        void clear() { fClusters.clear(); fFlashes.clear(); fMatchedFlash = recob::OpFlash(); fIsFlashMatched = false; } // clears the event
                    
                
    private:
        // fhicl::ParameterSet pset; // Create a parameter set (you may need to initialize it appropriately)
        std::vector<int> fClustersIdx; // indices of the clusters in the event as recovered from the original Blip/LowECluster/Pandora collection
        std::vector<solar::LowECluster> fClusters; // collection of clusters in the event
        std::vector<recob::OpFlash> fFlashes; // collection of PDS flash candidates in the event
        recob::OpFlash fMatchedFlash; // matched flash, if any
        bool fIsFlashMatched; // flag to indicate if a matched flash exists
    }; // class LowEEvent

    // Typedefs for convenience
    typedef std::vector<lowe::LowEEvent> LowEEventVector; // vector of LowEEvent objects
    typedef std::list<lowe::LowEEvent> LowEEventList; // list of LowEEvent objects
    typedef std::set<lowe::LowEEvent> LowEEventSet; // set of LowEEvent objects
    typedef std::map<size_t, lowe::LowEEvent> LowEEventMap; // map of LowEEvent objects indexed by size_t (e.g., event ID)

} // namespace lowe
#endif // LOWEEVENT_H
