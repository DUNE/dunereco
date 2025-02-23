/**
 *  @file   dunereco/FDSensOpt/ParticleSelectionAlg/ParticleSelectionAlg.h
*
*  @brief  Header file for the Particle Selection
*  Written by Henrique Souza (henrique.souza@mib.infn.it
*
*  $Log: $
*/
#ifndef DUNE_PARTICLE_SELECTION_ALG_H
#define DUNE_PARTICLE_SELECTION_ALG_H

//STL
#include <string>
#include <iostream>
//ART
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#pragma GCC diagnostic ignored "-Wunused-variable"
//LArSoft
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace dune
{
    /**
     *
     * @brief ParticleSelectionAlg class
     *
     */
    class ParticleSelectionAlg
    {
        public:
            /**
             * @brief  Constructor
             * @param  pset the FCL parameter set
             * @param  trackLabel the track label
             * @param  showerLabel the shower label
             * @param  hitLabel the hit label
             * @param  wireLabel the wire label
             * @param  trackToHitLabel the associated track-to-hit label
             * @param  showerToHitLabel the associated shower-to-hit label
             * @param  hitToSpacePointLabel the associated hit-to-space point label
             */
            ParticleSelectionAlg(
                    const fhicl::ParameterSet &pset,
                    const std::string &trackLabel,
                    const std::string &showerLabel,
                    const std::string &hitLabel,
                    const std::string &trackToHitLabel,
                    const std::string &showerToHitLabel,
                    const std::string &hitToSpacePointLabel);

            /**
             * @brief  Retrieve longest tracks from all recob::Tracks of the event
             *
             * @param  event the art event
             * @param  isTrackOnly Set to true as default. Only consider recob::Tracks which recob::PFParticle returns IsTrack as true
             *
             * @return the longest recob::Track
             */
            art::Ptr<recob::Track> GetLongestTrack(const art::Event &event, const bool isTrackOnly = true);

            /**
             * @brief  Retrieve longest tracks of the event using PIDA
             *
             * @param  event the art event
             * @param  tracks vector of recob::Tracks that are considered
             * @param  isTrackOnly Set to true as default. Only consider recob::Tracks which recob::PFParticle returns IsTrack as true
             *
             * @return the longest recob::Track
             */
            art::Ptr<recob::Track> GetLongestTrack(const art::Event &event,
                                                   const std::vector<art::Ptr<recob::Track>> &tracks,
                                                   const bool isTrackOnly);

            /**
             * @brief  Retrieve longest tracks from a vector of recob::Tracks
             * Store it in private member fMuTrack
             *
             * @param  event the art event
             *
             * @return the longest recob::Track
             */
            art::Ptr<recob::Track> GetLongestTrackPID(const art::Event& event);

            /**
             * @brief Apply filters and returns subset of recob::Tracks that
             * are candidates to be muons.
             *
             * @param  event the art event
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> GenMuonCandidates(const art::Event &event);

            std::vector<art::Ptr<recob::Track>> GenProtonCandidates(const art::Event &event);

            std::vector<art::Ptr<recob::Track>> GenPionCandidates(const art::Event &event);


            /**
             * @brief  Retrieve Highest Charge Shower from event
             *
             * @param  clockData 
             * @param  detProp
             * @param  event
             *
             * @return the highest charge recob::Shower
             */
            art::Ptr<recob::Shower> GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                           detinfo::DetectorPropertiesData const& detProp,
                                                           const art::Event &event);

            /**
             * @brief  Retrieve Highest Charge Shower from vector of showers
             *
             * @param  clockData 
             * @param  detProp
             * @param  event
             * @param  showers
             *
             * @return the highest charge recob::Shower
             */
            art::Ptr<recob::Shower> GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                           detinfo::DetectorPropertiesData const& detProp,
                                                           const art::Event &event,
                                                           const std::vector<art::Ptr<recob::Shower>> showers);

            /**
             * @brief Retrieve the PIDA score of the track from map
             *
             * @param  track the track

             * @return PIDA score
             */
            const double GetPIDAScore(const art::Ptr<recob::Track> &track);

            /**
             * @brief Retrieve the Calorimetry energy of track from map
             *
             * @param  track the track

             * @return Calo energy in GeV
             */
            const double GetTrackCalo(const art::Ptr<recob::Track> &track);

            /**
             * @brief Retrieve the momentum by range (proton) of the track from a map
             *
             * @param  track the track

             * @return Momentum in GeV/c
             */
            const double GetTrackMom(const art::Ptr<recob::Track> &track);

            const art::Ptr<recob::Track> &GetLepTrack(){ return fLepTrack; }
            const std::vector<art::Ptr<recob::Track>> &GetPrTracks(){ return fPrTrack; }
            const std::vector<art::Ptr<recob::Track>> &GetPiTracks(){ return fPiTrack; }


        private:
            calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm
            std::string fTrackLabel;                                 ///< the track label
            std::string fShowerLabel;                                ///< the shower label
            std::string fHitLabel;                                   ///< the hit label
            std::string fTrackToHitLabel;                            ///< the associated track-to-hit label
            std::string fShowerToHitLabel;                           ///< the associated shower-to-hit label
            std::string fHitToSpacePointLabel;                       ///< the associated hit-to-space point label

            double fMaxMuLenght;

            const unsigned int kNplanes = 3;

            double fRecombFactor;
            geo::View_t fPlane;

            std::string fParticleIDModuleLabel = "pandorapid";
            art::Handle< std::vector<recob::Track>> trackListHandle;
            std::unique_ptr<art::FindManyP<anab::ParticleID>> fmpid;

            double fDistanceToWallThreshold;                         ///< the min distance from a detector wall to be considered contained

            std::map<art::Ptr<recob::Track>::key_type, double> fTrkIdToPIDAMap;
            std::map<art::Ptr<recob::Track>::key_type, double> fTrkIdToCaloMap;
            std::map<art::Ptr<recob::Track>::key_type, double> fTrkIdToMomMap;
            std::map<art::Ptr<recob::Track>::key_type, bool> fTrkIdToContainmentMap;

            double fTotalCaloEvent = std::numeric_limits<double>::lowest();

            bool fAllTracksContained; 

            art::Ptr<recob::Track> fLepTrack;
            std::vector<art::Ptr<recob::Track>> fPrTrack;
            std::vector<art::Ptr<recob::Track>> fPiTrack;

            bool fPrDone = false; 

            /**
             * @brief Generates a map between track key and PIDA score 
             *
             * @param  event the art event
             * @param  plane plane to be used, options geo::kU, geo::kV,
             * geo::kW, geo::kUnknown. If set to geo::kUnknown, PIDA will be
             * taken from the plane with most valid entries.
             *
             * @return map of track.key and PIDA score
             */
            std::map<art::Ptr<recob::Track>::key_type, double> GenMapPIDAScore(
                    const art::Event &event);

            bool GenContainmentInfo(const art::Event &event);

            bool IsContained(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Event &event);

            bool IsPointContained(const double x, const double y, const double z);

            bool IsTrkContained(const art::Ptr<recob::Track> &track);

            double GetTotalCaloEvent(const art::Event &event);

            std::map<art::Ptr<recob::Track>::key_type, double> GenMapTracksCalo(const art::Event &event);

            std::map<art::Ptr<recob::Track>::key_type, double> GenMapTracksMom(const art::Event &event);

            double CalculateEnergyFromCharge(const double charge, const unsigned short plane );

            bool IsTrack(const art::Event &event, const art::Ptr<recob::Track> &track);

            bool IsTrack(const art::Event &event, const art::Ptr<recob::Shower> shower);

            const void applyMuLengthFilter(
                    std::vector<art::Ptr<recob::Track>>& tracks,
                    const double maxlen);

            const void applyMuMaxPIDA(
                    std::vector<art::Ptr<recob::Track>>& tracks,
                    const double maxPIDA);

            const void applyMuShowerCutCalo(
                    const art::Event &event,
                    std::vector<art::Ptr<recob::Track>>& tracks,
                    const double minCalo,
                    const double maxCalo);

            const void applyMuShowerCutPIDA(
                    const art::Event &event,
                    std::vector<art::Ptr<recob::Track>>& tracks,
                    const double minPIDA);

            const void applyMuCutContained(
                    std::vector<art::Ptr<recob::Track>>& tracks
                    );

            const void applyMuRemoveShowers(
                    const art::Event &event,
                    std::vector<art::Ptr<recob::Track>>& tracks
                    );
    };


}

#endif
