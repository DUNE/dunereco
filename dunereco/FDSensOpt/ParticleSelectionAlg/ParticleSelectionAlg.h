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
             * @param  pfpLabel the track label
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
                    const std::string &pfpLabel,
                    const std::string &trackLabel,
                    const std::string &showerLabel,
                    const std::string &hitLabel,
                    const std::string &trackToHitLabel,
                    const std::string &showerToHitLabel,
                    const std::string &hitToSpacePointLabel);


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
             * using PIDA
             *
             * @param  event the art event
             *
             * @return the longest recob::Track
             */
            art::Ptr<recob::Track> GetLongestTrackPID(const art::Event& event);


            /**
             * @brief  Default muon selection using filters
             *
             * @param  event the art event
             * @param  candidates initial tracks list to select
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> MuDefaultSelection(const art::Event &event, std::vector<art::Ptr<recob::Track>> &candidates);

            /**
             * @brief  Apply filters and returns subset of recob::Tracks that
             * are candidates to be muons.
             *
             * @param  event the art event
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> GenMuonCandidates(const art::Event &event);

            /**
             * @brief  Default electron selection using filters
             *
             * @param  event the art event
             * @param  candidates initial shower list to select
             *
             * @return vector of recob::Shower
             */
            std::vector<art::Ptr<recob::Shower>> EDefaultSelection(const art::Event &event, std::vector<art::Ptr<recob::Shower>> &candidates);

            /**
             * @brief  Apply filters and returns subset of recob::Showers that
             * are candidates to be electrons.
             *
             * @param  event the art event
             *
             * @return vector of recob::Shower
             */
            std::vector<art::Ptr<recob::Shower>> GenElectronCandidates(const art::Event &event);

            /**
             * @brief  Default proton selection
             *
             * @param  event the art event
             * @param  tracks initial tracks list to select
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> PrDefaultSelection(const art::Event &event, std::vector<art::Ptr<recob::Track>> &tracks);

            /**
             * @brief  Apply filters and returns subset of recob::Tracks that
             * are candidates to be protons.
             *
             * @param  event the art event
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> GenProtonCandidates(const art::Event &event);

            /**
             * @brief  Default pion selection
             *
             * @param  event the art event
             * @param  tracks initial tracks list to select
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> PionDefaultSelection(const art::Event &event, std::vector<art::Ptr<recob::Track>> &tracks);

            /**
             * @brief  Apply filters and returns subset of recob::Tracks that
             * are candidates to be pions.
             *
             * @param  event the art event
             *
             * @return vector of recob::Track
             */
            std::vector<art::Ptr<recob::Track>> GenPionCandidates(const art::Event &event);

            /**
             * @brief  Retrieve Highest Charge Shower from vector of showers.
             * The shower is also saved in private member fLepShower and the 
             * associated track (if available) is saved in fLepTrack
             *
             * @param  clockData 
             * @param  detProp
             * @param  event
             * @param  showers
             * @param  tPlane reference plane for getting charge
             *
             * @return the highest charge recob::Shower
             */
            art::Ptr<recob::Shower> GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                           detinfo::DetectorPropertiesData const& detProp,
                                                           const art::Event &event,
                                                           const std::vector<art::Ptr<recob::Shower>> &showers,
                                                           const geo::View_t tPlane = geo::kW);


            /**
             * @brief  Retrieve Highest Charge Shower from event using PIDA
             *
             * @param  event
             *
             * @return the highest charge recob::Shower
             */
            art::Ptr<recob::Shower> GetHighestChargeShowerPID(const art::Event &event);

            /**
             * @brief  Retrieve the PIDA score of the track from map
             *
             * @param  track the track

             * @return PIDA score
             */
            const double GetPIDAScore(const art::Ptr<recob::Track> &track);

            /**
             * @brief  Retrieve the PIDA score of the shower from map
             *
             * @param  shower the shower
             *
             * @return PIDA score
             */

            const double GetPIDAScore(const art::Ptr<recob::Shower> &shower);
            /**
             * @brief  Retrieve the Calorimetry energy of track from map
             *
             * @param  track the track

             * @return Calo energy in GeV
             */
            const double GetTrackCalo(const art::Ptr<recob::Track> &track);

            /**
             * @brief  Retrieve the momentum by range (proton) of the track from a map
             *
             * @param  track the track

             * @return Momentum in GeV/c
             */
            const double GetTrackMom(const art::Ptr<recob::Track> &track);

            const art::Ptr<recob::Track> &GetLepTrack(){ return fLepTrack; }
            const art::Ptr<recob::Shower> &GetLepShower(){ return fLepShower; }
            const std::vector<art::Ptr<recob::Track>> &GetPrTracks(){ return fPrTracks; }
            const std::vector<art::Ptr<recob::Track>> &GetPiTracks(){ return fPiTracks; }

            void SetLepTrack(const art::Ptr<recob::Track> &trk){ fLepTrack = trk; }
            void SetPrTracks(const std::vector<art::Ptr<recob::Track>> &tracks){ fPrTracks = tracks; }
            void SetPiTracks(const std::vector<art::Ptr<recob::Track>> &tracks){ fPiTracks = tracks; }

            /**
             * @brief  Retrieve the track assciated to the shower
             *
             * @param  event the art event
             * @param  shower the shower
             *
             * @return the track
             */
            art::Ptr<recob::Track> GetTrackFromShower(const art::Event &event, const art::Ptr<recob::Shower> &shower);

        private:
            calo::CalorimetryAlg fCalorimetryAlg;                    
            int    fMuSelectMethod;                                  ///< the method for muon selection
            double kMaxMuLength;                                     ///< the max length for muon candidates
            double kMaxMuPIDA;                                       ///< the max PIDA for muon candidates
            double kMinMuTotalCalo;                                  ///< the min calorimetric energy for muon candidates when they are showers
            double kMaxMuTotalCalo;                                  ///< the max calorimetric energy for muon candidates when they are showers
            double kMinMuPIDAShower;                                 ///< the min PIDA for muon candidates when they are showers
                                                                     ///(in combination with calorimetric energy)
            double kMaxMuPIDAAggressive;                             ///< the max PIDA for muon candidates (more aggressive, after all cleanings)
            double kMaxMuContainedCalo;                              ///< the max calorimetric energy for muon candidates.
                                                                     ///Used to exclude tracks that are contained


            int    fESelectMethod;                                  ///< the method for electron selection
            double kMaxETotalCaloForTracks;                          ///< the max calorimetric energy for electron candidates as track
            double kMaxEPIDA;                                        ///< the max PIDA for electron candidates as track

            int    fPrSelectMethod;                                  ///< the method for proton selection
            double kMinPrPIDATrack;                                  ///< the min PIDA for proton candidates as track
            double kMinPrPIDAShower;                                 ///< the min PIDA for proton candidates as shower
            double kMaxPrTrkCalo;                                    ///< the max calorimetric energy for proton candidates
            double kMaxPrTrkMom;                                     ///< the max momentum for proton candidates

            double kPrMomByRangeMinLength;                           ///< the min length for momemtum by range method
            double kPrMomByRangeMaxLength;                           ///< the max length for momemtum by range method

            int    fPionSelectMethod;                                ///< the method for pion selection
            size_t kMinPionNTrk;                                     ///< the min number of tracks for searching pion candidates
            double kMinPionTrkLength;                                ///< the min length for pion candidates
            double kMaxPionTrkLength;                                ///< the max length for pion candidates
            double kMinPionTrkLengthNotContained;                    ///< the min length for pion candidates not contained

            double kDistanceToWallThreshold;                         ///< the min distance from a detector wall to be considered contained

            double kRecombFactor;                                    ///< the recombination factor
            unsigned int kPlaneToUse;                                ///< the plane to use
                                                                     ///Set to geo::kUnknown to use the plane with most valid entries

            std::string fPFParticleLabel;                            ///< the pfp label
            std::string fTrackLabel;                                 ///< the track label
            std::string fShowerLabel;                                ///< the shower label
            std::string fHitLabel;                                   ///< the hit label
            std::string fTrackToHitLabel;                            ///< the associated track-to-hit label
            std::string fShowerToHitLabel;                           ///< the associated shower-to-hit label
            std::string fHitToSpacePointLabel;                       ///< the associated hit-to-space point label
            std::string fParticleIDModuleLabel;                      ///< the associated particle id label


            const unsigned int kNplanes = 3;
            geo::View_t kPlane; 


            // The following are the default values for the ParticleID module label
            // The code is much faster if FindManyP is used this way
            art::Handle< std::vector<recob::Hit> > hitListHandle; ///< handle to the hit list
            std::unique_ptr<art::FindManyP<recob::SpacePoint>> fmsp; ///< pointer to the space point


            std::map<art::Ptr<recob::Track>::key_type, double> fTrkIdToPIDAMap;            ///< map of track.key and PIDA score
            std::map<art::Ptr<recob::Track>::key_type, double> fTrkIdToCaloMap;            ///< map of track.key and calorimetric energy
            std::map<art::Ptr<recob::Track>::key_type, double> fTrkIdToMomMap;             ///< map of track.key and momentum
            std::map<art::Ptr<recob::Track>::key_type, bool> fTrkIdToContainmentMap;       ///< map of track.key and containment

            std::map<art::Ptr<recob::Shower>::key_type, double> fShwIdToPIDAMap;           ///< map of shower.key and PIDA score
            std::map<art::Ptr<recob::Shower>::key_type, geo::View_t> fShwIdToBestPlaneMap; ///< map of shower.key and best plane

            double fTotalCaloEvent = std::numeric_limits<double>::lowest();                ///< total calorimetric energy of the event

            art::Ptr<recob::Track> fLepTrack;             ///< the lepton track
            art::Ptr<recob::Shower> fLepShower;           ///< the lepton shower
            std::vector<art::Ptr<recob::Track>> fPrTracks; ///< proton candidates
            std::vector<art::Ptr<recob::Track>> fPiTracks; ///< pion candidates

            bool fPrDone = false;  ///< flag to check if proton candidates are already set


            // Setup of methods, for now, only default method implemented.
            enum class MuSelectMethod
            {
                kSimple = 0,
                kDefault
            };
            enum class ESelectMethod
            {
                kSimple = 0,
                kDefault
            };
            enum class PrSelectMethod
            {
                kDefault = 0,
            };
            enum class PionSelectMethod
            {
                kDefault = 0,
            };

            /**
             * @brief  Generates a map between track key and PIDA score 
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

            /**
             * @brief  Generates a map between shower key and PIDA score
             *
             * @param  event the art event
             *
             * @return map of shower.key and PIDA score
             */
            std::map<art::Ptr<recob::Shower>::key_type, double> GenMapPIDAScoreShw(const art::Event &event);

            /**
             * @brief  Generates a map between shower key and best plane. The
             * best plane is the one with most hits when tPlane is set to
             * geo::kUnknown, otherwise is fixed to tPlane
             *
             * @param  event the art event
             * @param  tPlane reference plane for getting charge
             *
             * @return map of shower.key and best plane
             */
            std::map<art::Ptr<recob::Shower>::key_type, geo::View_t> GenMapBestPlaneShower(const art::Event &event, const geo::View_t tPlane);

            /**
             * @brief  Generates a map between all tracks key and contaiment of SpacePoints
             * Returns true if all hits are contained, false otherwise
             *
             * @param  event the art event
             *
             * @return true if all hits are contained, false otherwise
             */
            bool GenContainmentInfo(const art::Event &event);

            /**
             * @brief  Generates a map between tracks key and contaiment of SpacePoints
             * Returns true if all hits are contained, false otherwise
             *
             * @param  event the art event
             * @param  tracks vector of recob::Tracks to be considered
             *
             * @return true if all hits are contained, false otherwise
             */
            bool GenContainmentInfo(const art::Event &event, const std::vector<art::Ptr<recob::Track>> &tracks);

            /**
             * @brief  Check if list of hits are contained
             *
             * @param  hits the list of hits
             * @param  event the art event
             *
             * @return true if all hits are contained, false otherwise
             */
            bool IsContained(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Event &event);

            /**
             * @brief  Check if point is contained
             *
             * @param  x the x coordinate
             * @param  y the y coordinate
             * @param  z the z coordinate
             *
             * @return true if point is contained, false otherwise
             */
            bool IsPointContained(const double x, const double y, const double z);

            /**
             * @brief  Check if track is contained. Information is retrieved
             * from map fTrkIdToContainmentMap
             *
             * @param  track the track
             *
             * @return true if track is contained, false otherwise
             */
            bool IsTrkContained(const art::Ptr<recob::Track> &track);

            /**
             * @brief  Calculate the total calorimetric energy of the event
             *
             * @param  event the art event
             *
             * @return the total calorimetric energy of the event
             */
            double GetTotalCaloEvent(const art::Event &event);

            /**
             * @brief  Generate a map between track key and their calorimetric energy.
             * The energy is computed using only charge information from hits.
             * The plane used is either the set by kPlaneToUse or the one with
             * most valid entries (when set to geo::kUnknown)
             *
             * @param  event the art event
             *
             * @return map of track.key and their calorimetric energy
             */
            std::map<art::Ptr<recob::Track>::key_type, double> GenMapTracksCalo(const art::Event &event);

            /**
             * @brief  Generate a map between track key and their momentum by range assuming they are protons
             *
             * @param  event the art event
             *
             * @return map of track.key and their momentum by range
             */
            std::map<art::Ptr<recob::Track>::key_type, double> GenMapTracksMom(const art::Event &event);

            /**
             * @brief  Calculate energy from charge for a given plaen
             *
             * @param  charge the charge
             * @param  plane the plane
             *
             * @return the energy in GeV
             */
            double CalculateEnergyFromCharge(const double charge, const unsigned short plane );

            /**
             * @brief  Check PFParticle associated with track is considered as Track by LArPandora
             *
             * @param  event the art event
             * @param  track the track
             *
             * @return true if PFP is considered as Track, false otherwise
             */
            bool IsTrack(const art::Event &event, const art::Ptr<recob::Track> &track);

            /**
             * @brief  Check PFParticle associated with shower is considered as Track by LArPandora
             *
             * @param  event the art event
             * @param  shower the shower
             *
             * @return true if PFP is considered as Shower, false otherwise
             */
            bool IsTrack(const art::Event &event, const art::Ptr<recob::Shower> &shower);


            geo::View_t GetBestPlane(const art::Event &event,  const std::vector<art::Ptr<recob::Hit>> &allhits);


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
                    const art::Event &event,
                    std::vector<art::Ptr<recob::Track>>& tracks
                    );

            const void applyMuRemoveShowers(
                    const art::Event &event,
                    std::vector<art::Ptr<recob::Track>>& tracks
                    );
    };


}

#endif
