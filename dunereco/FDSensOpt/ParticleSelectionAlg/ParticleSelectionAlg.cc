/**
 *  @file   dunereco/FDSensOpt/ParticleSelectionAlg/ParticleSelectionAlg.h
 *
 *  @brief  Implementation file for the Particle Selection algorithm.
 *  Written by Henrique Souza (henrique.souza@mib.infn.it
 *
 *  $Log: $
 */

//STL
#include <limits>
//ART
#include "canvas/Utilities/Exception.h"
//LArSoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

//DUNE
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

#include "dunereco/FDSensOpt/ParticleSelectionAlg/ParticleSelectionAlg.h"

namespace dune
{
    ParticleSelectionAlg::ParticleSelectionAlg(
            fhicl::ParameterSet const& pset,
            const std::string &pfpLabel,
            const std::string &trackLabel,
            const std::string &showerLabel,
            const std::string &hitLabel,
            const std::string &trackToHitLabel,
            const std::string &showerToHitLabel,
            const std::string &hitToSpacePointLabel) :
        fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
        kMaxMuLength(pset.get<double>("MaxMuLength")),
        kMaxMuPIDA(pset.get<double>("MaxMuPIDA")),
        kMinMuTotalCalo(pset.get<double>("MinMuTotalCalo")),
        kMaxMuTotalCalo(pset.get<double>("MaxMuTotalCalo")),
        kMinMuPIDAShower(pset.get<double>("MinMuPIDAShower")),
        kMaxMuPIDAAggressive(pset.get<double>("MaxMuPIDAAggressive")),
        kMaxMuContainedCalo(pset.get<double>("MaxMuContainedCalo")),
        kMaxETotalCaloForTracks(pset.get<double>("MaxETotalCaloForTracks")),
        kMaxEPIDA(pset.get<double>("MaxEPIDA")),
        kMinPrPIDATrack(pset.get<double>("MinPrPIDATrack")),
        kMinPrPIDAShower(pset.get<double>("MinPrPIDAShower")),
        kMaxPrTrkCalo(pset.get<double>("MaxPrTrkCalo")),
        kMaxPrTrkMom(pset.get<double>("MaxPrTrkMom")),
        kPrMomByRangeMinLength(pset.get<double>("PrMomByRangeMinLength")),
        kPrMomByRangeMaxLength(pset.get<double>("PrMomByRangeMaxLength")),
        kMinPionNTrk(pset.get<size_t>("MinPionNTrk")),
        kMinPionTrkLength(pset.get<double>("MinPionTrkLength")),
        kMaxPionTrkLength(pset.get<double>("MaxPionTrkLength")),
        kMinPionTrkLengthNotContained(pset.get<double>("MinPionTrkLengthNotContained")),
        kDistanceToWallThreshold(pset.get<double>("DistanceToWallThreshold")),
        kRecombFactor(pset.get<double>("RecombFactor")),
        kPlaneToUse(pset.get<unsigned int>("Plane")),
        fPFParticleLabel(pfpLabel),
        fTrackLabel(trackLabel),
        fShowerLabel(showerLabel),
        fHitLabel(hitLabel),
        fTrackToHitLabel(trackToHitLabel),
        fShowerToHitLabel(showerToHitLabel),
        fHitToSpacePointLabel(hitToSpacePointLabel),
        fParticleIDModuleLabel(pset.get<std::string>("ParticleIDModuleLabel"))
    {
        if (kPlaneToUse > 2)
            kPlaneToUse = static_cast<unsigned int>(geo::kUnknown);
        kPlane = static_cast<geo::View_t>(kPlaneToUse);


    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Track> ParticleSelectionAlg::GetLongestTrack(
            const art::Event &event,
            const std::vector<art::Ptr<recob::Track>> &tracks,
            const bool isTrackOnly = true)
    {
        art::Ptr<recob::Track> pTrack{};
        double longestLength(std::numeric_limits<double>::lowest());
        for (unsigned int iTrack = 0; iTrack < tracks.size(); ++iTrack)
        {
            if (isTrackOnly && !this->IsTrack(event, tracks[iTrack]))
                continue;
            const double length(tracks[iTrack]->Length());
            if (length-longestLength > std::numeric_limits<double>::epsilon())
            {
                longestLength = length;
                pTrack = tracks[iTrack];
            }
        }

        fLepTrack = pTrack;

        return pTrack;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Track> ParticleSelectionAlg::GetLongestTrack(const art::Event &event, const bool isTrackOnly)
    {
        art::Ptr<recob::Track> pTrack{};
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        if (0 == tracks.size())
            return pTrack;

        return ParticleSelectionAlg::GetLongestTrack(event, tracks, isTrackOnly);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Track> ParticleSelectionAlg::GetLongestTrackPID(const art::Event &event)
    {

        fTrkIdToPIDAMap = ParticleSelectionAlg::GenMapPIDAScore(event);
        fTotalCaloEvent = ParticleSelectionAlg::GetTotalCaloEvent(event);
        art::Ptr<recob::Track> pTrack{};

        std::vector<art::Ptr<recob::Track>> tracks(GenMuonCandidates(event));

        if (0 == tracks.size())
            return pTrack;

        const bool isTrackOnly = false;
        return ParticleSelectionAlg::GetLongestTrack(event, tracks, isTrackOnly);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<art::Ptr<recob::Track>> ParticleSelectionAlg::GenMuonCandidates(const art::Event &event)
    {

        std::vector<art::Ptr<recob::Track>> candidates(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        std::vector<art::Ptr<recob::Track>> prevCandidates = candidates;

        if (candidates.size() <= 1)
            return candidates;

        std::vector<std::function<void()>> filters = {
            // Remove tracks that are too big
            [&]() { ParticleSelectionAlg::applyMuLengthFilter(candidates, kMaxMuLength); },
            // Remove when PIDA is > threshold (if PIDA is not valid, do nothing)
            [&]() { ParticleSelectionAlg::applyMuMaxPIDA(candidates, kMaxMuPIDA); },
            // Removing PFPs that are assigned as shower in events with where Total Calorimetric Enegy is between user defined value
            [&]() { ParticleSelectionAlg::applyMuShowerCutCalo(event, candidates, kMinMuTotalCalo, kMaxMuTotalCalo); },
            // Removing PFPs that are shower and have PIDA score < threshold
            [&]() { ParticleSelectionAlg::applyMuShowerCutPIDA(event, candidates, kMinMuPIDAShower); },
            // If not all tracks are contained, remove the tracks that are contained
            [&]() { ParticleSelectionAlg::applyMuCutContained(event, candidates); },
            // Remove PFPs that are showers
            [&]() { ParticleSelectionAlg::applyMuRemoveShowers(event, candidates); },
            // Remove AGAIN when PIDA is > lower_threshold (if PIDA is not valid, do nothing)
            [&]() { ParticleSelectionAlg::applyMuMaxPIDA(candidates, kMaxMuPIDAAggressive); },
        };

        // After each filter, check if the number of candidates is 0 or 1
        // If there are no candidates after filtering, return the previous candidates
        // If there is only one candidate, return it
        for (auto &filter: filters)
        {
            prevCandidates = candidates;
            filter();
            if (candidates.empty()) return prevCandidates;
            if (candidates.size() == 1) return candidates;
        }

        return candidates;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Shower> ParticleSelectionAlg::GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                                         detinfo::DetectorPropertiesData const& detProp,
                                                                         const art::Event &event)
    {
        const std::vector<art::Ptr<recob::Shower>> showers(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));
        return ParticleSelectionAlg::GetHighestChargeShower(clockData, detProp, event, showers);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Shower> ParticleSelectionAlg::GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                                         detinfo::DetectorPropertiesData const& detProp,
                                                                         const art::Event &event,
                                                                         const std::vector<art::Ptr<recob::Shower>> showers,
                                                                         const geo::View_t tPlane)
    {
        art::Ptr<recob::Shower> pShower{};
        if (0 == showers.size())
            return pShower;

        // In case the plane is set to geo::Unknown, the plane with most hits will be used for each shower
        fShwIdToBestPlaneMap = ParticleSelectionAlg::GenMapBestPlaneShower(event, tPlane);

        double maxCharge(std::numeric_limits<double>::lowest());
        for (unsigned int iShower = 0; iShower < showers.size(); ++iShower)
        {
            const std::vector<art::Ptr<recob::Hit>>
                showerHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(showers[iShower],
                                event, fShowerToHitLabel), fShwIdToBestPlaneMap[showers[iShower].key()]));
            const double showerCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, showerHits));
            if (showerCharge-maxCharge > std::numeric_limits<double>::epsilon())
            {
                maxCharge = showerCharge;
                pShower = showers[iShower];
            }
        }

        // Store the shower and the associated track (if available)
        fLepShower = pShower;
        fLepTrack = ParticleSelectionAlg::GetTrackFromShower(event, pShower);

        return pShower;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Shower> ParticleSelectionAlg::GetHighestChargeShowerPID(const art::Event &event)
    {
        std::vector<art::Ptr<recob::Shower>> showers(GenElectronCandidates(event));
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
        return ParticleSelectionAlg::GetHighestChargeShower(clockData, detProp, event, showers, kPlane);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<art::Ptr<recob::Shower>> ParticleSelectionAlg::GenElectronCandidates(const art::Event &event)
    {

        fTotalCaloEvent = ParticleSelectionAlg::GetTotalCaloEvent(event);
        fShwIdToPIDAMap = ParticleSelectionAlg::GenMapPIDAScoreShw(event);

        std::vector<art::Ptr<recob::Shower>> candidates(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));

        // Removing showers if they are considered as Tracks AND the total calorimetric energy of the event is above a threshold
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                    [this, &event](const art::Ptr<recob::Shower> &s) {
                    return (this->IsTrack(event, s) && fTotalCaloEvent >= this->kMaxETotalCaloForTracks);
                    }), candidates.end());

        // removing showers if they are considered as Tracks AND the PIDA score is above a threshold
        candidates.erase(std::remove_if(candidates.begin(), candidates.end(),
                    [this, &event](const art::Ptr<recob::Shower> &s) {
                    return (this->IsTrack(event, s) && this->GetPIDAScore(s) > this->kMaxEPIDA);
                    }), candidates.end());

        return candidates;

    }
    //------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<art::Ptr<recob::Track>> ParticleSelectionAlg::GenProtonCandidates(const art::Event &event)
    {
        std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));

        if (0 == tracks.size())
            return tracks;

        // If lepton was not searched, these maps are empty, so we need to generate them
        if (fTrkIdToPIDAMap.empty())
            fTrkIdToPIDAMap = ParticleSelectionAlg::GenMapPIDAScore(event);
        if (fTotalCaloEvent == std::numeric_limits<double>::lowest())
            fTotalCaloEvent = ParticleSelectionAlg::GetTotalCaloEvent(event);
        if (fTrkIdToCaloMap.empty())
            fTrkIdToCaloMap = ParticleSelectionAlg::GenMapTracksCalo(event);
        if (fTrkIdToMomMap.empty())
            fTrkIdToMomMap = ParticleSelectionAlg::GenMapTracksMom(event);


        // If there is a lepton track, remove it from the list of tracks
        if (this->GetLepTrack().isAvailable())
        {
            tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                        [this](const art::Ptr<recob::Track> &t) {
                        return t.key() == this->GetLepTrack().key();
                        }), tracks.end());
        }

        // Removing tracks that are not considered as tracks by LArPandora when their PIDA score is below a threshold
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t) {
                    return (this->IsTrack(event, t) && this->GetPIDAScore(t) <= this->kMinPrPIDATrack);
                    }), tracks.end());

        // Removing tracks that are considered as showers by LArPandora when their PIDA score is below a threshold
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrack(event, t) && this->GetPIDAScore(t) <= this->kMinPrPIDAShower);
                    }), tracks.end());

        // Removing tracks that have calorimetric energy above a threshold
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return (this->GetTrackCalo(t) >= this->kMaxPrTrkCalo);
                    }), tracks.end());

        // Removing tracks that have momentum above a threshold
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return (this->GetTrackMom(t) >= this->kMaxPrTrkMom);
                    }), tracks.end());

        fPrTrack = tracks;
        fPrDone = true;
        return tracks;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<art::Ptr<recob::Track>> ParticleSelectionAlg::GenPionCandidates(const art::Event &event)
    {
        std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));

        // Removing lepton candidate if any
        if (this->GetLepTrack().isAvailable())
        {
            tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                        [&](const art::Ptr<recob::Track> &t) {
                        return t.key() == this->GetLepTrack().key();
                        }), tracks.end());
        }

        // If there are not enough tracks, no pion candidates
        if (tracks.size() <= kMinPionNTrk)
        {
            std::vector<art::Ptr<recob::Track>> empty;
            fPiTrack = empty;
            return empty;
        }

        // Enforce the presence of proton candidates.
        // This is done because pion selection is rather weak
        std::vector<art::Ptr<recob::Track>> trkprotons;
        if (fPrDone)
            trkprotons = ParticleSelectionAlg::GenProtonCandidates(event);
        else
            trkprotons = ParticleSelectionAlg::GetPrTracks();


        // Removing proton candidates
        if (!trkprotons.empty())
        {
            std::unordered_set<art::Ptr<recob::Track>::key_type> keys_to_remove;
            for (const art::Ptr<recob::Track> &t : trkprotons) keys_to_remove.insert(t.key());

            tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                        [&keys_to_remove](const art::Ptr<recob::Track> &t) {
                        return keys_to_remove.count(t.key());
                        }), tracks.end());
        }

        // Removing tracks that are not considered as tracks by LArPandora
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrack(event, t));
                    }), tracks.end());

        // Removing tracks that are not in the length range
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return !(t->Length() >= this->kMinPionTrkLength && t->Length() <= this->kMaxPionTrkLength);
                    }), tracks.end());

        // Removing tracks that are not contained and are not long enough
        std::optional<bool> allContained(GenContainmentInfo(event, tracks));
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrkContained(t) && t->Length() < this->kMinPionTrkLengthNotContained);
                    }), tracks.end());

        fPiTrack = tracks;

        return tracks;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Track>::key_type, double> ParticleSelectionAlg::GenMapPIDAScore(const art::Event &event)
    {

        std::map<art::Ptr<recob::Track>::key_type, double> tTrkIdToPIDAMap;

        art::Handle< std::vector<recob::Track>> trackListHandle = event.getHandle< std::vector<recob::Track>>(fTrackLabel);
        std::unique_ptr<art::FindManyP<anab::ParticleID>> fmpid;
        fmpid = std::make_unique<art::FindManyP<anab::ParticleID>>(trackListHandle, event, fParticleIDModuleLabel);
        if(!fmpid->isValid())
        {
            mf::LogWarning("ParticleSelectionAlg") << " Could not find any ParticleID association" << std::endl;
            return tTrkIdToPIDAMap;
        }

        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));

        for (const art::Ptr<recob::Track> &track : tracks)
        {

            const std::vector< art::Ptr<anab::ParticleID>> pids = fmpid->at(track.key());
            if (pids.empty())
            {
                tTrkIdToPIDAMap[track.key()] = 0;
                continue;
            }

            int biggestNdf(std::numeric_limits<int>::lowest());
            geo::View_t tplane = kPlane;
            // if plane is set to unknown, search the plane with most points
            // when PIDA was evaluated
            if (tplane == geo::kUnknown)
            {
                for (size_t ipid = 0; ipid < pids.size(); ++ipid)
                {
                    if (!pids[ipid]->PlaneID().isValid) continue;
                    geo::PlaneID::PlaneID_t planenum = pids[ipid]->PlaneID().Plane;
                    if (planenum<0||planenum>2) continue;
                    const std::vector< anab::sParticleIDAlgScores> pidScore = pids[ipid]->ParticleIDAlgScores();
                    for (const anab::sParticleIDAlgScores &pScore : pidScore)
                    {
                        if (pScore.fAssumedPdg == 0 && pScore.fNdf > biggestNdf)
                        {
                            biggestNdf = pScore.fNdf;
                            tplane = static_cast<geo::View_t>(ipid);
                        }
                    }
                }
                // If failed to find a plane, set to PIDA score to zero
                if (tplane == geo::kUnknown)
                {
                    tTrkIdToPIDAMap[track.key()] = 0;
                    continue;
                }
            }

            double pidaScore = 0;
            if(pids[tplane]->PlaneID().isValid)
            {
                const std::vector< anab::sParticleIDAlgScores> pidScore = pids[tplane]->ParticleIDAlgScores();
                for (const anab::sParticleIDAlgScores &pScore : pidScore)
                {
                    if (pScore.fAssumedPdg == 0)
                    {
                        pidaScore = pScore.fValue;
                    }
                }
            }
            tTrkIdToPIDAMap[track.key()] = pidaScore;
        }
        return tTrkIdToPIDAMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Shower>::key_type, double> ParticleSelectionAlg::GenMapPIDAScoreShw(const art::Event &event)
    {
        std::map<art::Ptr<recob::Shower>::key_type, double> tShwIdToPIDAMap;
        const std::vector<art::Ptr<recob::Shower>> showers(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));
        if (fTrkIdToPIDAMap.size() == 0)
        {
            fTrkIdToPIDAMap = ParticleSelectionAlg::GenMapPIDAScore(event);
        }
        for (const art::Ptr<recob::Shower> &shower : showers)
        {
            art::Ptr<recob::Track> shwTrk = ParticleSelectionAlg::GetTrackFromShower(event, shower);
            if (!shwTrk.isAvailable())
                tShwIdToPIDAMap[shower.key()] = 0;
            else
                tShwIdToPIDAMap[shower.key()] = ParticleSelectionAlg::GetPIDAScore(shwTrk);
        }

        return tShwIdToPIDAMap;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::GenContainmentInfo(const art::Event &event)
    {
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        return ParticleSelectionAlg::GenContainmentInfo(event, tracks);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::GenContainmentInfo(const art::Event &event, const std::vector<art::Ptr<recob::Track>> &tracks)
    {
        bool allContained = true;
        hitListHandle    = event.getHandle< std::vector<recob::Hit> >(fHitLabel);
        fmsp  = std::make_unique<art::FindManyP<recob::SpacePoint>>(hitListHandle, event, fHitToSpacePointLabel);
        for (const art::Ptr<recob::Track> &track : tracks)
        {
            if (fTrkIdToContainmentMap.count(track.key()) == 0)
            {
                std::vector<art::Ptr<recob::Hit>> trkhits = (dune_ana::DUNEAnaTrackUtils::GetHits(track, event, fTrackLabel));
                bool trkIsContained = ParticleSelectionAlg::IsContained(trkhits, event);
                fTrkIdToContainmentMap[track.key()] = trkIsContained;
            }
            allContained &= fTrkIdToContainmentMap[track.key()];
        }
        return allContained;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetPIDAScore(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToPIDAMap.count(track.key()))
        {
            return fTrkIdToPIDAMap[track.key()];
        }
        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetPIDAScore(const art::Ptr<recob::Shower> &shower)
    {
        if (fShwIdToPIDAMap.count(shower.key()))
        {
            return fShwIdToPIDAMap[shower.key()];
        }
        return 0;
    }


    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetTrackCalo(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToCaloMap.count(track.key()))
        {
            return fTrkIdToCaloMap[track.key()];
        }
        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetTrackMom(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToMomMap.count(track.key()))
        {
            return fTrkIdToMomMap[track.key()];
        }
        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    double ParticleSelectionAlg::GetTotalCaloEvent(const art::Event &event)
    {
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

        geo::View_t tplane = kPlane;
        if (tplane == geo::kUnknown)
        {
            size_t maxNhits(std::numeric_limits<size_t>::lowest());
            for (size_t ipl = 0; ipl < kNplanes; ++ipl)
            {
                const std::vector<art::Ptr<recob::Hit>>
                    eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event,
                                    fHitLabel),ipl));
                const size_t nhits = eventHits.size();
                if (nhits-maxNhits > std::numeric_limits<size_t>::epsilon())
                {
                    maxNhits = nhits;
                    tplane = static_cast<geo::View_t>(ipl);
                }
            }
        }

        const std::vector<art::Ptr<recob::Hit>>
            eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event,
                            fHitLabel), tplane));
        const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits));

        const double CaloEnergy(this->CalculateEnergyFromCharge(eventObservedCharge, tplane));
        return CaloEnergy;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Track> ParticleSelectionAlg::GetTrackFromShower(const art::Event &event, const art::Ptr<recob::Shower> &shower)
    {
        const art::Ptr<recob::PFParticle> pfp(dune_ana::DUNEAnaShowerUtils::GetPFParticle(shower, event, fShowerLabel));
        if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, event, fPFParticleLabel, fTrackLabel))
        {
            return dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, event, fPFParticleLabel, fTrackLabel);

        }
        return art::Ptr<recob::Track>{};
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Track>::key_type, double> ParticleSelectionAlg::GenMapTracksCalo(const art::Event &event)
    {
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        std::map<art::Ptr<recob::Track>::key_type, double> tTrkIdToCaloMap;

        for (const art::Ptr<recob::Track> &track : tracks)
        {

            geo::View_t tplane = kPlane;

            if (tplane == geo::kUnknown)
            {
                size_t maxNhits(std::numeric_limits<size_t>::lowest());
                for (size_t ipl = 0; ipl < kNplanes; ++ipl)
                {
                    const std::vector<art::Ptr<recob::Hit>>
                        trkHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(track,
                                        event, fTrackLabel), ipl));
                    const size_t nhits = trkHits.size();
                    if (nhits-maxNhits > std::numeric_limits<size_t>::epsilon())
                    {
                        maxNhits = nhits;
                        tplane = static_cast<geo::View_t>(ipl);
                    }
                }
            }

            const std::vector<art::Ptr<recob::Hit>>
                trkHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(track,
                                event, fTrackLabel), tplane));
            const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, trkHits));

            const double CaloEnergy(this->CalculateEnergyFromCharge(eventObservedCharge, tplane));

            tTrkIdToCaloMap[track.key()] = CaloEnergy;

        }

        return tTrkIdToCaloMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Track>::key_type, double> ParticleSelectionAlg::GenMapTracksMom(const art::Event &event)
    {
        std::map<art::Ptr<recob::Track>::key_type, double> tTrkIdToMomMap;
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        trkf::TrackMomentumCalculator trkmrange{kPrMomByRangeMinLength, kPrMomByRangeMaxLength};
        for (const art::Ptr<recob::Track> &track : tracks)
        {
            double trkmomrangepr = trkmrange.GetTrackMomentum(track->Length(),2212);
            tTrkIdToMomMap[track.key()] = trkmomrangepr;
        }
        return tTrkIdToMomMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Shower>::key_type, geo::View_t> ParticleSelectionAlg::GenMapBestPlaneShower(const art::Event &event, const geo::View_t tPlane)
    {
        std::map<art::Ptr<recob::Shower>::key_type, geo::View_t> tShwIdToBestPlaneMap;
        const std::vector<art::Ptr<recob::Shower>> showers(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));
        for (const art::Ptr<recob::Shower> &shower : showers)
        {
            geo::View_t tplane = tPlane;
            if (tplane == geo::kUnknown)
            {
                size_t maxNhits(std::numeric_limits<size_t>::lowest());
                for (size_t ipl = 0; ipl < kNplanes; ++ipl)
                {
                    const std::vector<art::Ptr<recob::Hit>>
                        shwHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(shower,
                                        event, fShowerLabel), ipl));
                    const size_t nhits = shwHits.size();
                    if (nhits-maxNhits > std::numeric_limits<size_t>::epsilon())
                    {
                        maxNhits = nhits;
                        tplane = static_cast<geo::View_t>(ipl);
                    }
                }
            }
            tShwIdToBestPlaneMap[shower.key()] = tplane;
        }

        return tShwIdToBestPlaneMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    double ParticleSelectionAlg::CalculateEnergyFromCharge(const double charge, const unsigned short plane )
    {
        return fCalorimetryAlg.ElectronsFromADCArea(charge, plane)*1./kRecombFactor/util::kGeVToElectrons;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::IsTrack(const art::Event &event, const art::Ptr<recob::Track> &track)
    {
        const art::Ptr< recob::PFParticle > pfp(dune_ana::DUNEAnaTrackUtils::GetPFParticle(track, event, fTrackLabel));
        return lar_pandora::LArPandoraHelper::IsTrack(pfp);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::IsTrack(const art::Event &event, const art::Ptr<recob::Shower> shower)
    {
        const art::Ptr< recob::PFParticle > pfp(dune_ana::DUNEAnaShowerUtils::GetPFParticle(shower, event, fShowerLabel));
        return lar_pandora::LArPandoraHelper::IsTrack(pfp);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::IsPointContained(const double x, const double y, const double z)
    {
        geo::Point_t const position{x,y,z};
        art::ServiceHandle<geo::Geometry> fGeometry;

        geo::TPCID tpcID(fGeometry->FindTPCAtPosition(position));
        if (!(fGeometry->HasTPC(tpcID)))
            return false;

        double minX(std::numeric_limits<double>::max());
        double maxX(std::numeric_limits<double>::lowest());
        double minY(std::numeric_limits<double>::max());
        double maxY(std::numeric_limits<double>::lowest());
        double minZ(std::numeric_limits<double>::max());
        double maxZ(std::numeric_limits<double>::lowest());
        for (auto const& tpc : fGeometry->Iterate<geo::TPCGeo>())
        {
            minX = std::min(minX,tpc.MinX());
            maxX = std::max(maxX,tpc.MaxX());
            minY = std::min(minY,tpc.MinY());
            maxY = std::max(maxY,tpc.MaxY());
            minZ = std::min(minZ,tpc.MinZ());
            maxZ = std::max(maxZ,tpc.MaxZ());
        }

        minX += kDistanceToWallThreshold;
        maxX -= kDistanceToWallThreshold;
        minY += kDistanceToWallThreshold;
        maxY -= kDistanceToWallThreshold;
        minZ += kDistanceToWallThreshold;
        maxZ -= kDistanceToWallThreshold;

        if (x - minX < -1.*std::numeric_limits<double>::epsilon() ||
                x - maxX > std::numeric_limits<double>::epsilon())
            return false;
        if (y - minY < -1.*std::numeric_limits<double>::epsilon() ||
                y - maxY > std::numeric_limits<double>::epsilon())
            return false;
        if (z - minZ < -1.*std::numeric_limits<double>::epsilon() ||
                z - maxZ > std::numeric_limits<double>::epsilon())
            return false;

        return true;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::IsContained(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Event &event)
    {

        for (unsigned int iHit = 0; iHit < hits.size(); ++iHit)
        {
            std::vector<art::Ptr<recob::SpacePoint> > spacePoints;
            if(fmsp->isValid())
            {
                spacePoints = fmsp->at(hits[iHit].key());
            }
            for (unsigned int iSpacePoint = 0; iSpacePoint < spacePoints.size(); ++iSpacePoint)
            {
                const art::Ptr<recob::SpacePoint> spacePoint(spacePoints[iSpacePoint]);
                if (!(ParticleSelectionAlg::IsPointContained(spacePoint->XYZ()[0],spacePoint->XYZ()[1],spacePoint->XYZ()[2])))
                    return false;
            }
        }
        return true;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::IsTrkContained(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToContainmentMap.count(track.key()))
        {
            return fTrkIdToContainmentMap[track.key()];
        }
        else
        {
            throw art::Exception(art::errors::LogicError) << "The map between \
                track id and containment information does not have this track \
                id, this shouldn't happen.\n";

        }
        return false;
    }


    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuLengthFilter(
            std::vector<art::Ptr<recob::Track>>& tracks,
            const double maxlen)
    {
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [&maxlen](const art::Ptr<recob::Track> &t) {
                    return t->Length() > maxlen;
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuMaxPIDA(
            std::vector<art::Ptr<recob::Track>>& tracks,
            const double maxPIDA)
    {
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &maxPIDA](const art::Ptr<recob::Track> &t) {
                    return this->GetPIDAScore(t) > maxPIDA;
                    }), tracks.end());

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuShowerCutCalo(
            const art::Event &event,
            std::vector<art::Ptr<recob::Track>>& tracks,
            const double minCalo,
            const double maxCalo)
    {
        // Nothing to be filtered in this event
        if ( !(this->fTotalCaloEvent >= minCalo && this->fTotalCaloEvent <= maxCalo) )
            return;
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrack(event, t));
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuShowerCutPIDA(
            const art::Event &event,
            std::vector<art::Ptr<recob::Track>>& tracks,
            const double minPIDA)
    {
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event, &minPIDA](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrack(event, t) && this->GetPIDAScore(t) < minPIDA);
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuCutContained(
            const art::Event &event,
            std::vector<art::Ptr<recob::Track>>& tracks
            )
    {
        if (tracks.size()<=1) return; // Shouldn't get here...

        const bool allContained{GenContainmentInfo(event, tracks)};
        if (allContained) return; // no filter needed

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return (this->IsTrkContained(t) && fTotalCaloEvent >= kMaxMuContainedCalo);
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuRemoveShowers(
            const art::Event &event,
            std::vector<art::Ptr<recob::Track>>& tracks
            )
    {
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrack(event, t));
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

}

