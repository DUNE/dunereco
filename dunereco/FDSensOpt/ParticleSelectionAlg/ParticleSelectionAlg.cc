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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

//DUNE
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"

#include "dunereco/FDSensOpt/ParticleSelectionAlg/ParticleSelectionAlg.h"

#pragma GCC diagnostic ignored "-Wunused-variable"
namespace dune
{
    ParticleSelectionAlg::ParticleSelectionAlg(
            fhicl::ParameterSet const& pset,
            const std::string &trackLabel,
            const std::string &showerLabel,
            const std::string &hitLabel,
            const std::string &trackToHitLabel,
            const std::string &showerToHitLabel,
            const std::string &hitToSpacePointLabel) :
        fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
        fTrackLabel(trackLabel),
        fShowerLabel(showerLabel),
        fHitLabel(hitLabel),
        fTrackToHitLabel(trackToHitLabel),
        fShowerToHitLabel(showerToHitLabel),
        fHitToSpacePointLabel(hitToSpacePointLabel)
        // fDistanceToWallThreshold(pset.get<double>("DistanceToWallThreshold"))
    {
        fPlane = geo::kUnknown;
        fRecombFactor = 0.63;
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
        fAllTracksContained = ParticleSelectionAlg::GenContainmentInfo(event);
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
            [&]() { ParticleSelectionAlg::applyMuLengthFilter(candidates, fMaxMuLenght); },
            // Remove when PIDA is > threshold (if PIDA is not valid, do nothing)
            [&]() { ParticleSelectionAlg::applyMuMaxPIDA(candidates, 13); },
            // Removing PFPs that are assigned as shower in events with where Total Calorimetric Enegy is between user defined value
            [&]() { ParticleSelectionAlg::applyMuShowerCutCalo(event, candidates, 0.3, 2); },
            // Removing PFPs that are shower and have PIDA score < threshold
            [&]() { ParticleSelectionAlg::applyMuShowerCutPIDA(event, candidates, 5); },
            // If not all tracks are contained, remove the tracks that are contained
            [&]() { ParticleSelectionAlg::applyMuCutContained(candidates); },
            // Remove PFPs that are showers
            [&]() { ParticleSelectionAlg::applyMuRemoveShowers(event, candidates); },
            // Remove AGAIN when PIDA is > lower_threshold (if PIDA is not valid, do nothing)
            [&]() { ParticleSelectionAlg::applyMuMaxPIDA(candidates, 10); },
        };

        for (auto &filter: filters){
            prevCandidates = candidates;
            filter();
            if (candidates.empty()) return prevCandidates;
            if (candidates.size() == 1) return candidates;
        }

        return candidates;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::vector<art::Ptr<recob::Track>> ParticleSelectionAlg::GenProtonCandidates(const art::Event &event)
    {
        if (fTrkIdToPIDAMap.size() == 0)
            fTrkIdToPIDAMap = ParticleSelectionAlg::GenMapPIDAScore(event);
        if (fTotalCaloEvent == std::numeric_limits<double>::lowest())
            fTotalCaloEvent = ParticleSelectionAlg::GetTotalCaloEvent(event);
        if (fTrkIdToContainmentMap.size() == 0)
            fAllTracksContained = ParticleSelectionAlg::GenContainmentInfo(event);
        if (fTrkIdToCaloMap.size() == 0)
            fTrkIdToCaloMap = ParticleSelectionAlg::GenMapTracksCalo(event);
        if (fTrkIdToMomMap.size() == 0)
            fTrkIdToMomMap = ParticleSelectionAlg::GenMapTracksMom(event);

        std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));

        if (0 == tracks.size())
            return tracks;
        
        if (this->GetLepTrack().isAvailable()){
            tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                        [this](const art::Ptr<recob::Track> &t){
                        return t.key() == this->GetLepTrack().key(); 
                        }), tracks.end());
        }

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t){
                    return (this->IsTrack(event, t) && this->GetPIDAScore(t) <= 10);
                    }), tracks.end());

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t){
                    return (!this->IsTrack(event, t) && this->GetPIDAScore(t) <= 13);
                    }), tracks.end());

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t){
                    return (this->GetTrackCalo(t) >= 1.5);
                    }), tracks.end());

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t){
                    return (this->GetTrackMom(t) >= 2.0);
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
        if (this->GetLepTrack().isAvailable()){
            tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                        [&](const art::Ptr<recob::Track> &t) {
                        return t.key() == this->GetLepTrack().key(); 
                        }), tracks.end());
        }

        if (tracks.size() <= 2){
            std::vector<art::Ptr<recob::Track>> empty;
            fPiTrack = empty;
            return empty;
        }

        std::vector<art::Ptr<recob::Track>> trkprotons;
        if (fPrDone)
            trkprotons = ParticleSelectionAlg::GenProtonCandidates(event);
        else
            trkprotons = ParticleSelectionAlg::GetPrTracks();


        if (trkprotons.size() > 0){
            std::unordered_set<art::Ptr<recob::Track>::key_type> keys_to_remove;
            for (const art::Ptr<recob::Track> &t : trkprotons) keys_to_remove.insert(t.key());

            tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                        [&keys_to_remove](const art::Ptr<recob::Track> &t) {
                        return keys_to_remove.count(t.key());  // O(1) lookup
                        }), tracks.end());
        }

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrack(event, t));
                    }), tracks.end());

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return !(t->Length() >= 3.5 && t->Length() <= 350);
                    }), tracks.end());

        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t) {
                    return (!this->IsTrkContained(t) && t->Length() < 10);
                    }), tracks.end());

        fPiTrack = tracks;

        return tracks;

    }


    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Shower> ParticleSelectionAlg::GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                                         detinfo::DetectorPropertiesData const& detProp,
                                                                         const art::Event &event)
    {
        art::Ptr<recob::Shower> pShower{};
        const std::vector<art::Ptr<recob::Shower>> showers(dune_ana::DUNEAnaEventUtils::GetShowers(event, fShowerLabel));
        return ParticleSelectionAlg::GetHighestChargeShower(clockData, detProp, event, showers);

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    art::Ptr<recob::Shower> ParticleSelectionAlg::GetHighestChargeShower(detinfo::DetectorClocksData const& clockData,
                                                                         detinfo::DetectorPropertiesData const& detProp,
                                                                         const art::Event &event,
                                                                         const std::vector<art::Ptr<recob::Shower>> showers
                                                                         )
    {
        art::Ptr<recob::Shower> pShower{};
        if (0 == showers.size())
            return pShower;

        double maxCharge(std::numeric_limits<double>::lowest());
        for (unsigned int iShower = 0; iShower < showers.size(); ++iShower)
        {
            const std::vector<art::Ptr<recob::Hit>>
                showerHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(showers[iShower],
                                event,fShowerToHitLabel),2));
            const double showerCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, showerHits));
            if (showerCharge-maxCharge > std::numeric_limits<double>::epsilon())
            {
                maxCharge = showerCharge;
                pShower = showers[iShower];
            }
        }
        return pShower;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Track>::key_type, double> ParticleSelectionAlg::GenMapPIDAScore(const art::Event &event)
    {

        std::map<art::Ptr<recob::Track>::key_type, double> tTrkIdToPIDAMap;

        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        trackListHandle = event.getHandle< std::vector<recob::Track>>(fTrackLabel);
        fmpid = std::make_unique<art::FindManyP<anab::ParticleID>>(trackListHandle, event, fParticleIDModuleLabel);
        if(!fmpid->isValid()) {
            mf::LogWarning("ParticleSelectionAlg") << " Could not find any ParticleID association" << std::endl;
            return tTrkIdToPIDAMap;
        }

        for (const art::Ptr<recob::Track> &track : tracks){

            const std::vector< art::Ptr<anab::ParticleID>> pids = this->fmpid->at(track.key());

            int biggestNdf(std::numeric_limits<int>::lowest());
            geo::View_t kplane = fPlane;
            // if plane is set to unknown, search the plane with most points
            // when PIDA was evaluated
            if (kplane == geo::kUnknown){
                for (size_t ipid = 0; ipid < pids.size(); ++ipid){
                    if (!pids[ipid]->PlaneID().isValid) continue;
                    geo::PlaneID::PlaneID_t planenum = pids[ipid]->PlaneID().Plane;
                    if (planenum<0||planenum>2) continue;
                    auto pidScore = pids[ipid]->ParticleIDAlgScores();
                    anab::sParticleIDAlgScores pScore = pidScore.back(); 
                    if (pScore.fAssumedPdg != 0){
                        throw art::Exception(art::errors::LogicError) << "PIDA is not the last PID algorith! \n";
                    }
                    const int pidndf = pScore.fNdf;
                    if (pidndf-biggestNdf > std::numeric_limits<int>::epsilon())
                    {
                        biggestNdf = pidndf;
                        kplane = static_cast<geo::View_t>(ipid);
                    }
                }
            }
            auto pidScore = pids[kplane]->ParticleIDAlgScores();
            anab::sParticleIDAlgScores pScore = pidScore.back(); 
            if (pScore.fAssumedPdg != 0){
                throw art::Exception(art::errors::LogicError) << "PIDA is not the last PID algorith! \n";
            }
            tTrkIdToPIDAMap[track.key()] = pScore.fValue;
        }
        return tTrkIdToPIDAMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool ParticleSelectionAlg::GenContainmentInfo(const art::Event &event)
    {
        bool allContained = true;
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        for (const art::Ptr<recob::Track> &track : tracks){

            std::vector<art::Ptr<recob::Hit>> trkhits = (dune_ana::DUNEAnaTrackUtils::GetHits(track, event, fTrackLabel));
            bool trkIsContained = ParticleSelectionAlg::IsContained(trkhits, event);
            allContained &= trkIsContained;
            fTrkIdToContainmentMap[track.key()] = trkIsContained; 
        }

        return allContained;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetPIDAScore(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToPIDAMap.count(track.key())){
            return fTrkIdToPIDAMap[track.key()];
        }
        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetTrackCalo(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToCaloMap.count(track.key())){
            return fTrkIdToCaloMap[track.key()];
        }
        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const double ParticleSelectionAlg::GetTrackMom(const art::Ptr<recob::Track> &track)
    {
        if (fTrkIdToMomMap.count(track.key())){
            return fTrkIdToMomMap[track.key()];
        }
        return 0;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    double ParticleSelectionAlg::GetTotalCaloEvent(const art::Event &event)
    {
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

        geo::View_t kplane = fPlane;
        if (kplane == geo::kUnknown){
            size_t maxNhits(std::numeric_limits<size_t>::lowest());
            for (size_t ipl = 0; ipl < kNplanes; ++ipl){
                const std::vector<art::Ptr<recob::Hit>>
                    eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event,
                                    fHitLabel),ipl));
                const size_t nhits = eventHits.size();
                if (nhits-maxNhits > std::numeric_limits<size_t>::epsilon())
                {
                    maxNhits = nhits;
                    kplane = static_cast<geo::View_t>(ipl);
                }
            }
        }

        const std::vector<art::Ptr<recob::Hit>>
            eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event,
                            fHitLabel), kplane));
        const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits));

        const double CaloEnergy(this->CalculateEnergyFromCharge(eventObservedCharge, kplane));
        return CaloEnergy;

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Track>::key_type, double> ParticleSelectionAlg::GenMapTracksCalo(const art::Event &event)
    {
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        std::map<art::Ptr<recob::Track>::key_type, double> tTrkIdToCaloMap; 
        
        for (const art::Ptr<recob::Track> &track : tracks){

            geo::View_t kplane = fPlane;

            if (kplane == geo::kUnknown){
                size_t maxNhits(std::numeric_limits<size_t>::lowest());
                for (size_t ipl = 0; ipl < kNplanes; ++ipl){
                    const std::vector<art::Ptr<recob::Hit>>
                        trkHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(track,
                                        event, fTrackLabel), ipl));
                    const size_t nhits = trkHits.size();
                    if (nhits-maxNhits > std::numeric_limits<size_t>::epsilon())
                    {
                        maxNhits = nhits;
                        kplane = static_cast<geo::View_t>(ipl);
                    }
                }
            }

            const std::vector<art::Ptr<recob::Hit>>
                trkHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(track,
                                event, fTrackLabel), kplane));
            const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, trkHits));

            const double CaloEnergy(this->CalculateEnergyFromCharge(eventObservedCharge, kplane));

            tTrkIdToCaloMap[track.key()] = CaloEnergy;

        }
    
        return tTrkIdToCaloMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::map<art::Ptr<recob::Track>::key_type, double> ParticleSelectionAlg::GenMapTracksMom(const art::Event &event)
    {
        std::map<art::Ptr<recob::Track>::key_type, double> tTrkIdToMomMap;
        const std::vector<art::Ptr<recob::Track>> tracks(dune_ana::DUNEAnaEventUtils::GetTracks(event, fTrackLabel));
        trkf::TrackMomentumCalculator trkmrange{0,3000};
        for (const art::Ptr<recob::Track> &track : tracks){
            double trkmomrangepr = trkmrange.GetTrackMomentum(track->Length(),2212);
            tTrkIdToMomMap[track.key()] = trkmomrangepr;
        }
        return tTrkIdToMomMap;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    double ParticleSelectionAlg::CalculateEnergyFromCharge(const double charge, const unsigned short plane )
    {
        return fCalorimetryAlg.ElectronsFromADCArea(charge, plane)*1./fRecombFactor/util::kGeVToElectrons;
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

        minX += fDistanceToWallThreshold;
        maxX -= fDistanceToWallThreshold;
        minY += fDistanceToWallThreshold;
        maxY -= fDistanceToWallThreshold;
        minZ += fDistanceToWallThreshold;
        maxZ -= fDistanceToWallThreshold;

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
            std::vector<art::Ptr<recob::SpacePoint>>
                spacePoints(dune_ana::DUNEAnaHitUtils::GetSpacePoints(hits[iHit],
                            event, fHitLabel, fHitToSpacePointLabel));
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
        if (fTrkIdToContainmentMap.count(track.key())){
            return fTrkIdToContainmentMap[track.key()];
        }
        else{
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
                    [this, &maxPIDA](const art::Ptr<recob::Track> &t){
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
                    [this, &event](const art::Ptr<recob::Track> &t){
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
                    [this, &event, &minPIDA](const art::Ptr<recob::Track> &t){
                    return (!this->IsTrack(event, t) && this->GetPIDAScore(t) < minPIDA);
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuCutContained(
            std::vector<art::Ptr<recob::Track>>& tracks
            )
    {
        if (fAllTracksContained)
            return;
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this](const art::Ptr<recob::Track> &t){
                    return (this->IsTrkContained(t));
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const void ParticleSelectionAlg::applyMuRemoveShowers(
            const art::Event &event,
            std::vector<art::Ptr<recob::Track>>& tracks
            )
    {
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                    [this, &event](const art::Ptr<recob::Track> &t){
                    return (!this->IsTrack(event, t));
                    }), tracks.end());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

}

