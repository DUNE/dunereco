/**
*  @file   dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoPIDAlg.cc

*  @brief  Implementation file for the neutrino energy reconstruction algorithm with usage of PID
*  Written by Henrique Souza (hsouza@mib.infn.it)
*
*  $Log: $
*/

//ART
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
//DUNE
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaHitUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"

#include "dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"

namespace dune
{

    //-----------------------------------------------------------------------------------------------------------------------------------------
    dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyPID(const art::Ptr<recob::Track> &pMuonTrack,
            const art::Event &event,
            const int fLongestTrackMethod,
            const std::vector<art::Ptr<recob::Track>> &protonTracks,
            const std::vector<art::Ptr<recob::Track>> &pionTracks)
    {

        if (!pMuonTrack.isAvailable() || pMuonTrack.isNull())
        {
            mf::LogWarning("NeutrinoEnergyRecoPIDAlg") 
                << " Cannot access the muon track which is needed for this energy reconstruction method.\n"
                << "Swapping to energy reconstruction method " << kProtonPionHadronic << " for this calculation." << std::endl;
            return this->CalculateNeutrinoEnergyPID(event, protonTracks, pionTracks);
        }

        Point_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

        const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel),2));

        bool isContained(this->IsContained(muonHits, event));
        const double uncorrectedMuonMomentumMCS(this->CalculateUncorrectedMuonMomentumByMCS(pMuonTrack));
        const double muonMomentumMCS(this->CalculateLinearlyCorrectedValue(uncorrectedMuonMomentumMCS, fGradTrkMomMCS, fIntTrkMomMCS));
        const double muonMomentumRange(CalculateMuonMomentumByRange(pMuonTrack));
        const MuonContainmentStatus MuContainment = (isContained) ? MuonContainmentStatus::kIsContained : MuonContainmentStatus::kIsExiting;

        if (( !isContained && fLongestTrackMethod != LongestTrackMethod::kbyRange)
                || fLongestTrackMethod == LongestTrackMethod::kbyMCS) {
            if (uncorrectedMuonMomentumMCS > std::numeric_limits<double>::epsilon())
            {
                EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                        this->CalculateParticle4Momentum(kMuonMass, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                        kMuonAndHadronic, kMCS, MuContainment, fGradNuMuHadEnExit, fIntNuMuHadEnExit);
                return this->CalculateNeutrinoEnergyPID(muonHits, event, energyRecoInputHolder, protonTracks, pionTracks);
            }
            else
            {
                return this->CalculateNeutrinoEnergy(event);
            }
        }
        else
        {
            if ( fLongestTrackMethod != LongestTrackMethod::kbyRange
                    && uncorrectedMuonMomentumMCS > std::numeric_limits<double>::epsilon() 
                    && muonMomentumRange/muonMomentumMCS - fMuonRangeToMCSThreshold < -1.*std::numeric_limits<double>::epsilon() )
            {
                EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                        this->CalculateParticle4Momentum(kMuonMass, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                        kMuonAndHadronic, kMCS, MuContainment, fGradNuMuHadEnExit, fIntNuMuHadEnExit);

                return this->CalculateNeutrinoEnergyPID(muonHits, event, energyRecoInputHolder, protonTracks, pionTracks);
            }
            else
            {
                EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                        this->CalculateParticle4Momentum(kMuonMass, muonMomentumRange, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                        kMuonAndHadronic, kContained, MuContainment, fGradNuMuHadEnCont, fIntNuMuHadEnCont);

                return this->CalculateNeutrinoEnergyPID(muonHits, event, energyRecoInputHolder, protonTracks, pionTracks);
            }
        }

        throw art::Exception(art::errors::LogicError) << "Unable to determine how to calculate neutrino energy using muon track! \n";

    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyPID(const art::Ptr<recob::Shower> &pElectronShower, 
            const art::Event &event,
            const std::vector<art::Ptr<recob::Track>> &protonTracks,
            const std::vector<art::Ptr<recob::Track>> &pionTracks)
    {
        if (!pElectronShower.isAvailable() || pElectronShower.isNull())
        {
            mf::LogWarning("NeutrinoEnergyRecoPIDAlg") 
                << " Cannot access the electron shower which is needed for this energy reconstructio method.\n"
                << "Swapping to energy reconstruction method " << kProtonPionHadronic << " for this calculation." << std::endl;
            return this->CalculateNeutrinoEnergyPID(event, protonTracks, pionTracks);
        }


        Point_t vertex(pElectronShower->ShowerStart().X(), pElectronShower->ShowerStart().Y(), pElectronShower->ShowerStart().Z());

        const std::vector<art::Ptr<recob::Hit> > electronHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(pElectronShower, event, fShowerToHitLabel),2));
        const double electronEnergy(this->CalculateElectronEnergy(pElectronShower, event));
        const double electronMomentum = std::sqrt(electronEnergy*(electronEnergy + 2*kElectronMass));

        const Momentum4_t electron4Momentum(this->CalculateParticle4Momentum(kElectronMass, electronMomentum,
                    pElectronShower->Direction().X(), pElectronShower->Direction().Y(), pElectronShower->Direction().Z()));

        EnergyRecoInputHolder energyRecoInputHolder(vertex, electron4Momentum, 
                kElectronAndHadronic, kTrackMethodNotSet, kContainmentNotSet, fGradNuEHadEn, fIntNuEHadEn);

        return this->CalculateNeutrinoEnergyPID(electronHits, event, energyRecoInputHolder, protonTracks, pionTracks);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyPID(const std::vector<art::Ptr<recob::Hit> > &leptonHits, 
            const art::Event &event,
            const EnergyRecoInputHolder &energyRecoInputHolder,
            const std::vector<art::Ptr<recob::Track>> &protonTracks,
            const std::vector<art::Ptr<recob::Track>> &pionTracks)
    {


        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
        const double leptonObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, leptonHits));

        const int npions = pionTracks.size();
        const std::pair<double, double> protonsEnergyCharge(CalculateTracksObservedEnergyPID(event, protonTracks, pProton));
        const std::pair<double, double> pionsEnergyCharge(CalculateTracksObservedEnergyPID(event, pionTracks, pPion));

        const double protonsKin = protonsEnergyCharge.first;
        const double pionsKin = pionsEnergyCharge.first;

        const double protonsObservedCharge = protonsEnergyCharge.second;
        const double pionsObservedCharge = pionsEnergyCharge.second;

        const std::vector<art::Ptr<recob::Hit> > eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event, fHitLabel),2));
        const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits));

        const double hadronicObservedCharge(eventObservedCharge-(leptonObservedCharge+protonsObservedCharge + pionsObservedCharge));
        const double uncorrectedHadronicEnergy = this->CalculateEnergyFromCharge(hadronicObservedCharge) + protonsKin + pionsKin + kPionMass*npions;
        const double correctedHadronicEnergy(
                this->CalculateLinearlyCorrectedValue(uncorrectedHadronicEnergy,energyRecoInputHolder.fHadronicCorrectionGradient,
                    energyRecoInputHolder.fHadronicCorrectionIntercept));

        const double neutrinoEnergy(energyRecoInputHolder.fLeptonMomentum.E()+correctedHadronicEnergy);

        dune::EnergyRecoOutput output;
        output.recoMethodUsed = static_cast<int>(energyRecoInputHolder.fEnergyRecoMethod);
        output.fRecoVertex = energyRecoInputHolder.fVertex;
        output.fNuLorentzVector.SetE(neutrinoEnergy);
        output.fLepLorentzVector = energyRecoInputHolder.fLeptonMomentum;
        output.fHadLorentzVector.SetE(correctedHadronicEnergy);
        output.longestTrackContained = static_cast<int>(energyRecoInputHolder.fMuonContainmentStatus);
        output.trackMomMethod = static_cast<int>(energyRecoInputHolder.fMuonTrackMethod);
        return output;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyPID(const art::Event &event,
            const std::vector<art::Ptr<recob::Track>> &protonTracks,
            const std::vector<art::Ptr<recob::Track>> &pionTracks)
    {
        auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
        const double triggerTme(clockData.TriggerTime());

        const std::vector<art::Ptr<recob::Wire> > wires(dune_ana::DUNEAnaEventUtils::GetWires(event, fWireLabel));
        double wireCharge(0);

        for (unsigned int iWire = 0; iWire < wires.size(); ++iWire)
        {
            if (wireReadout.SignalType(wires[iWire]->Channel()) != geo::kCollection)
                continue;

            const recob::Wire::RegionsOfInterest_t& signalROI(wires[iWire]->SignalROI());
            for (const lar::sparse_vector<float>::datarange_t& range : signalROI.get_ranges())
            {
                const std::vector<float>& signal(range.data());
                const raw::TDCtick_t binFirstTickROI(range.begin_index());
                for (unsigned int iSignal = 0; iSignal < signal.size(); ++iSignal)
                    wireCharge += signal[iSignal]*dune_ana::DUNEAnaHitUtils::LifetimeCorrection(clockData, detProp, iSignal+binFirstTickROI, triggerTme);
            }
        }

        const int npions = pionTracks.size();
        const std::pair<double, double> protonsEnergyCharge(CalculateTracksObservedEnergyPID(event, protonTracks, pProton));
        const std::pair<double, double> pionsEnergyCharge(CalculateTracksObservedEnergyPID(event, pionTracks, pPion));

        const double protonsKin = protonsEnergyCharge.first;
        const double pionsKin = pionsEnergyCharge.first;

        const double protonsObservedCharge = protonsEnergyCharge.second;
        const double pionsObservedCharge = pionsEnergyCharge.second;
        const double totalEnergy(this->CalculateEnergyFromCharge(wireCharge - (protonsObservedCharge + pionsObservedCharge)) + protonsKin + pionsKin + kPionMass*npions);

        dune::EnergyRecoOutput output;
        output.recoMethodUsed = static_cast<int>(kAllCharges);
        output.fRecoVertex.SetXYZ(0.,0.,0.);
        output.fNuLorentzVector.SetXYZT(0.,0.,0.,totalEnergy);
        output.fLepLorentzVector.SetXYZT(0.,0.,0.,0.);
        output.fHadLorentzVector.SetXYZT(0.,0.,0.,0.);
        output.longestTrackContained = static_cast<int>(kContainmentNotSet);
        output.trackMomMethod = static_cast<int>(kTrackMethodNotSet);
        return output;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    std::pair<double, double> NeutrinoEnergyRecoAlg::CalculateTracksObservedEnergyPID(const art::Event &event,
            const std::vector<art::Ptr<recob::Track>> &tracks,
            const std::pair<int, double> &pParticle)
        {

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);

        trkf::TrackMomentumCalculator TrackMomCalc{0,3000}; // tracks should already be filtered, no reason to limit the lenght

        const size_t ntrks = tracks.size();
        const double pMass = pParticle.second;
        const int pPdgId = pParticle.first;

        double totalKinFromRange = 0;
        double totalObservedCharge = 0;

        for (size_t i = 0; i < ntrks; i++)
        {
            const double pMomRange = TrackMomCalc.GetTrackMomentum(tracks[i]->Length(), pPdgId);
            double pKinFromRange = 0;
            
            if (pMomRange > 0)
            {
                pKinFromRange = NeutrinoEnergyRecoAlg::KfromP(pMomRange,pMass);
            }
            else{
                std::cout << pKinFromRange << std::endl;
                // throw art::Exception(art::errors::LogicError) << "WHAT!?!?! \n";
            }
                std::cout << "... ok " << pKinFromRange << " " << pMomRange << " " << pPdgId << std::endl;
            const std::vector<art::Ptr<recob::Hit> > trkHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(tracks[i], event, fTrackToHitLabel),2));
            const double trkObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, trkHits));
            if (pParticle == pPion)
            {
                const double trkObservedEnergy(this->CalculateEnergyFromCharge(trkObservedCharge));
                if (trkObservedEnergy > pKinFromRange)
                    pKinFromRange = trkObservedEnergy;
            }
            totalObservedCharge=trkObservedCharge;
            totalKinFromRange+=pKinFromRange;
        }
        return std::make_pair(totalKinFromRange, totalObservedCharge);
    }

} //dune
