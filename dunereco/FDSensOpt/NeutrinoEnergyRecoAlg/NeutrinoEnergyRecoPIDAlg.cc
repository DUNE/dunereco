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

    //------------------------------------------------------------------------------------------------------------------------------------------

    dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyPID(const std::vector<art::Ptr<recob::Hit> > &leptonHits, 
            const art::Event &event,
            const EnergyRecoInputHolder &energyRecoInputHolder)
    {


        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
        const double leptonObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, leptonHits));

        const int npions = NeutrinoEnergyRecoAlg::fPionTracks.size();
        const std::pair<double, double> protonsEnergyCharge(CalculateTracksObservedEnergyPID(event, NeutrinoEnergyRecoAlg::fProtonTracks, pProton));
        const std::pair<double, double> pionsEnergyCharge(CalculateTracksObservedEnergyPID(event, NeutrinoEnergyRecoAlg::fPionTracks, pPion));

        const double protonsKin = protonsEnergyCharge.first;
        const double pionsKin = pionsEnergyCharge.first;

        const double protonsObservedCharge = protonsEnergyCharge.second;
        const double pionsObservedCharge = pionsEnergyCharge.second;

        const std::vector<art::Ptr<recob::Hit>>
            eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event, fHitLabel),2));
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

    dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyPID(const double wireCharge,
            const art::Event &event)
    {

        const int npions = NeutrinoEnergyRecoAlg::fPionTracks.size();
        const std::pair<double, double> protonsEnergyCharge(CalculateTracksObservedEnergyPID(event, NeutrinoEnergyRecoAlg::fProtonTracks, pProton));
        const std::pair<double, double> pionsEnergyCharge(CalculateTracksObservedEnergyPID(event, NeutrinoEnergyRecoAlg::fPionTracks, pPion));

        const double protonsKin = protonsEnergyCharge.first;
        const double pionsKin = pionsEnergyCharge.first;

        const double protonsObservedCharge = protonsEnergyCharge.second;
        const double pionsObservedCharge = pionsEnergyCharge.second;
        const double totalEnergy(this->CalculateEnergyFromCharge(wireCharge - (protonsObservedCharge + pionsObservedCharge)) + protonsKin
                + pionsKin + kPionMass*npions);

        dune::EnergyRecoOutput output;
        output.recoMethodUsed = static_cast<int>(kProtonPionHadronic);
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

        trkf::TrackMomentumCalculator TrackMomCalc{0,3000}; // tracks should already be filtered, no reason to limit the length

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
                mf::LogError("NeutrinoEnergyRecoPIDAlg") << "Momentum by range was lower than one. Setting Kinetic energy to zero";
                pKinFromRange = 0;
            }
            const std::vector<art::Ptr<recob::Hit>>
                trkHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(tracks[i], event, fTrackToHitLabel),2));

            const double trkObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, trkHits));
            if (pParticle == pPion)
            {
                const double trkObservedEnergy(this->CalculateEnergyFromCharge(trkObservedCharge));
                if (trkObservedEnergy > pKinFromRange)
                    pKinFromRange = trkObservedEnergy;
            }
            totalObservedCharge+=trkObservedCharge;
            totalKinFromRange+=pKinFromRange;
        }
        return std::make_pair(totalKinFromRange, totalObservedCharge);
    }

} //dune
