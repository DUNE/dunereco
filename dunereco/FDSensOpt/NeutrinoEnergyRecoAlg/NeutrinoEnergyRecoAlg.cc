/**
*  @file   dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.cc

*  @brief  Implementation file for the neutrino energy reconstruction algorithm.  A heavily refactored version of Nick Grant's module
*
*  $Log: $
*/

//STL
#include <limits>
#include <algorithm>
//ROOT
//ART
#include "canvas/Utilities/Exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//LArSoft
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
NeutrinoEnergyRecoAlg::NeutrinoEnergyRecoAlg(fhicl::ParameterSet const& pset, const std::string &trackLabel, 
    const std::string &showerLabel, const std::string &hitLabel, const std::string wireLabel, 
    const std::string &trackToHitLabel, const std::string &showerToHitLabel, const std::string &hitToSpacePointLabel) :
    fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")),
    fGradTrkMomRange(pset.get<double>("GradTrkMomRange")),
    fIntTrkMomRange(pset.get<double>("IntTrkMomRange")),
    fGradTrkMomMCS(pset.get<double>("GradTrkMomMCS")),
    fIntTrkMomMCS(pset.get<double>("IntTrkMomMCS")),
    fGradNuMuHadEnCont(pset.get<double>("GradNuMuHadEnCont")),
    fIntNuMuHadEnCont(pset.get<double>("IntNuMuHadEnCont")),
    fGradNuMuHadEnExit(pset.get<double>("GradNuMuHadEnExit")),
    fIntNuMuHadEnExit(pset.get<double>("IntNuMuHadEnExit")),
    fGradShwEnergy(pset.get<double>("GradShwEnergy")),
    fIntShwEnergy(pset.get<double>("IntShwEnergy")),
    fGradNuEHadEn(pset.get<double>("GradNuEHadEn")),
    fIntNuEHadEn(pset.get<double>("IntNuEHadEn")),
    fDistanceToWallThreshold(pset.get<double>("DistanceToWallThreshold")),
    fMuonRangeToMCSThreshold(pset.get<double>("MuonRangeToMCSThreshold")),
    fRecombFactor(pset.get<double>("RecombFactor")),
    fTrackLabel(trackLabel),
    fShowerLabel(showerLabel),
    fHitLabel(hitLabel),
    fWireLabel(wireLabel),
    fTrackToHitLabel(trackToHitLabel),
    fShowerToHitLabel(showerToHitLabel),
    fHitToSpacePointLabel(hitToSpacePointLabel)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event)
{
    if (!pMuonTrack.isAvailable() || pMuonTrack.isNull())
    {
        mf::LogWarning("NeutrinoEnergyRecoAlg") << " Cannot access the muon track which is needed for this energy reconstructio method.\n"
        << "Swapping to energy reconstruction method " << kAllCharges << " for this calculation." << std::endl;
        return this->CalculateNeutrinoEnergy(event);
    }

    Point_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel),2));
    bool isContained(this->IsContained(muonHits, event));
    const double uncorrectedMuonMomentumMCS(this->CalculateUncorrectedMuonMomentumByMCS(pMuonTrack));
    const double muonMomentumMCS(this->CalculateLinearlyCorrectedValue(uncorrectedMuonMomentumMCS, fGradTrkMomMCS, fIntTrkMomMCS));
    if (!isContained)
    {
        if (uncorrectedMuonMomentumMCS > std::numeric_limits<double>::epsilon())
        {
            EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                this->CalculateParticle4Momentum(kMuonMass, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                kMuonAndHadronic, kMCS, kIsExiting, fGradNuMuHadEnExit, fIntNuMuHadEnExit);
            return this->CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
        }
        else
        {
            return this->CalculateNeutrinoEnergy(event);
        }
    }
    else
    {
        const double muonMomentumRange(CalculateMuonMomentumByRange(pMuonTrack));
        if (uncorrectedMuonMomentumMCS > std::numeric_limits<double>::epsilon() 
            && muonMomentumRange/muonMomentumMCS - fMuonRangeToMCSThreshold < -1.*std::numeric_limits<double>::epsilon())
        {
            EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                this->CalculateParticle4Momentum(kMuonMass, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                kMuonAndHadronic, kMCS, kIsContained, fGradNuMuHadEnExit, fIntNuMuHadEnExit);

            return this->CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
        }
        else
        {
            EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                this->CalculateParticle4Momentum(kMuonMass, muonMomentumRange, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                kMuonAndHadronic, kContained, kIsContained, fGradNuMuHadEnCont, fIntNuMuHadEnCont);

            return this->CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
        }
    }

    throw art::Exception(art::errors::LogicError) << "Unable to determine how to calculate neutrino energy using muon track! \n";
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Ptr<recob::Shower> &pElectronShower, 
    const art::Event &event)
{
    if (!pElectronShower.isAvailable() || pElectronShower.isNull())
    {
        mf::LogWarning("NeutrinoEnergyRecoAlg") 
        << " Cannot access the electron shower which is needed for this energy reconstructio method.\n"
        << "Swapping to energy reconstruction method " << kAllCharges << " for this calculation." << std::endl;
        return this->CalculateNeutrinoEnergy(event);
    }


    Point_t vertex(pElectronShower->ShowerStart().X(), pElectronShower->ShowerStart().Y(), pElectronShower->ShowerStart().Z());

    const std::vector<art::Ptr<recob::Hit> > electronHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(pElectronShower, event, fShowerToHitLabel),2));
    const double electronEnergy(this->CalculateElectronEnergy(pElectronShower, event));
    //ATTN yep this line is bugged.  It's deliberately bugged to maintain how the original code functioned
    //Because of the small electron mass vs total energy deposition, this bug will have a very very small effect on the 4vector
    const Momentum4_t electron4Momentum(this->CalculateParticle4Momentum(kElectronMass, electronEnergy, 
        pElectronShower->Direction().X(), pElectronShower->Direction().Y(), pElectronShower->Direction().Z()));

    EnergyRecoInputHolder energyRecoInputHolder(vertex, electron4Momentum, 
    kElectronAndHadronic, kTrackMethodNotSet, kContainmentNotSet, fGradNuEHadEn, fIntNuEHadEn);

    return this->CalculateNeutrinoEnergy(electronHits, event, energyRecoInputHolder);
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Event &event)
{
    art::ServiceHandle<geo::Geometry> fGeometry;
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
    const double triggerTme(clockData.TriggerTime());

    const std::vector<art::Ptr<recob::Wire> > wires(dune_ana::DUNEAnaEventUtils::GetWires(event, fWireLabel));
    double wireCharge(0);

    for (unsigned int iWire = 0; iWire < wires.size(); ++iWire)
    {
        if (fGeometry->SignalType(wires[iWire]->Channel()) != geo::kCollection)
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
    const double totalEnergy(this->CalculateEnergyFromCharge(wireCharge));

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

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyViaMuonRanging(const art::Ptr<recob::Track> &pMuonTrack,
    const art::Event &event)
{
    Point_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel),2));
    bool isContained(this->IsContained(muonHits, event));
    const double muonMomentumRange(this->CalculateMuonMomentumByRange(pMuonTrack));

    EnergyRecoInputHolder energyRecoInputHolder(vertex, 
        this->CalculateParticle4Momentum(kMuonMass, muonMomentumRange, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
        kMuonAndHadronic, kContained, static_cast<MuonContainmentStatus>(isContained), fGradNuMuHadEnCont, fIntNuMuHadEnCont);

    return this->CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergyViaMuonMCS(const art::Ptr<recob::Track> &pMuonTrack,
    const art::Event &event)
{
    Point_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel),2));
    bool isContained(this->IsContained(muonHits, event));
    const double muonMomentumMCS(this->CalculateMuonMomentumByMCS(pMuonTrack));

    EnergyRecoInputHolder energyRecoInputHolder(vertex, 
        this->CalculateParticle4Momentum(kMuonMass, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
        kMuonAndHadronic, kMCS, static_cast<MuonContainmentStatus>(isContained), fGradNuMuHadEnExit, fIntNuMuHadEnExit);

    return this->CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
}

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoEnergyRecoAlg::Momentum4_t NeutrinoEnergyRecoAlg::CalculateParticle4Momentum(const double mass, const double momentum, 
    const double directionX, const double directionY, const double directionZ)
{
    const double E(std::sqrt(momentum*momentum + mass*mass));
    const double pX(directionX * momentum);
    const double pY(directionY * momentum);
    const double pZ(directionZ * momentum);
    return Momentum4_t(pX, pY, pZ, E);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double NeutrinoEnergyRecoAlg::CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack)
{
    return this->CalculateLinearlyCorrectedValue(pMuonTrack->Length(), fGradTrkMomRange, fIntTrkMomRange);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double NeutrinoEnergyRecoAlg::CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack)
{
    const double uncorrectedMomentum(this->CalculateUncorrectedMuonMomentumByMCS(pMuonTrack));
    return this->CalculateLinearlyCorrectedValue(uncorrectedMomentum, fGradTrkMomMCS, fIntTrkMomMCS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double NeutrinoEnergyRecoAlg::CalculateElectronEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event)
{
    const std::vector<art::Ptr<recob::Hit> > electronHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(pElectronShower, event, fShowerToHitLabel),2));
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
    const double electronObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, electronHits));
    const double uncorrectedElectronEnergy(this->CalculateEnergyFromCharge(electronObservedCharge));

    return this->CalculateLinearlyCorrectedValue(uncorrectedElectronEnergy, fGradShwEnergy, fIntShwEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double NeutrinoEnergyRecoAlg::CalculateEnergyFromCharge(const double charge)
{
    return fCalorimetryAlg.ElectronsFromADCArea(charge,2)*1./fRecombFactor/util::kGeVToElectrons;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoEnergyRecoAlg::IsContained(const std::vector<art::Ptr<recob::Hit> > &hits, const art::Event &event)
{
    for (unsigned int iHit = 0; iHit < hits.size(); ++iHit)
    {
        const std::vector<art::Ptr<recob::SpacePoint> > spacePoints(dune_ana::DUNEAnaHitUtils::GetSpacePoints(hits[iHit], 
            event, fHitLabel, fHitToSpacePointLabel));
        for (unsigned int iSpacePoint = 0; iSpacePoint < spacePoints.size(); ++iSpacePoint)
        {
            const art::Ptr<recob::SpacePoint> spacePoint(spacePoints[iSpacePoint]);
            if (!(this->IsPointContained(spacePoint->XYZ()[0],spacePoint->XYZ()[1],spacePoint->XYZ()[2])))
                return false;
        }
    }
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double NeutrinoEnergyRecoAlg::CalculateLinearlyCorrectedValue(const double value, const double correctionGradient,
    const double correctionIntercept)
{
    return (value - correctionIntercept) / correctionGradient;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double NeutrinoEnergyRecoAlg::CalculateUncorrectedMuonMomentumByMCS(const art::Ptr<recob::Track> &pMuonTrack)
{
    trkf::TrackMomentumCalculator TrackMomCalc;
    return (TrackMomCalc.GetMomentumMultiScatterChi2(pMuonTrack, true));
}

//------------------------------------------------------------------------------------------------------------------------------------------

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const std::vector<art::Ptr<recob::Hit> > &leptonHits, 
    const art::Event &event, const EnergyRecoInputHolder &energyRecoInputHolder)
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event, clockData);
    const double leptonObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, leptonHits));

    const std::vector<art::Ptr<recob::Hit> > eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event, fHitLabel),2));
    const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockData, detProp, eventHits));

    const double hadronicObservedCharge(eventObservedCharge-leptonObservedCharge);
    const double uncorrectedHadronicEnergy(this->CalculateEnergyFromCharge(hadronicObservedCharge));
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

bool NeutrinoEnergyRecoAlg::IsPointContained(const double x, const double y, const double z)
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

} //dune
