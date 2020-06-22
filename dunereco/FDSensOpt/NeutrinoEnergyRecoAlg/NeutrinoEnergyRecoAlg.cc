//STL
#include <limits>
#include <algorithm>
//ROOT
#include "TDatabasePDG.h"
//ART
#include "canvas/Utilities/Exception.h"
//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
//DUNE
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"
#include "dune/AnaUtils/DUNEAnaShowerUtils.h"
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"

#include "dune/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"

namespace dune
{
NeutrinoEnergyRecoAlg::NeutrinoEnergyRecoAlg(fhicl::ParameterSet const& pset, const std::string &trackLabel, 
    const std::string &showerLabel, const std::string &hitLabel, const std::string &trackToHitLabel, const std::string &showerToHitLabel, 
    const std::string &hitToSpacePointLabel) :
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
    fRecombFactor(pset.get<double>("RecombFactor")),
    fTrackLabel(trackLabel),
    fShowerLabel(showerLabel),
    fHitLabel(hitLabel),
    fTrackToHitLabel(trackToHitLabel),
    fShowerToHitLabel(showerToHitLabel),
    fHitToSpacePointLabel(hitToSpacePointLabel)
{
}

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event)
{
    Point_t vertex(pMuonTrack->Start().X(), pMuonTrack->Start().Y(), pMuonTrack->Start().Z());

    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel),2));
    bool isContained(IsContained(muonHits, event));
    const double uncorrectedMuonMomentumMCS(this->CalculateUncorrectedMuonMomentumByMCS(pMuonTrack));
    const double muonMomentumMCS(this->CalculateLinearlyCorrectedValue(uncorrectedMuonMomentumMCS, fGradTrkMomMCS, fIntTrkMomMCS));
    if (!isContained)
    {
        EnergyRecoInputHolder energyRecoInputHolder(vertex, 
            this->CalculateParticle4Momentum(13, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
            kMuonAndHadronic, kMCS, kIsExiting, fGradNuMuHadEnExit, fIntNuMuHadEnExit);
        return CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
    }
    else
    {
        const double muonMomentumRange(CalculateMuonMomentumByRange(pMuonTrack));
        if (uncorrectedMuonMomentumMCS > std::numeric_limits<double>::epsilon() 
            && muonMomentumRange/muonMomentumMCS - 0.7 < -1.*std::numeric_limits<double>::epsilon())
        {
            EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                this->CalculateParticle4Momentum(13, muonMomentumMCS, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                kMuonAndHadronic, kMCS, kIsContained, fGradNuMuHadEnExit, fIntNuMuHadEnExit);

            return CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
        }
        else
        {
            EnergyRecoInputHolder energyRecoInputHolder(vertex, 
                this->CalculateParticle4Momentum(13, muonMomentumRange, pMuonTrack->VertexDirection().X(), pMuonTrack->VertexDirection().Y(), pMuonTrack->VertexDirection().Z()), 
                kMuonAndHadronic, kContained, kIsContained, fGradNuMuHadEnCont, fIntNuMuHadEnCont);

            return CalculateNeutrinoEnergy(muonHits, event, energyRecoInputHolder);
        }
    }

    throw art::Exception(art::errors::LogicError) << "Unable to determine how to calculate neutrino energy using muon track! \n";
}

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event)
{
    Point_t vertex(pElectronShower->ShowerStart().X(), pElectronShower->ShowerStart().Y(), pElectronShower->ShowerStart().Z());

    const std::vector<art::Ptr<recob::Hit> > electronHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(pElectronShower, event, fShowerToHitLabel),2));
    const double electronEnergy(this->CalculateElectronEnergy(pElectronShower, event));
    const Momentum4_t electron4Momentum(this->CalculateParticle4Momentum(11, electronEnergy, 
        pElectronShower->Direction().X(), pElectronShower->Direction().Y(), pElectronShower->Direction().Z()));

    EnergyRecoInputHolder energyRecoInputHolder(vertex, electron4Momentum, 
    kElectronAndHadronic, kTrackMethodNotSet, kContainmentNotSet, fGradNuEHadEn, fIntNuEHadEn);

    return CalculateNeutrinoEnergy(electronHits, event, energyRecoInputHolder);
}

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Event &event)
{

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


NeutrinoEnergyRecoAlg::Momentum4_t NeutrinoEnergyRecoAlg::CalculateParticle4Momentum(const int pdg, const double momentum, 
    const double directionX, const double directionY, const double directionZ)
{
    static TDatabasePDG databasePDG;
    TParticlePDG* particlePDG(databasePDG.GetParticle(pdg));
    if (particlePDG) 
    {
        const double E(std::sqrt(momentum*momentum + particlePDG->Mass()*particlePDG->Mass()));
        const double pX(directionX * momentum);
        const double pY(directionY * momentum);
        const double pZ(directionZ * momentum);
        return Momentum4_t(pX, pY, pZ, E);
    }

    throw art::Exception(art::errors::NotFound) << "Unable to retrieve TParticlePDG for paticle with pdg " << pdg << "\n";
}

double NeutrinoEnergyRecoAlg::CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack)
{
    return this->CalculateLinearlyCorrectedValue(pMuonTrack->Length(), fGradTrkMomRange, fIntTrkMomRange);
}

double NeutrinoEnergyRecoAlg::CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack)
{
    const double uncorrectedMomentum(this->CalculateUncorrectedMuonMomentumByMCS(pMuonTrack));
    return this->CalculateLinearlyCorrectedValue(uncorrectedMomentum, fGradTrkMomMCS, fIntTrkMomMCS);
}

double NeutrinoEnergyRecoAlg::CalculateElectronEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event)
{
    const std::vector<art::Ptr<recob::Hit> > electronHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(pElectronShower, event, fShowerToHitLabel),2));
    const double electronObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(electronHits));
    const double uncorrectedElectronEnergy(this->CalculateEnergyFromCharge(electronObservedCharge));

    return this->CalculateLinearlyCorrectedValue(uncorrectedElectronEnergy, fGradShwEnergy, fIntShwEnergy);
}

double NeutrinoEnergyRecoAlg::CalculateEnergyFromCharge(const double charge)
{
    return fCalorimetryAlg.ElectronsFromADCArea(charge,2)*1./fRecombFactor/util::kGeVToElectrons;
}

bool NeutrinoEnergyRecoAlg::IsContained(const art::Ptr<recob::Track> pMuonTrack, const art::Event &event)
{
    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel),2));
    return this->IsContained(muonHits, event);
}

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

double NeutrinoEnergyRecoAlg::CalculateLinearlyCorrectedValue(const double value, const double correctionGradient,
    const double correctionIntercept)
{
    return (value - correctionIntercept) / correctionGradient;
}

double NeutrinoEnergyRecoAlg::CalculateUncorrectedMuonMomentumByMCS(const art::Ptr<recob::Track> &pMuonTrack)
{
    trkf::TrackMomentumCalculator TrackMomCalc;
    return (TrackMomCalc.GetMomentumMultiScatterChi2(pMuonTrack));
}

dune::EnergyRecoOutput NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const std::vector<art::Ptr<recob::Hit> > &leptonHits, 
    const art::Event &event, const EnergyRecoInputHolder &energyRecoInputHolder)
{
    const double leptonObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(leptonHits));

    const std::vector<art::Ptr<recob::Hit> > eventHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaEventUtils::GetHits(event, fHitLabel),2));
    const double eventObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(eventHits));

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

//Initalise (mode) (contained)
//Lepton momentum
//


bool NeutrinoEnergyRecoAlg::IsPointContained(const double x, const double y, const double z)
{
    double position[3] = {x,y,z};
    art::ServiceHandle<geo::Geometry> fGeometry;

    geo::TPCID tpcID(fGeometry->FindTPCAtPosition(position));
    if (!(fGeometry->HasTPC(tpcID)))
        return false;

    double minX(std::numeric_limits<double>::max());
    double maxX(std::numeric_limits<double>::min());
    double minY(std::numeric_limits<double>::max());
    double maxY(std::numeric_limits<double>::min());
    double minZ(std::numeric_limits<double>::max());
    double maxZ(std::numeric_limits<double>::min());
    for (unsigned int iCryostat = 0; iCryostat < fGeometry->Ncryostats(); ++iCryostat)
    {
        const geo::CryostatGeo& cryostat(fGeometry->Cryostat(iCryostat));
        for (unsigned int iTPC = 0; iTPC < cryostat.NTPC(); ++iTPC)
        {
            const geo::TPCGeo& tpc(cryostat.TPC(iTPC));
            minX = std::min(minX,tpc.MinX());
            maxX = std::max(maxX,tpc.MaxX());
            minY = std::min(minY,tpc.MinY());
            maxY = std::max(maxY,tpc.MaxY());
            minZ = std::min(minZ,tpc.MinZ());
            maxZ = std::max(maxZ,tpc.MaxZ());
        }
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

}

