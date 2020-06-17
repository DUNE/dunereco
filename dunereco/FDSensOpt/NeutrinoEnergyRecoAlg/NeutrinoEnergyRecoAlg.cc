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
#include "dune/AnaUtils/DUNEAnaTrackUtils.h"

#include "dune/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"

namespace dune
{
NeutrinoEnergyRecoAlg::NeutrinoEnergyRecoAlg(fhicl::ParameterSet const& pset, const std::string &trackLabel, 
    const std::string &hitLabel, const std::string &trackToHitLabel, const std::string &hitToSpacePointLabel) :
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
    fHitLabel(hitLabel),
    fTrackToHitLabel(trackToHitLabel),
    fHitToSpacePointLabel(hitToSpacePointLabel)
{
    std::cout<<"This thing runs"<<std::endl;
}

double NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event)
{
    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel));
    const double muonHitCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(muonHits));

    const std::vector<art::Ptr<recob::Hit> > eventHits(dune_ana::DUNEAnaEventUtils::GetHits(event, fHitLabel));
    const double totalHitCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(eventHits));

    const double hadronicCharge(totalHitCharge-muonHitCharge);

    const double muonMomentum(this->CalculateAndPickMuonMomentum(pMuonTrack,event));

    return hadronicCharge+muonMomentum;
}

double NeutrinoEnergyRecoAlg::CalculateAndPickMuonMomentum(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event)
{
    bool isContained(this->IsContained(pMuonTrack,event));

    double momentumByMCS(this->CalculateMuonMomentumByMCS(pMuonTrack));
    if(!isContained)    
        return momentumByMCS;

    double momentumByRange(this->CalculateMuonMomentumByRange(pMuonTrack));

    if (momentumByRange/momentumByMCS - 0.7 > std::numeric_limits<double>::epsilon())
        return momentumByRange;
    else
        return momentumByMCS;

    throw art::Exception(art::errors::LogicError) << "Error when deciding which muon momentum to use! \n";
}

double NeutrinoEnergyRecoAlg::CalculateMuonTotalEnergyFromMomentum(const double muonMomentum)
{
    static TDatabasePDG pdg;
    TParticlePDG* muPDG(pdg.GetParticle(13));
    if (muPDG) return std::sqrt(muonMomentum*muonMomentum + muPDG->Mass()*muPDG->Mass());

    throw art::Exception(art::errors::NotFound) << "Unable to retrieve TParticlePDG for the muon \n";
}

double NeutrinoEnergyRecoAlg::CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack)
{
    return pMuonTrack->Length() - fIntTrkMomRange / fGradTrkMomRange;
}

double NeutrinoEnergyRecoAlg::CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack)
{
    trkf::TrackMomentumCalculator TrackMomCalc;
    return TrackMomCalc.GetMomentumMultiScatterChi2(pMuonTrack) - fIntTrkMomMCS / fGradTrkMomMCS;
}

bool NeutrinoEnergyRecoAlg::IsContained(const art::Ptr<recob::Track> pMuonTrack, const art::Event &event)
{
    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, fTrackToHitLabel));
    for (unsigned int iHit = 0; iHit < muonHits.size(); ++iHit)
    {
        const std::vector<art::Ptr<recob::SpacePoint> > muonSpacePoints(dune_ana::DUNEAnaHitUtils::GetSpacePoints(muonHits[iHit], 
            event, fHitLabel, fHitToSpacePointLabel));
        for (unsigned int iSpacePoint = 0; iSpacePoint < muonSpacePoints.size(); ++iSpacePoint)
        {
            const art::Ptr<recob::SpacePoint> spacePoint(muonSpacePoints[iSpacePoint]);
            if (!(this->IsPointContained(spacePoint->XYZ()[0],spacePoint->XYZ()[1],spacePoint->XYZ()[2])))
                return false;
        }
    }
    return true;
}

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

