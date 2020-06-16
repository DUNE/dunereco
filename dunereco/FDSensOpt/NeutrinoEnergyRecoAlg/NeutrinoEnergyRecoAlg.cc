//STL
//ROOT
//ART
//LArSoft
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
NeutrinoEnergyRecoAlg::NeutrinoEnergyRecoAlg(fhicl::ParameterSet const& pset) :
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
    fRecombFactor(pset.get<double>("RecombFactor"))
{
    std::cout<<"This thing runs"<<std::endl;
}

double NeutrinoEnergyRecoAlg::CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event, 
    const std::string &hitLabel, const std::string &trackToHitLabel, const std::string &hitToSpacePointLabel)
{
    const std::vector<art::Ptr<recob::Hit> > muonHits(dune_ana::DUNEAnaTrackUtils::GetHits(pMuonTrack, event, trackToHitLabel));
    const double muonHitCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(muonHits));

    const std::vector<art::Ptr<recob::Hit> > eventHits(dune_ana::DUNEAnaEventUtils::GetHits(event, hitLabel));
    const double totalHitCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(eventHits));

    const double hadronicCharge(totalHitCharge-muonHitCharge);

    return hadronicCharge;
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

}

