//STL
//ROOT
//ART
//LArSoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
//DUNE
#include "dune/AnaUtils/DUNEAnaEventUtils.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"

#include "dune/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"

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

double NeutrinoEnergyRecoAlg::GetLifetimeCorrectedTotalHitCharge(const art::Event &event, const std::string &hitLabel)
{
    const detinfo::DetectorProperties *detprop(lar::providerFrom<detinfo::DetectorPropertiesService>());
    const double t0(detprop->TriggerOffset());


    std::vector<art::Ptr<recob::Hit> > hits(dune_ana::DUNEAnaEventUtils::GetHits(event, hitLabel));
    double lifetimeCorreectedEventCharge(0.);

    for (unsigned int i_hit = 0; i_hit < hits.size(); ++i_hit)
    {
        art::Ptr<recob::Hit> pHit(hits[i_hit]);

        if (2 == pHit->WireID().Plane)
            lifetimeCorreectedEventCharge += pHit->Integral() * fCalorimetryAlg.LifetimeCorrection(pHit->PeakTime(), t0);
        std::cout<<"Corrections: " << fCalorimetryAlg.LifetimeCorrection(pHit->PeakTime(), t0) << "  " << dune_ana::DUNEAnaHitUtils::GetLifetimeCorrection(pHit) << std::endl; 
    }

    return lifetimeCorreectedEventCharge;

}
