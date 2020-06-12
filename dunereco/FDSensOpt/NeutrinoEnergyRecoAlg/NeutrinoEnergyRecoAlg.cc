//STL
//ROOT
//ART
//LArSoft
//DUNE
#include "dune/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h"

NeutrinoEnergyRecoAlg::NeutrinoEnergyRecoAlg(fhicl::ParameterSet const& pset)
    fGradTrkMomRange(pset.get<double>("GradTrkMomRange")),
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
