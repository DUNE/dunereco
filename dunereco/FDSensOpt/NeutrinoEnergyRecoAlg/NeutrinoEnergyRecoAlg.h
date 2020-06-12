//STL
#include <iostream>
//ROOT
//ART
#include "fhiclcpp/ParameterSet.h" 
//LArSoft
//DUNE

class NeutrinoEnergyRecoAlg 
{
    public:
        NeutrinoEnergyRecoAlg(fhicl::ParameterSet const& pset);
    private:
        double fGradTrkMomRange;
        double fIntTrkMomRange;
        double fGradTrkMomMCS;
        double fIntTrkMomMCS;
        double fGradNuMuHadEnCont;
        double fIntNuMuHadEnCont;
        double fGradNuMuHadEnExit;
        double fIntNuMuHadEnExit;
        double fGradShwEnergy;
        double fIntShwEnergy;
        double fGradNuEHadEn; 
        double fIntNuEHadEn;
        double fRecombFactor;
};
