//STL
#include <string>
#include <iostream>
//ROOT
//ART
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
//LArSoft
#include "larreco/Calorimetry/CalorimetryAlg.h"

//DUNE

class NeutrinoEnergyRecoAlg 
{
    public:
        NeutrinoEnergyRecoAlg(const fhicl::ParameterSet &pset);

        double GetLifetimeCorrectedTotalHitCharge(const art::Event &event, const std::string &hitLabel);

    private:
        calo::CalorimetryAlg fCalorimetryAlg;

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
