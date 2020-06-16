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

namespace dune
{
class NeutrinoEnergyRecoAlg 
{
    public:
        NeutrinoEnergyRecoAlg(const fhicl::ParameterSet &pset);

        double CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event, 
            const std::string &hitLabel, const std::string &TrackToHitLabel, const std::string &hitToSpacePointLabel);

        
        double CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack);

        double CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack);


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
}
