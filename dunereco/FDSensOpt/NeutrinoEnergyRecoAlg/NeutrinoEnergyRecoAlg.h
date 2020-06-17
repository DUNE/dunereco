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
        NeutrinoEnergyRecoAlg(const fhicl::ParameterSet &pset, const std::string &trackLabel, 
            const std::string &hitLabel, const std::string &trackToHitLabel, const std::string &hitToSpacePointLabel);

        double CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event);
        
        double CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack);

        double CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack);

        bool IsContained(const art::Ptr<recob::Track> pMuonTrack, const art::Event &event);


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

        std::string fTrackLabel;
        std::string fHitLabel;
        std::string fTrackToHitLabel;
        std::string fHitToSpacePointLabel;
};
}
