//STL
#include <string>
#include <iostream>
//ROOT
#include "Math/GenVector/LorentzVector.h" 
//ART
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
//LArSoft
#include "larreco/Calorimetry/CalorimetryAlg.h"
//DUNE
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

namespace dune
{
class NeutrinoEnergyRecoAlg 
{
    public:
        NeutrinoEnergyRecoAlg(const fhicl::ParameterSet &pset, const std::string &trackLabel, const std::string &fShowerLabel, 
            const std::string &hitLabel, const std::string wireLabel, const std::string &trackToHitLabel, 
            const std::string &showerToHitLabel, const std::string &hitToSpacePointLabel);

        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event);

        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event);

        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const art::Event &event);

        double CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack);

        double CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack);

        double CalculateElectronEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event);

        double CalculateEnergyFromCharge(const double charge);

        bool IsContained(const art::Ptr<recob::Track> pMuonTrack, const art::Event &event);

        bool IsContained(const std::vector<art::Ptr<recob::Hit> > &hits, const art::Event &event);
       
    private:

        typedef Position4_t Momentum4_t;

        enum EnergyRecoMethod
        {
            kRecoMethodNotSet = -1,
            kMuonAndHadronic = 1,
            kElectronAndHadronic,
            kAllCharges
        };
        enum MuonTrackMethod
        {
            kTrackMethodNotSet = -1,
            kMCS,
            kContained
        };
        enum MuonContainmentStatus
        {
            kContainmentNotSet = -1,
            kIsExiting,
            kIsContained
        };

        struct EnergyRecoInputHolder
        {
            EnergyRecoInputHolder(const Point_t &vertex, const Momentum4_t &leptonMomentum, const EnergyRecoMethod energyRecoMethod,
                const MuonTrackMethod muonTrackMethod, const MuonContainmentStatus muonContainmentStatus, 
                const double hadronicCorrectionGradient, const double hadronicCorrectionIntercept) :
                fVertex(vertex),
                fLeptonMomentum(leptonMomentum),
                fEnergyRecoMethod(energyRecoMethod),
                fMuonTrackMethod(muonTrackMethod),
                fMuonContainmentStatus(muonContainmentStatus),
                fHadronicCorrectionGradient(hadronicCorrectionGradient),
                fHadronicCorrectionIntercept(hadronicCorrectionIntercept) {};

            const Point_t fVertex;
            const Momentum4_t fLeptonMomentum;
            const EnergyRecoMethod fEnergyRecoMethod;
            const MuonTrackMethod fMuonTrackMethod;
            const MuonContainmentStatus fMuonContainmentStatus;
            const double fHadronicCorrectionGradient;
            const double fHadronicCorrectionIntercept;
        };


        Momentum4_t CalculateParticle4Momentum(const int pdg, const double momentum, 
            const double directionX, const double directionY, const double directionZ);

        double CalculateLinearlyCorrectedValue(const double value, const double correctionGradient,
            const double correctionIntercept);

        double CalculateUncorrectedMuonMomentumByMCS(const art::Ptr<recob::Track> &pMuonTrack);

        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const std::vector<art::Ptr<recob::Hit> > &leptonHits, const art::Event &event, 
            const EnergyRecoInputHolder &energyRecoInputHolder);


        bool IsPointContained(const double x, const double y, const double z);

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
        double fDistanceToWallThreshold;
        double fRecombFactor;

        std::string fTrackLabel;
        std::string fShowerLabel;
        std::string fHitLabel;
        std::string fWireLabel;
        std::string fTrackToHitLabel;
        std::string fShowerToHitLabel;
        std::string fHitToSpacePointLabel;
};
}
