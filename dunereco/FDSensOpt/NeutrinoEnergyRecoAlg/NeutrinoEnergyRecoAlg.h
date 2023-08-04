/**
*  @file   dunereco/FDSensOpt/NeutrinoEnergyRecoAlg/NeutrinoEnergyRecoAlg.h
*
*  @brief  Header file for the neutrino energy reconstruction algorithm.  A heavily refactored version of Nick Grant's module
*
*  $Log: $
*/
#ifndef DUNE_NEUTRINO_ENERGY_RECO_ALG_H
#define DUNE_NEUTRINO_ENERGY_RECO_ALG_H

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
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

namespace dune
{
/**
 *
 * @brief NeutrinoEnergyRecoAlg class
 *
*/
class NeutrinoEnergyRecoAlg 
{
    public:
        /**
        * @brief  Constructor
        *
        * @param  pset the FCL parameter set
        * @param  trackLabel the track label
        * @param  showerLabel the shower label
        * @param  hitLabel the hit label
        * @param  wireLabel the wire label
        * @param  trackToHitLabel the associated track-to-hit label
        * @param  showerToHitLabel the associated shower-to-hit label
        * @param  hitToSpacePointLabel the associated hit-to-space point label
        */
        NeutrinoEnergyRecoAlg(const fhicl::ParameterSet &pset, const std::string &trackLabel, const std::string &showerLabel, 
            const std::string &hitLabel, const std::string wireLabel, const std::string &trackToHitLabel, 
            const std::string &showerToHitLabel, const std::string &hitToSpacePointLabel);

        /**
        * @brief  Calculates neutrino energy using a muon track (the muon track may be ignored if it isn't of a suitable quality)
        *
        * @param  pMuonTrack the muon track
        * @param  event the art event
        *
        * @return the neutrino energy summary object
        */
        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event);

        /**
        * @brief  Calculates neutrino energy using an electron shower(the electron may be ignored if it isn't of a suitable quality)
        *
        * @param  pElectronShower the electron shower
        * @param  event the art event
        *
        * @return the neutrino energy summary object
        */
        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event);

        /**
        * @brief  Calculates neutrino energy by summing wire charges
        *
        * @param  event the art event
        *
        * @return the neutrino energy summary object
        */
        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const art::Event &event);

        /**
        * @brief  Calculates neutrino energy explicitly using muon momentum by range
        *
        * @param  pMuonTrack the muon track
        * @param  event the art event
        *
        * @return the neutrino energy summary object
        */
        dune::EnergyRecoOutput CalculateNeutrinoEnergyViaMuonRanging(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event);

        /**
        * @brief  Calculates neutrino energy explicitly using muon multiple scattering
        *
        * @param  pMuonTrack the muon track
        * @param  event the art event
        *
        * @return the neutrino energy summary object
        */
        dune::EnergyRecoOutput CalculateNeutrinoEnergyViaMuonMCS(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event);

    private:

        typedef Position4_t Momentum4_t;

        double kMuonMass = 0.1056583745;                          ///< the muon mass (hardcoded unfortunately)
        double kElectronMass = 0.0005109989461;                   ///< the electron mass (hardcoded unfortunately);

        ///The energy reconstruction method
        enum EnergyRecoMethod
        {
            kRecoMethodNotSet = -1,                               ///< method not set
            kMuonAndHadronic = 1,                                 ///< muon momentum and hadronic deposited energy method
            kElectronAndHadronic,                                 ///< electron deposited energy and hadronic deposited energy method
            kAllCharges                                           ///< summed wire charges
        };
        ///The muon momentum reconstruction method
        enum MuonTrackMethod
        {
            kTrackMethodNotSet = -1,                              ///< method not set
            kMCS,                                                 ///< muon momentum by multiple scattering
            kContained                                            ///< muon momentum by range
        };
        ///The muon containment status
        enum MuonContainmentStatus
        {
            kContainmentNotSet = -1,                              ///< Containment not set
            kIsExiting,                                           ///< Muon exits
            kIsContained                                          ///< Muon is contained
        };

        /**
        *
        * @brief EnergyRecoInputHolder struct
        *
        */
        struct EnergyRecoInputHolder
        {
            /**
            * @brief  Constructor
            *
            * @param  vertex the reconstructed vertex
            * @param  leptonMomentum the reconstructed lepton momentum
            * @param  energyRecoMethod the neutrino energy reconstruction method
            * @param  muonTrackMethod the muon momentum reconstruction method
            * @param  muonContainmentStatus the muon containment status
            * @param  hadronicCorrectionGradient the linear correction gradient for reconstructed hadronic energy
            * @param  hadronicCorrectionIntercept the linear correction intercept for reconstructed hadronic energy
            */
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

            const Point_t fVertex;                                   ///< the reconstructed vertex
            const Momentum4_t fLeptonMomentum;                       ///< the reconstructed lepton four-momentum
            const EnergyRecoMethod fEnergyRecoMethod;                ///< the neutrino energy reconstruction method
            const MuonTrackMethod fMuonTrackMethod;                  ///< the muon momentum reconstruction method
            const MuonContainmentStatus fMuonContainmentStatus;      ///< the muon containment status
            const double fHadronicCorrectionGradient;                ///< the hadronic energy reconstruction linear correction gradient
            const double fHadronicCorrectionIntercept;               ///< the hadronic energy reconstruction linear correction intercept
        };

        /**
        * @brief  Calculates muon momentum by range
        *
        * @param  pMuonTrack the muon track
        *
        * @return the reconstructed muon momentum
        */
        double CalculateMuonMomentumByRange(const art::Ptr<recob::Track> pMuonTrack);

        /**
        * @brief  Calculates muon momentum by multiple coulomb scattering
        *
        * @param  pMuonTrack the muon track
        *
        * @return the reconstructed muon momentum
        */
        double CalculateMuonMomentumByMCS(const art::Ptr<recob::Track> pMuonTrack);

        /**
        * @brief  Calculates an electron shower's deposited energy by converting its deposited charge
        *
        * @param  pElectronShower the electron shower
        * @param  event the art event
        *
        * @return the reconstructed electron energy
        */
        double CalculateElectronEnergy(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event);

        /**
        * @brief  Converts deposited charge into energy by converting to number of electrons and correcting for average recombination
        *
        * @param  charge the deposited charge
        *
        * @return the reconstructed deposited energy
        */
        double CalculateEnergyFromCharge(const double charge);

        /**
        * @brief  Checks if a set of track hits are contained within a central volume of the detector 
        *
        * @param  hits the track hits
        *
        * @return an is contained bool
        */
        bool IsContained(const std::vector<art::Ptr<recob::Hit> > &hits, const art::Event &event);

        /**
        * @brief  Calculates a particle's four-momentum vector
        *
        * @param  mass the particle mass
        * @param  momentum the particle momentum
        * @param  directionX direction X component
        * @param  directionY direction X component
        * @param  directionZ direction X component
        *
        * @return the particle's four-momenutm vector
        */
        Momentum4_t CalculateParticle4Momentum(const double mass, const double momentum, 
            const double directionX, const double directionY, const double directionZ);

        /**
        * @brief  Linearly corrects a value
        *
        * @param  value the raw value
        * @param  correctionGradient the linear correction gradient
        * @param  correctionIntercept the linear correction intercept
        *
        * @return the linearly corrected value
        */
        double CalculateLinearlyCorrectedValue(const double value, const double correctionGradient,
            const double correctionIntercept);

        /**
        * @brief  Calculates the raw muon momentum by continuous-slowing-down approximation (CSDA) table
        *
        * @param  trkrange the muon track range (in cm)
        *
        * @return the uncorrected reconstructed muon momentum (in GeV)
        */
        double CalculateUncorrectedMuonMomentumByRange(const art::Ptr<recob::Track> &pMuonTrack);

        /**
        * @brief  Calculates the raw muon momentum by multiple coulomb scattering
        *
        * @param  pMuonTrack the muon track
        *
        * @return the uncorrected reconstructed muon momentum
        */
        double CalculateUncorrectedMuonMomentumByMCS(const art::Ptr<recob::Track> &pMuonTrack);

        /**
        * @brief  Calculates neutrino energy by summing hadronic deposited energy and lepton energy
        *
        * @param  leptonHits the lepton hits
        * @param  event the art event
        * @param  energyRecoInputHolder the holder object holding pre-calculated or pre-existing information
        *
        * @return the neutrino energy summary object
        */
        dune::EnergyRecoOutput CalculateNeutrinoEnergy(const std::vector<art::Ptr<recob::Hit> > &leptonHits, const art::Event &event, 
            const EnergyRecoInputHolder &energyRecoInputHolder);

        /**
        * @brief  Check's if a point is contained within a central detector volume
        *
        * @param  x the x component of the position
        * @param  y the y component of the position
        * @param  z the z component of the position
        *
        * @return an is contained bool
        */
        bool IsPointContained(const double x, const double y, const double z);

        calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm

        double fGradTrkMomRange;                                 ///< the correction gradient for muon momentum by range 
        double fIntTrkMomRange;                                  ///< the correction intercept for muon momentum by range
        double fGradTrkMomMCS;                                   ///< the correction gradient for muom momentum by MCS
        double fIntTrkMomMCS;                                    ///< the correction intercept for muon momentum by MCS
        double fGradNuMuHadEnCont;                               ///< the hqdronic energy correction gradient for numu+contained muon
        double fIntNuMuHadEnCont;                                ///< the hadronic energy correction intercept for numu+contained muon
        double fGradNuMuHadEnExit;                               ///< the hadronic energy correction gradient for numu+exiting muon
        double fIntNuMuHadEnExit;                                ///< the hadronic energy correction intercept for numu+exiting muon
        double fGradShwEnergy;                                   ///< the electron shower energy correction gradient
        double fIntShwEnergy;                                    ///< the electron shower energy correction intercept
        double fGradNuEHadEn;                                    ///< the hadronic energy correction gradient for nue
        double fIntNuEHadEn;                                     ///< the hadronic energy correction intercept for nue
        double fDistanceToWallThreshold;                         ///< the min distance from a detector wall to be considered contained
        double fMuonRangeToMCSThreshold;                         ///< the ratio threshold at which MCS is used for contained muons
        double fRecombFactor;                                    ///< the average reccombination factor

        std::string fTrackLabel;                                 ///< the track label
        std::string fShowerLabel;                                ///< the shower label
        std::string fHitLabel;                                   ///< the hit label
        std::string fWireLabel;                                  ///< the wire label
        std::string fTrackToHitLabel;                            ///< the associated track-to-hit label
        std::string fShowerToHitLabel;                           ///< the associated shower-to-hit label
        std::string fHitToSpacePointLabel;                       ///< the associated hit-to-space point label
};
} //namespace dune_ana
#endif //DUNE_NEUTRINO_ENERGY_RECO_ALG_H
