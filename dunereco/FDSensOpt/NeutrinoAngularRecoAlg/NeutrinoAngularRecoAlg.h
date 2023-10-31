/**
*  @file   dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h
*
*  @brief  Header file for the neutrino angular reconstruction algorithm.
*
*  $Log: $
*/
#ifndef DUNE_NEUTRINO_ANGULAR_RECO_ALG_H
#define DUNE_NEUTRINO_ANGULAR_RECO_ALG_H

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
#include "dunereco/FDSensOpt/FDSensOptData/AngularRecoOutput.h"

namespace dune
{
/**
 *
 * @brief NeutrinoAngularRecoAlg class
 *
*/
class NeutrinoAngularRecoAlg 
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
        NeutrinoAngularRecoAlg(const fhicl::ParameterSet &pset, const std::string &trackLabel, const std::string &showerLabel, 
            const std::string &hitLabel, const std::string wireLabel, const std::string &trackToHitLabel, 
            const std::string &showerToHitLabel, const std::string &hitToSpacePointLabel);

        /**
        * @brief  Calculates neutrino angle assuming it follow the longest track 
        *
        * @param  pMuonTrack the muon track
        * @param  event the art event
        *
        * @return the neutrino direction summary object
        */
        dune::AngularRecoOutput CalculateNeutrinoAngle(const art::Ptr<recob::Track> &pMuonTrack, const art::Event &event);

        /**
        * @brief  Calculates neutrino angle using an electron shower 
        *
        * @param  pElectronShower the electron shower
        * @param  event the art event
        *
        * @return the neutrino direction summary object
        */
        dune::AngularRecoOutput CalculateNeutrinoAngle(const art::Ptr<recob::Shower> &pElectronShower, const art::Event &event);




    private:


        ///The angular reconstruction method
        enum AngularRecoMethod
        {
            kRecoMethodNotSet = -1,                               ///< method not set
            kMuon = 1,                                 ///< muon momentum direction method
            kElectron,                                 ///< electron direction method
        };

        /**
        *
        * @brief AngularRecoInputHolder struct
        *
        */
        struct AngularRecoInputHolder
        {
            /**
            * @brief  Constructor
            *
            * @param  vertex the reconstructed vertex
            * @param  nuDirection the reconstructed nu direction
            * @param  angularRecoMethod the neutrino direction reconstruction method
            */
            AngularRecoInputHolder(const Point_t &vertex, const Direction_t &nuDirection, const AngularRecoMethod angularRecoMethod) :
                fVertex(vertex),
                fNuDirection(nuDirection),
                fAngularRecoMethod(angularRecoMethod) {};

            const Point_t fVertex;                                   ///< the reconstructed vertex
            const Direction_t fNuDirection;                       ///< the reconstructed neutrino direction
            const AngularRecoMethod fAngularRecoMethod;                ///< the neutrino angle reconstruction method
        };

        /**
        * @brief  Return neutrino angle reconstruction output
        *
        * @return AngularRecoOutput 
        */
        dune::AngularRecoOutput ReturnNeutrinoAngle(const AngularRecoInputHolder &angularRecoInputHolder);

        calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm

        std::string fTrackLabel;                                 ///< the track label
        std::string fShowerLabel;                                ///< the shower label
        std::string fHitLabel;                                   ///< the hit label
        std::string fWireLabel;                                  ///< the wire label
        std::string fTrackToHitLabel;                            ///< the associated track-to-hit label
        std::string fShowerToHitLabel;                           ///< the associated shower-to-hit label
        std::string fHitToSpacePointLabel;                       ///< the associated hit-to-space point label
};
} //namespace dune_ana
#endif //DUNE_NEUTRINO_ANGULAR_RECO_ALG_H
