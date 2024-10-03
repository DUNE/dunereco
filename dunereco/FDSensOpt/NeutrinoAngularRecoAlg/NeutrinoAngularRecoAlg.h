/**
*  @file   dunereco/FDSensOpt/NeutrinoAngularRecoAlg/NeutrinoAngularRecoAlg.h
*
*  @brief  Header file for the neutrino angular reconstruction algorithm.
*  Written by Pierre Granger (pierre.granger@cern.ch) & Henrique Souza (hvsouza@apc.in2p3.fr)
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
#include "lardataobj/AnalysisBase/Calorimetry.h"

//DUNE
#include "dunereco/FDSensOpt/FDSensOptData/AngularRecoOutput.h"

namespace dune
{
/**
 *
 * @brief NeutrinoAngularRecoAlg class
 *
*/

struct PolarFitOutput
{
    bool success;
    float pdrift;
    float pview;
    int nhits;
    int view;
};

class CalorimetricDirectionFitter {
    public:
        CalorimetricDirectionFitter(std::vector<PolarFitOutput>& observed_data, double thetaU, double thetaV);
        double Chi2(double const* x) const;
        int NViews() const;
    private:
        double const _thetaU;
        double const _thetaV;
        double _pU = 0, _pV = 0, _pW = 0;
        double _pxU = 0, _pxV = 0, _pxW = 0;
        int _nhitsU = 0, _nhitsV = 0, _nhitsW = 0;
        const std::vector<double> _calib_consts = {
            1.5e-5, //pxU
            1.5e-5, //pxV
            1e-5,   //pxW
            1.3e-5, //pU
            1.3e-5, //pV
            1e-5    //pW
        }; //TODO: Not ideal at all... I guess these could be obtained from calculations
};

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
        * @brief  Calculates neutrino angle assuming it follows the longest track 
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

        /**
        * @brief  Calculates neutrino angle using all the tracks and showers
        *
        * @param  pMuonTrack the muon track
        * @param  event the art event
        *
        * @return the neutrino direction summary object
        */
        dune::AngularRecoOutput CalculateNeutrinoAngle(const std::vector<art::Ptr<recob::Track>> &pTracks,
                                                       const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                                       const std::vector<art::Ptr<recob::Shower>> &pShowers,
                                                       const art::Event &event);

        /**
        * @brief  Calculates neutrino angle using all the tracks and showers, assuming the longest track is a muon
        *
        * @param  pMuonTrack the muon track
        * @param  event the art event
        *
        * @return the neutrino direction summary object
        */
        dune::AngularRecoOutput CalculateNeutrinoAngle(const art::Ptr<recob::Track> &pMuonTrack,
                                                       const std::vector<art::Ptr<recob::Track>> &pTracks,
                                                       const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                                       const std::vector<art::Ptr<recob::Shower>> &pShowers,
                                                       const art::Event &event);
                                                    
        /**
        * @brief  Calculates neutrino angle using only the hits information
        *
        * @param  event the art event
        * @param  vertex the reconstructed interaction vertex
        *
        * @return the neutrino direction summary object
        */
        dune::AngularRecoOutput CalculateNeutrinoAngle(const art::Event &event, const Point_t& vertex);


    private:


        ///The angular reconstruction method
        enum AngularRecoMethod
        {
            kRecoMethodNotSet = -1,                               ///< method not set
            kMuon = 1,                                 ///< muon momentum direction method
            kElectron,                                 ///< electron direction method
            kRecoParticles,                            ///< using all the reco particles       
            kMuonRecoParticles,                        ///< using all the reco particles nut assuming longest track is muon
            kHits                                      ///< using only the hits information      
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

        Momentum_t ComputeShowersMomentum(const std::vector<art::Ptr<recob::Shower>> &pShowers) const;
        Momentum_t ComputeTracksMomentum(const std::vector<art::Ptr<recob::Track>> &pTracks,
                                         const std::map<art::Ptr<recob::Track>, int> &tracksPID,
                                         const art::Event &event) const;
        float GetShowerEnergy(const art::Ptr<recob::Shower> &pShower) const;
        float GetTrackKE(const std::vector<art::Ptr<anab::Calorimetry>> &calos) const;
        bool IsTrackContained(const art::Ptr<recob::Track> &pTrack) const;

        PolarFitOutput FitViewHits(const art::Event &event, const geo::View_t &view, const std::vector<art::Ptr<recob::Hit>> &hits, const Point_t &vertex) const;
        PolarFitOutput PolarFit(const std::vector<float> &xvec, const std::vector<float> &yvec, const std::vector<float> &weights) const;
        geo::View_t GetTargetView(const art::Ptr<recob::Hit> &hit) const;
        float GetViewTheta(geo::View_t view) const;
        
        calo::CalorimetryAlg fCalorimetryAlg;                    ///< the calorimetry algorithm

        std::string fTrackLabel;                                 ///< the track label
        std::string fShowerLabel;                                ///< the shower label
        std::string fHitLabel;                                   ///< the hit label
        std::string fWireLabel;                                  ///< the wire label
        std::string fTrackToHitLabel;                            ///< the associated track-to-hit label
        std::string fShowerToHitLabel;                           ///< the associated shower-to-hit label
        std::string fHitToSpacePointLabel;                       ///< the associated hit-to-space point label
        std::string fCalorimetryLabel;                           ///< the calorimetry label
        float fDistanceToWallThreshold;                          ///< margin to consider wether a track is contained
        art::ServiceHandle<geo::Geometry const> fGeometry;       ///< handle to the geometry service

        const float fPION_MASS = 139.57; //MeV
};
} //namespace dune_ana
#endif //DUNE_NEUTRINO_ANGULAR_RECO_ALG_H
