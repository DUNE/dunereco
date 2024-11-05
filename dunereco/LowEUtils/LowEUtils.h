// ========================================================================================
// LowEUtils.h
// This utils library works with the SolarNuAna and LowEAna modules.
// It is used to find adjacent hits in time and space to create clusters in the context of
// Low Energy reconstruction and DUNE's solar neutrino analysis.
//
// @authors     : Sergio Manthey Corchado
// @created     : Nov, 2024
//=========================================================================================

#ifndef LowETool_h
#define LowETool_h

#include <iostream>
#include <vector>
#include <fcntl.h>

#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1I.h"
#include "TH1F.h"

namespace solar
{
    class LowEUtils
    {
    public:
        explicit LowEUtils(fhicl::ParameterSet const &p);
        void CalcAdjHits(std::vector<recob::Hit> MyVec, std::vector<std::vector<recob::Hit>> &Clusters, TH1I *MyHist, TH1F *ADCIntHist, bool HeavDebug);
        void FillClusterVariables(
            std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
            std::vector<std::vector<int>> &ClNHits,
            std::vector<std::vector<float>> &ClT,
            std::vector<std::vector<float>> &ClCharge,
            bool HeavDebug);
        void MatchClusters(
            std::vector<std::vector<std::vector<recob::Hit>>> MatchedClusters,
            std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
            std::vector<std::vector<int>> &ClNHits,
            std::vector<std::vector<float>> &ClT,
            std::vector<std::vector<float>> &ClCharge,
            bool HeavDebug);
        static std::vector<double> ComputeRecoY(
            int Event,
            std::vector<int> &IndTPC,
            std::vector<double> &Z,
            std::vector<double> &T,
            std::vector<double> &IndZ,
            std::vector<double> &IndY,
            std::vector<double> &IndT,
            std::vector<double> &IndDir,
            bool HeavDebug);
    private:
        // From fhicl configuration
        const double fClusterAlgoTime;
        const int fClusterAlgoAdjChannel;
        const int fClusterMatchNHit;
        const double fClusterMatchCharge;
        const double fClusterInd0MatchTime;
        const double fClusterInd1MatchTime;
        const double fClusterMatchTime;
    };
}
#endif