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
#include "lardataobj/RecoBase/Cluster.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1I.h"
#include "TH1F.h"

namespace solar
{
    class LowEUtils
    {
    public:
        struct LowEClusterInfo
        {
            float StartWire;
            float SigmaStartWire;
            float StartTick;
            float SigmaStartTick;
            float StartCharge;
            float StartAngle;
            float StartOpeningAngle;
            float EndWire;
            float SigmaEndWire;
            float EndTick;
            float SigmaEndTick;
            float EndCharge;
            float EndAngle;
            float EndOpeningAngle;
            float Integral;
            float IntegralStdDev;
            float SummedADC;
            float SummedADCstdDev;
            int NHit;
            float MultipleHitDensity;
            float Width;
            int ID;
            geo::View_t View;
            geo::PlaneID Plane;
            float WireCoord;
            float SigmaWireCoord;
            float TickCoord;
            float SigmaTickCoord;
            float IntegralAverage;
            float SummedADCaverage;
            float Charge;
            float ChargeStdDev;
            float ChargeAverage;
        };
        struct SolarCluster
        {
            int NHit;
            float X;
            float Y;
            float YWidth;
            float Z;
            float ZWidth;
            float Charge;
            float Integral;
            float SummedADC;
            int ID;
            geo::View_t View;
            geo::PlaneID Plane;
        };
        explicit LowEUtils(fhicl::ParameterSet const &p);
        void CalcAdjHits(std::vector<recob::Hit> MyVec, std::vector<std::vector<recob::Hit>> &Clusters, TH1I *MyHist, TH1F *ADCIntHist, bool HeavDebug);
        void CalcAdjHits(std::vector<art::Ptr<recob::Hit>> MyVec, std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters, std::vector<std::vector<int>> &ClusterIdx);
        void MakeClusterVector(std::vector<LowEClusterInfo> &ClusterVec, std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters, art::Event const &evt);
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
        const std::string fHitLabel;
        const std::string fGeometry;
        const double fDetectorSizeX;
        const double fClusterAlgoTime;
        const int fClusterAlgoAdjChannel;
        const std::string fClusterChargeVariable;
        const int fClusterMatchNHit;
        const double fClusterMatchCharge;
        const double fClusterInd0MatchTime;
        const double fClusterInd1MatchTime;
        const double fClusterMatchTime;
    };
}
#endif