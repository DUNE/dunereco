// ========================================================================================
// LowEUtils.h
// This utils library works with the SolarNuAna and LowEAna modules.
// It is used to find adjacent hits in time and space to create clusters in the context of
// Low Energy reconstruction and DUNE's solar neutrino analysis.
//
// @authors     : Sergio Manthey Corchado
// @created     : Nov, 2024
//=========================================================================================

#ifndef LOWEUTILS_H
#define LOWEUTILS_H

// LArSoft includes
#include "dunecore/ProducerUtils/ProducerUtils.h"
#include "dunereco/LowEUtils/LowECluster.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1I.h"
#include "TH1F.h"

// C++ Includes
#include <iostream>
#include <fcntl.h>
#include <string>
#include <memory>
#include <limits>
#include <vector>

using namespace producer;

namespace lowe
{
    class LowEUtils
    {
        public:
            explicit LowEUtils(fhicl::ParameterSet const &p);

            struct RawPerPlaneCluster
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
            
            struct RawSolarCluster
            {
                std::vector<float> Position;
                int MainID;
                int NHits;
                int TPC;
                int MainChannel;
                float TotalCharge;
                float AveragePeakTime;
                float Purity;
                float Completeness;
                std::vector<int> ClusterIdxVec;
            };
            
            void CalcAdjHits(
                std::vector<recob::Hit> MyVec,
                std::vector<std::vector<recob::Hit>> &Clusters,
                TH1I *MyHist,
                TH1F *ADCIntHist,
                art::Event const &evt,
                bool debug);
            
            void CalcAdjHits(
                std::vector<art::Ptr<recob::Hit>> MyVec,
                std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters,
                std::vector<std::vector<int>> &ClusterIdx,
                art::Event const &evt);
            
            void MakeClusterVector(
                std::vector<RawPerPlaneCluster> &ClusterVec,
                std::vector<std::vector<art::Ptr<recob::Hit>>> &Clusters,
                art::Event const &evt);

            void FillClusterVariables(
                std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
                std::vector<std::vector<int>> &ClNHits,
                std::vector<std::vector<float>> &ClT,
                std::vector<std::vector<float>> &ClCharge,
                art::Event const &evt,
                bool debug);            
            
            void FillClusterVariables(
                std::set<int> SignalTrackIDs,
                std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
                std::vector<std::vector<int>> &ClMainID,
                std::vector<std::vector<int>> &ClNHits,
                std::vector<std::vector<int>> &ClTPC,
                std::vector<std::vector<int>> &ClChannel,
                std::vector<std::vector<float>> &ClT,
                std::vector<std::vector<float>> &ClY,
                std::vector<std::vector<float>> &ClZ,
                std::vector<std::vector<float>> &ClDir,
                std::vector<std::vector<float>> &ClCharge,
                std::vector<std::vector<float>> &ClPurity,
                std::vector<std::vector<float>> &ClCompleteness,
                detinfo::DetectorClocksData const &clockData,
                bool debug);
            
            void FillClusterHitVectors(
                std::vector<recob::Hit> Cluster,
                std::vector<int> &TPC,
                std::vector<int> &Channel,
                std::vector<double> &Charge,
                std::vector<double> &Time,
                double &DeltaTime,
                double &SigmaTima,
                std::vector<double> &Y,
                std::vector<double> &Z,
                std::vector<double> &Dir,
                detinfo::DetectorClocksData clockData);

            void MatchClusters(
                std::vector<std::vector<std::vector<recob::Hit>>> &MatchedClusters,
                std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
                std::vector<std::vector<int>> &ClNHits,
                std::vector<std::vector<float>> &ClT,
                std::vector<std::vector<float>> &ClCharge,
                art::Event const &evt,
                bool debug = false);
            
            void MatchClusters(
                std::set<int> SignalTrackIDs,
                std::vector<std::vector<int>> &MatchedClustersIdx,
                std::vector<std::vector<std::vector<recob::Hit>>> &MatchedClusters,
                std::vector<std::vector<int>> ClustersIdx,
                std::vector<std::vector<std::vector<recob::Hit>>> Clusters,
                std::vector<std::vector<int>> &ClMainID,
                std::vector<std::vector<int>> &ClNHits,
                std::vector<std::vector<int>> &ClTPC,
                std::vector<std::vector<int>> &ClChannel,
                std::vector<std::vector<float>> &ClT,
                std::vector<std::vector<float>> &ClY,
                std::vector<std::vector<float>> &ClZ,
                std::vector<std::vector<float>> &ClDir,
                std::vector<std::vector<float>> &ClCharge,
                std::vector<std::vector<float>> &ClPurity,
                std::vector<std::vector<float>> &ClCompleteness,
                detinfo::DetectorClocksData const &clockData,
                bool debug = false);
            
            static std::vector<double> ComputeRecoY(
                int Event,
                std::vector<int> &IndTPC,
                std::vector<double> &Z,
                std::vector<double>  &T,
                std::vector<double> &IndZ,
                std::vector<double> &IndY,
                std::vector<double> &IndT,
                std::vector<double> &IndDir,
                bool debug = false);
            
            double STD(
                const std::vector<double>& Vec);

            void FindPrimaryClusters(
                const std::vector<art::Ptr<solar::LowECluster>> &SolarClusterVector,
                std::vector<bool> &EventCandidateFound,
                std::vector<std::vector<art::Ptr<solar::LowECluster>>> &EventCandidateVector,
                std::vector<std::vector<int>> &EventCandidateIdx,
                const detinfo::DetectorClocksData &clockData,
                const art::Event &evt);

            int MatchPDSFlash(
                const std::vector<art::Ptr<solar::LowECluster>> &SolarClusterVector,
                const std::vector<art::Ptr<recob::OpFlash>> &PDSFlashes,
                const detinfo::DetectorClocksData &clockData,
                const art::Event &evt,
                bool debug = false);

            bool SelectPDSFlashPE(
                const float &TPCDriftTime,
                const float &MatchedDriftTime,
                const float &ClusterCharge,
                const float &OpFlashPE);

            bool CutPDSFlashMinPE(
                const float &TPCDriftTime,
                const float &MatchedDriftTime,
                const float &ClusterCharge,
                const float &OpFlashPE);

            bool CutPDSFlashMaxPE(
                const float &TPCDriftTime,
                const float &MatchedDriftTime,
                const float &ClusterCharge,
                const float &OpFlashPE);

            void GetLightMapParameters(
                const std::string &LightMapType,
                const float &ClusterCharge,
                float &a,
                float &b,
                float &c);

            bool SelectPDSFlash(
                const float &TPCDriftTime,
                const float &ClusterTime,
                const float &ClusterCharge,
                const float &RefOpFlashTime,
                const float &RefOpFlashPe,
                const float &OpFlashTime,
                const float &OpFlashPE);

        // Declare member data here.
        private:
            // From fhicl configuration
            const std::string fHitLabel;
            const double fClusterAlgoTime;
            const int fClusterAlgoAdjChannel;
            const std::string fClusterChargeVariable;
            const int fClusterMatchNHit;
            const double fClusterMatchCharge;
            const double fClusterInd0MatchTime;
            const double fClusterInd1MatchTime;
            const double fClusterMatchTime;
            const bool fClusterPreselectionSignal;             // If true, only consider clusters with signal track IDs
            const int fClusterPreselectionNHits;
            const float fAdjClusterRad;                        // Radius in cm to search for adjacent clusters
            const bool fAdjClusterSingleMatch;                 // If true, match adjacent clusters to the first primary cluster only
            const int fAdjOpFlashMinNHitCut;                   // Minimum number of hits for adjacent flash
            const double fAdjOpFlashX;                         // Maximum X distance for adjacent flash matching [cm]
            const double fAdjOpFlashY;                         // Maximum Y distance for adjacent flash matching [cm]
            const double fAdjOpFlashZ;                         // Maximum Z distance for adjacent flash matching [cm]
            const double fAdjOpFlashMaxPERatioCut;             // Maximum photoelectrons ratio for adjacent flash
            const double fAdjOpFlashMinPECut;                  // Minimum photoelectrons for adjacent flash
            const double fAdjOpFlashMaxPECut;                  // Maximum photoelectrons for adjacent flash
            const bool fAdjOpFlashMembraneProjection;          // If true, project the TPC reco onto the membrane
            const bool fAdjOpFlashEndCapProjection;            // If true, project the TPC reco onto the end cap
            const double fAdjOpFlashMinPEAttenuation;          // Attenuation factor for minimum PE cut based on drift time [us]
            const double fAdjOpFlashMaxPEAttenuation;          // Attenuation factor for maximum PE cut based on drift time [us]
            const std::string fAdjOpFlashMinPEAttenuate;       // Type of attenuation for minimum PE cut ("light_map", "asymptotic", "linear" or "flat")
            const std::string fAdjOpFlashMaxPEAttenuate;       // Type of attenuation for maximum PE cut ("light_map", "asymptotic", "linear" or "flat")
            const int fAdjOpFlashMinPEAttenuationStrength;     // Strength of the asymptotic attenuation for minimum PE cut (in powers of 10)
            const int fAdjOpFlashMaxPEAttenuationStrength;     // Strength of the asymptotic attenuation for maximum PE cut (in powers of 10)
            const std::vector<std::pair<std::string, std::vector<double>>> fAdjOpFlashMinPELightMap; // Light map file and histogram name for light map attenuation
            const std::vector<std::pair<std::string, std::vector<double>>> fAdjOpFlashMaxPELightMap; // Light map file and histogram name for light map attenuation
            const std::vector<std::pair<std::string, std::vector<double>>> fAdjOpFlashPELightMap;    // Light map file and histogram name for PE attenuation
            const std::string fFlashMatchBy;                    // Method to match flashes ("maximum" or "light_map")
            std::unique_ptr<producer::ProducerUtils> producer; // Pointer to the ProducerUtils instance
    };
} // namespace lowe
#endif