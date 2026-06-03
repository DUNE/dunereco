#ifndef PANDRIZZLEALG_H_SEEN
#define PANDRIZZLEALG_H_SEEN

///////////////////////////////////////////////
// PandrizzleAlg.h (D. Brailsford)
//
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/PtrVector.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "Pandora/PandoraInternal.h"
// c++
#include <map>
#include <memory>
#include <vector>

// ROOT
#include "TMVA/Reader.h"
#include "TTree.h"

namespace FDSelection
{
  class PandrizzleAlg 
  {
    public:  
      enum Vars{
        kEvalRatio = 0,
        kConcentration,
        kCoreHaloRatio,
        kConicalness,
        kdEdxBestPlane,
        kDisplacement,
        kDCA,
        kWideness,
        kEnergyDensity,
        kPathwayLengthMin,
        kMaxShowerStartPathwayScatteringAngle2D,
        kMaxNPostShowerStartHits,
        kMaxPostShowerStartScatterAngle,
        kMaxPostShowerStartNuVertexEnergyAsymmetry,
        kMaxPostShowerStartShowerStartEnergyAsymmetry,
        kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance,
        kMinPostShowerStartShowerStartMoliereRadius,
        kMaxPostShowerStartOpeningAngle,
        kMaxFoundHitRatio,
        kMaxInitialGapSize,
        kMinLargestProjectedGapSize,
        kNViewsWithAmbiguousHits,
        kAmbiguousHitMaxUnaccountedEnergy,
        kModularShowerPathwayLengthMin,
        kModularShowerMaxNShowerHits,
        kModularShowerMaxNuVertexChargeWeightedMeanRadialDistance,
        kBDTMethod,
        kEnhancedPandrizzleScore,
        kBackupPandrizzleScore,
        kTerminatingValue //terminates the enum and not an actual variable
      };

      using InputVarsToReader = std::map<Vars, std::unique_ptr<Float_t>>;

      class Record {
        public:
          Record(const InputVarsToReader &inputVars, const Float_t mvaScore, const bool isFilled);

          Float_t GetVar(const FDSelection::PandrizzleAlg::Vars var);
          bool IsFilled();
          Float_t GetMVAScore();

        private:
          using InputVars = std::map<Vars, Float_t>;
          InputVars fInputs;
          Float_t fMVAScore;
          bool fIsFilled;
      };

      PandrizzleAlg(const fhicl::ParameterSet& pset);

      void Run(const art::Event& evt);
      Record RunPID(const art::Ptr<recob::Shower> pShower, const art::Event& evt);

      private:
      void InitialiseTrees();
      void BookTreeInt(TTree *tree, std::string branch_name);
      void BookTreeFloat(TTree *tree, std::string branch_name);
      void BookTreeBool(TTree *tree, std::string branch_name);
      void ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
      void FillTruthInfo(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
      void FillMVAInfo(const art::Ptr<recob::Shower> pShower, const art::Event& evt);
      void FillPandrizzleInfo(const art::Ptr<recob::Shower> pShower, const art::Event& evt);
      void FillEnhancedPandrizzleInfo(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
      void FillBackupPandrizzleInfo(const art::Ptr<recob::PFParticle> pfp, const art::Ptr<recob::Shower> pShower, const art::Event& evt);
      void GetPathwayVariables(art::Ptr<recob::Track> &trackStub, float &modularShowerPathwayLengthMin, float &modularShowerPathwayKink3D);
      float GetLargest3DPathwayKink(art::Ptr<recob::Track> &trackStub);
      void SetInitialRegionVariables(TVector3 &nuVertexPosition, art::Ptr<recob::Track> &trackStub);
      void GetShowerRegionVariables(const TVector3 &nuVertex, const art::Ptr<recob::PFParticle> &pfp, art::Ptr<recob::Track> &trackStub, float &modularShowerMaxNShowerHits,
        float &modularShowerMaxFoundHitRatio, float &modularShowerMaxOpeningAngle, float &modularShowerMaxScatterAngle, float &modularShowerMaxNuVertexChargeAsymmetry,
        float &modularShowerMaxShowerStartChargeAsymmetry, float &modularShowerMaxNuVertexChargeWeightedMeanRadialDistance, float &modularShowerMinShowerStartMoliereRadius, 
        const art::Event& evt);
      float GetShowerOpeningAngle(const geo::Vector_t &showerStart, const pandora::CartesianVector &fittedShowerDirection,
        pandora::CartesianPointVector &cartesianPointVector);
      void GetShowerChargeDistributionVariables(const pandora::CartesianVector nuVertexPosition, const pandora::CartesianVector connectionPathwayDirection,
        const pandora::CartesianVector &fittedShowerStart, const pandora::CartesianVector &fittedShowerDirection, std::vector<art::Ptr<recob::Hit>> &showerHits, float &nuVertexChargeAsymmetry,
        float &showerStartChargeAsymmetry, float &nuVertexChargeWeightedMeanRadialDistance, float &showerStartMoliereRadius, const art::Event& evt);
      void FillTree();
      void ResetTreeVariables();
      double YZtoU(double y, double z);
      double YZtoV(double y, double z);
      double YZtoW(double y, double z);
      const geo::Vector_t GetPandoraHitPosition(art::Ptr<recob::Hit> pHit, const art::Event& evt);
      const geo::View_t GetPandoraHitView(art::Ptr<recob::Hit> pHit);

      Float_t* GetVarPtr(const FDSelection::PandrizzleAlg::Vars var);
      void SetVar(const FDSelection::PandrizzleAlg::Vars var, const Float_t value);
      Record ReturnEmptyRecord();

      std::string fRecoModuleLabel;
      std::string fShowerModuleLabel;
      std::string fPIDModuleLabel;
      std::string fPandrizzleWeightFileName;
      std::string fEnhancedPandrizzleWeightFileName;

      TMVA::Reader fReader;
      TMVA::Reader fEnhancedReader;
      InputVarsToReader fInputsToReader;

      // Tree things
      bool fMakeSelectionTrainingTrees;
      TTree *fSignalShowerTree;
      TTree *fBackgroundShowerTree;

      // Configuration for BDTs
      bool fUseConcentration;
      bool fUseDisplacement;
      bool fUseDCA;
      bool fUseBDTVariables;
      bool fUseModularShowerVariables;
      int fEnhancedPandrizzleHitCut;
      int fBackupPandrizzleHitCut;

      //Hold all variables to be handed to the trees
      struct VarHolder
      {
        std::map<std::string, int> IntVars;
        std::map<std::string, float> FloatVars;
        std::map<std::string, bool> BoolVars;

        std::vector<float> m_trajPositionX;
        std::vector<float> m_trajPositionY;
        std::vector<float> m_trajPositionZ;
      };

      VarHolder fVarHolder;

      //services
      art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
      art::ServiceHandle<art::TFileService> tfs;
  };
}

Float_t FDSelection::PandrizzleAlg::Record::GetVar(const FDSelection::PandrizzleAlg::Vars var)
{
  return (fInputs.at(var));
}

bool FDSelection::PandrizzleAlg::Record::IsFilled()
{
  return fIsFilled;
}

Float_t FDSelection::PandrizzleAlg::Record::GetMVAScore()
{
  return fMVAScore;
}

Float_t* FDSelection::PandrizzleAlg::GetVarPtr(const FDSelection::PandrizzleAlg::Vars var)
{
  return fInputsToReader.at(var).get();
}

void FDSelection::PandrizzleAlg::SetVar(const FDSelection::PandrizzleAlg::Vars var, const Float_t value)
{
  std::map<Vars, std::unique_ptr<Float_t> >::iterator itr(fInputsToReader.find(var));
  if (itr != fInputsToReader.end())
  {
    *(itr->second.get()) = value;
  }

  return;
}

#endif
