#ifndef PANDIZZLEALG_H_SEEN
#define PANDIZZLEALG_H_SEEN

///////////////////////////////////////////////
// PandizzleAlg.h
//
// Muon PID
// I Mawby, D Brailsford & M Wallbank, May 2023
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
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
#include "larreco/RecoAlg/ShowerEnergyAlg.h"

// c++
#include <vector>

// ROOT
#include "TTree.h"
#include "TMVA/Reader.h"

namespace FDSelection {
  class PandizzleAlg;
}

class FDSelection::PandizzleAlg {

 public:
  enum Vars{
    kMichelNHits = 0,
    kMichelElectronMVA,
    kMichelRecoEnergyPlane2,
    kTrackDeflecAngleSD,
    kTrackLength,
    kEvalRatio,
    kConcentration,
    kCoreHaloRatio,
    kConicalness,
    kdEdxStart,
    kdEdxEnd,
    kdEdxEndRatio,
    kTerminatingValue //terminates the enum and not an actual variable
  };

  using InputVarsToReader = std::map<Vars, std::unique_ptr<Float_t>>;

  class Record {
    public:
      Record(const InputVarsToReader &inputVars, const Float_t mvaScore, const bool isFilled);

      Float_t GetVar(const FDSelection::PandizzleAlg::Vars var);
      bool IsFilled();
      Float_t GetMVAScore();

    private:
      using InputVars = std::map<Vars, Float_t>;
      InputVars fInputs;
      Float_t fMVAScore;
      bool fIsFilled;
  };

  PandizzleAlg(const fhicl::ParameterSet& pset);

  void Run(const art::Event& evt);
  Record RunPID(const art::Ptr<recob::Track> pTrack, const art::Event& evt);

 private:
  void InitialiseTrees();
  void BookTreeInt(TTree *tree, std::string branch_name);
  void BookTreeFloat(TTree *tree, std::string branch_name);
  void BookTreeBool(TTree *tree, std::string branch_name);
  void ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  int CountPFPWithPDG(const std::vector<art::Ptr<recob::PFParticle> > & pfps, int pdg);
  void FillTruthInfo(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  void FillMVAVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  void FillMichelElectronVariables(const art::Ptr<recob::PFParticle> mu_pfp, const art::Event& evt);
  void FillTrackVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  void CalculateTrackDeflection(const art::Ptr<recob::Track> track);
  void CalculateTrackLengthVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  void FillTree();
  void ResetTreeVariables();

  Float_t* GetVarPtr(const FDSelection::PandizzleAlg::Vars var);
  void SetVar(const FDSelection::PandizzleAlg::Vars var, const Float_t value);
  Record ReturnEmptyRecord();

  //Algs
  shower::ShowerEnergyAlg fShowerEnergyAlg;

  // module labels
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPIDModuleLabel;
  std::string fRecoModuleLabel;

  //Input params
  double fIPMichelCandidateDistance;

  // tree
  bool fMakeSelectionTrainingTrees;
  bool fReducedTreeMode;
  TTree* fSignalTrackTree;
  TTree *fBackgroundTrackTree;

  std::string fPandizzleWeightFileName;
  TMVA::Reader fPandizzleReader;
  InputVarsToReader fInputsToReader;

  struct VarHolder
  {
    std::map<std::string, int> IntVars;
    std::map<std::string, float> FloatVars;
    std::map<std::string, bool> BoolVars;
  };

  VarHolder fVarHolder;

  //Services
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<art::TFileService> tfs;
};

Float_t FDSelection::PandizzleAlg::Record::GetVar(const FDSelection::PandizzleAlg::Vars var)
{
  return (fInputs.at(var));
}

bool FDSelection::PandizzleAlg::Record::IsFilled()
{
  return fIsFilled;
}

Float_t FDSelection::PandizzleAlg::Record::GetMVAScore()
{
  return fMVAScore;
}

Float_t* FDSelection::PandizzleAlg::GetVarPtr(const FDSelection::PandizzleAlg::Vars var)
{
  return fInputsToReader.at(var).get();
}

void FDSelection::PandizzleAlg::SetVar(const FDSelection::PandizzleAlg::Vars var, const Float_t value)
{
  std::map<Vars, std::unique_ptr<Float_t> >::iterator itr(fInputsToReader.find(var));

  if (itr != fInputsToReader.end())
    *(itr->second.get()) = value;

  return;
}

#endif
