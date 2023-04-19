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

namespace FDSelection {
  class PandizzleAlg;
}

class FDSelection::PandizzleAlg {

 public:
  PandizzleAlg(const fhicl::ParameterSet& pset);

  void Run(const art::Event& evt);
  void ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);

  int GetIntVar(std::string name) { return fVarHolder.IntVars[name]; };
  float GetFloatVar(std::string name) { return fVarHolder.FloatVars[name]; };
  bool GetBoolVar(std::string name) { return fVarHolder.BoolVars[name]; };

 private:
  void InitialiseTrees();
  void BookTreeInt(TTree *tree, std::string branch_name);
  void BookTreeFloat(TTree *tree, std::string branch_name);
  void BookTreeBool(TTree *tree, std::string branch_name);
  void FillTruthInfo(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  int CountPFPWithPDG(const std::vector<art::Ptr<recob::PFParticle> > & pfps, int pdg);
  void FillMichelElectronVariables(const art::Ptr<recob::PFParticle> mu_pfp, const art::Event& evt);
  void FillTrackVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  void CalculateTrackDeflection(const art::Ptr<recob::Track> track);
  void CalculateTrackLengthVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);
  void FillTree();
  void ResetTreeVariables();

  // module labels
  std::string fTrackModuleLabel;
  std::string fShowerModuleLabel;
  std::string fPIDModuleLabel;
  std::string fPFParticleModuleLabel;
  std::string fSpacePointModuleLabel;
  std::string fClusterModuleLabel;

  // tree
  bool fMakeTree;
  TTree* fSignalTrackTree;
  TTree *fBackgroundTrackTree;

  //Input params
  double fIPMichelCandidateDistance;

  //Algs
  shower::ShowerEnergyAlg fShowerEnergyAlg;

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

#endif
