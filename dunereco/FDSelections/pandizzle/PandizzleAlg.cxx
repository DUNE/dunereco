///////////////////////////////////////////////
// PandizzleAlg.cxx
//
// Muon PID
// I Mawby, D Brailsford & M Wallbank, May 2023
///////////////////////////////////////////////

#include "PandizzleAlg.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"

namespace
{
  constexpr Float_t kDefValue(std::numeric_limits<Float_t>::lowest());

  using namespace FDSelection;

  void Reset(PandizzleAlg::InputVarsToReader &inputVarsToReader)
  {
    for (PandizzleAlg::Vars var = PandizzleAlg::kMichelNHits; var < PandizzleAlg::kTerminatingValue; var=static_cast<PandizzleAlg::Vars>(static_cast<int>(var)+1))
    {
      auto [itr, inserted] = inputVarsToReader.try_emplace(var, std::make_unique<Float_t>(kDefValue));
      if (!inserted)
        *(itr->second.get()) = kDefValue;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandizzleAlg::Record::Record(const InputVarsToReader &inputVarsToReader, const Float_t mvaScore, const bool isFilled) :
  fMVAScore(mvaScore),
  fIsFilled(isFilled)
{
  for (PandizzleAlg::Vars var = PandizzleAlg::kMichelNHits; var < PandizzleAlg::kTerminatingValue; var=static_cast<PandizzleAlg::Vars>(static_cast<int>(var)+1)) 
    fInputs.try_emplace(var, *(inputVarsToReader.at(var)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandizzleAlg::PandizzleAlg(const fhicl::ParameterSet& pset) :
  fShowerEnergyAlg(pset.get<fhicl::ParameterSet>("ShowerEnergyAlg")),
  fTrackModuleLabel(pset.get<std::string>("ModuleLabels.TrackModuleLabel")),
  fShowerModuleLabel(pset.get<std::string>("ModuleLabels.ShowerModuleLabel")),
  fPIDModuleLabel(pset.get<std::string>("ModuleLabels.PIDModuleLabel")),
  fPFParticleModuleLabel(pset.get<std::string>("ModuleLabels.PFParticleModuleLabel")),
  fSpacePointModuleLabel( pset.get<std::string>("ModuleLabels.SpacePointModuleLabel")),
  fClusterModuleLabel(pset.get<std::string>("ModuleLabels.ClusterModuleLabel")),
  fIPMichelCandidateDistance(pset.get<double>("MichelCandidateDistance")),
  fMakeSelectionTrainingTrees(pset.get<bool>("MakeSelectionTrainingTrees")),
  fPandizzleWeightFileName(pset.get< std::string > ("PandizzleWeightFileName")),
  fPandizzleReader("", 0)
{
  Reset(fInputsToReader);

  fPandizzleReader.AddVariable("PFPMichelNHits", GetVarPtr(kMichelNHits));
  fPandizzleReader.AddVariable("PFPMichelElectronMVA", GetVarPtr(kMichelElectronMVA)); 
  fPandizzleReader.AddVariable("PFPMichelRecoEnergyPlane2", GetVarPtr(kMichelRecoEnergyPlane2));
  fPandizzleReader.AddVariable("PFPTrackDeflecAngleSD", GetVarPtr(kTrackDeflecAngleSD));
  fPandizzleReader.AddVariable("PFPTrackLength", GetVarPtr(kTrackLength));
  fPandizzleReader.AddVariable("PFPTrackEvalRatio", GetVarPtr(kEvalRatio));
  fPandizzleReader.AddVariable("PFPTrackConcentration", GetVarPtr(kConcentration));
  fPandizzleReader.AddVariable("PFPTrackCoreHaloRatio", GetVarPtr(kCoreHaloRatio));
  fPandizzleReader.AddVariable("PFPTrackConicalness", GetVarPtr(kConicalness));
  fPandizzleReader.AddVariable("PFPTrackdEdxStart", GetVarPtr(kdEdxStart));
  fPandizzleReader.AddVariable("PFPTrackdEdxEnd", GetVarPtr(kdEdxEnd));
  fPandizzleReader.AddVariable("PFPTrackdEdxEndRatio", GetVarPtr(kdEdxEndRatio));

  const std::string weightFileName(fPandizzleWeightFileName);
  std::cout << "fPandizzleWeightFileName: " << fPandizzleWeightFileName << std::endl;
  std::string weightFilePath;
  cet::search_path sP("FW_SEARCH_PATH");
  std::cout << "sP: " << sP << std::endl;
  sP.find_file(weightFileName, weightFilePath);
  fPandizzleReader.BookMVA("BDTG", weightFilePath);

  if (fMakeSelectionTrainingTrees)
  {
      InitialiseTrees();
      ResetTreeVariables();
  }
}

////////////////////////

void FDSelection::PandizzleAlg::InitialiseTrees() {
  fSignalTrackTree = tfs->make<TTree>("DizzleSigTrackTree", "Pandizzle Signal Track Tree");
  fBackgroundTrackTree = tfs->make<TTree>("DizzleBgTrackTree", "Pandizzle Background Track Tree");

  for (TTree* tree : {fSignalTrackTree, fBackgroundTrackTree})
  {
    BookTreeInt(tree, "Event");
    BookTreeInt(tree, "Run");
    BookTreeInt(tree, "SubRun");
    BookTreeInt(tree, "PFPPDG");
    BookTreeInt(tree, "PFPNHits");
    BookTreeInt(tree, "PFPTruePDG");
    BookTreeInt(tree, "PFPTrueMotherID");
    BookTreeFloat(tree, "PFPTrueMomT");
    BookTreeFloat(tree, "PFPTrueMomX");
    BookTreeFloat(tree, "PFPTrueMomY");
    BookTreeFloat(tree, "PFPTrueMomZ");
    BookTreeBool(tree,"PFPTrueIsMichelDecay");
    BookTreeInt(tree, "PFPNChildren");
    BookTreeInt(tree, "PFPNShowerChildren");
    BookTreeInt(tree, "PFPNTrackChildren");
    BookTreeFloat(tree, "PFPMichelDist");
    BookTreeInt(tree, "PFPMichelNHits");
    BookTreeInt(tree, "PFPMichelTrueID");
    BookTreeInt(tree, "PFPMichelTrueMotherID");
    BookTreeInt(tree, "PFPMichelTruePDG");
    BookTreeFloat(tree, "PFPMichelElectronMVA");
    BookTreeFloat(tree, "PFPMichelRecoEnergyPlane0");
    BookTreeFloat(tree, "PFPMichelRecoEnergyPlane1");
    BookTreeFloat(tree, "PFPMichelRecoEnergyPlane2");
    BookTreeFloat(tree,"PFPTrackDeflecAngleMean");
    BookTreeFloat(tree,"PFPTrackDeflecAngleVar");
    BookTreeFloat(tree,"PFPTrackDeflecAngleSD");
    BookTreeInt(tree,"PFPTrackDeflecNAngles");
    BookTreeBool(tree, "MVAVarsFilled");
    BookTreeFloat(tree,"PFPTrackEvalRatio");
    BookTreeFloat(tree,"PFPTrackConcentration");
    BookTreeFloat(tree,"PFPTrackCoreHaloRatio");
    BookTreeFloat(tree,"PFPTrackConicalness");
    BookTreeFloat(tree,"PFPTrackdEdxStart");
    BookTreeFloat(tree,"PFPTrackdEdxEnd");
    BookTreeFloat(tree,"PFPTrackdEdxEndRatio");
    BookTreeFloat(tree,"PFPTrackMuonMVA");
    BookTreeFloat(tree,"PFPTrackProtonMVA");
    BookTreeFloat(tree,"PFPTrackPionMVA");
    BookTreeFloat(tree,"PFPTrackPhotonMVA");
    BookTreeFloat(tree,"PFPTrackElectronMVA");
    BookTreeFloat(tree,"PFPTrackLengthRatio");
    BookTreeFloat(tree,"PFPTrackLength");
    BookTreeFloat(tree,"PFPTrackOtherLengths");
    BookTreeInt(tree,"PFPTrackIsLongest");
  }
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::BookTreeInt(TTree *tree, std::string branch_name){
  fVarHolder.IntVars[branch_name] = -9999;
  tree->Branch(branch_name.c_str(), &(fVarHolder.IntVars[branch_name]));
  return;
}

void FDSelection::PandizzleAlg::BookTreeFloat(TTree *tree, std::string branch_name){
  fVarHolder.FloatVars[branch_name] = -9999.0;
  tree->Branch(branch_name.c_str(), &(fVarHolder.FloatVars[branch_name]));
  return;
}

void FDSelection::PandizzleAlg::BookTreeBool(TTree *tree, std::string branch_name){
  fVarHolder.BoolVars[branch_name] = 0;
  tree->Branch(branch_name.c_str(), &(fVarHolder.BoolVars[branch_name]));
  return;
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::Run(const art::Event& evt) 
{
  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return;
  
  art::Ptr<recob::PFParticle> neutrinoPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(neutrinoPFP, evt, fPFParticleModuleLabel);

  for (art::Ptr<recob::PFParticle> childPFP : childPFPs)
  {
    if (childPFP->PdgCode() != 13)
      continue;

    ProcessPFParticle(childPFP, evt);
    FillTree();
    ResetTreeVariables();
  }
  
  return;
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::ProcessPFParticle(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
  // Get event, subrun, run
  fVarHolder.IntVars["Run"] = evt.run();
  fVarHolder.IntVars["SubRun"] = evt.subRun();
  fVarHolder.IntVars["Event"] = evt.event();

  // PDG
  fVarHolder.IntVars["PFPPDG"] = pfp->PdgCode();

  // Get the PFP hits
  std::vector<art::Ptr<recob::Hit>> pfp_hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fClusterModuleLabel);
  fVarHolder.IntVars["PFPNHits"] = pfp_hits.size();

  //Count all of the children
  std::vector<art::Ptr<recob::PFParticle>> child_pfps = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(pfp, evt, fPFParticleModuleLabel);
  fVarHolder.IntVars["PFPNChildren"] = pfp->NumDaughters();
  fVarHolder.IntVars["PFPNShowerChildren"] = CountPFPWithPDG(child_pfps, 11);
  fVarHolder.IntVars["PFPNTrackChildren"] = CountPFPWithPDG(child_pfps, 13);
  
  FillTruthInfo(pfp, evt);

  FillMVAVariables(pfp, evt);

  if (child_pfps.size() != 0)
    FillMichelElectronVariables(pfp, evt);

  FillTrackVariables(pfp, evt);

  return;
}


///////////////////////////////////////////////////////

int FDSelection::PandizzleAlg::CountPFPWithPDG(const std::vector<art::Ptr<recob::PFParticle> > & pfps, int pdg)
{
  int NPFP = 0;

  for (art::Ptr<recob::PFParticle> pfp : pfps)
  {
    int pfp_pdg = pfp->PdgCode();
    if (pfp_pdg == pdg) NPFP++;
  }

  return NPFP;
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::FillTruthInfo(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
  std::vector<art::Ptr<recob::Hit> > pfp_hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(pfp, evt, fClusterModuleLabel);

  simb::MCParticle* matchedMCParticle = nullptr;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, pfp_hits, 1); 
  if (g4id > 0)
    matchedMCParticle = pi_serv->ParticleList().at(g4id);

  if (matchedMCParticle)
  {
    fVarHolder.IntVars["PFPTruePDG"] = matchedMCParticle->PdgCode();
    fVarHolder.IntVars["PFPTrueMotherID"] = matchedMCParticle->Mother();
    fVarHolder.FloatVars["PFPTrueMomX"] = matchedMCParticle->Momentum(0).X();
    fVarHolder.FloatVars["PFPTrueMomY"] = matchedMCParticle->Momentum(0).Y();
    fVarHolder.FloatVars["PFPTrueMomZ"] = matchedMCParticle->Momentum(0).Z();
    fVarHolder.FloatVars["PFPTrueMomT"] = matchedMCParticle->Momentum(0).T();

    if (std::abs(fVarHolder.IntVars["PFPTruePDG"]) == 13)
    {
      bool has_electron = 0;
      bool has_numu = 0;
      bool has_nue = 0;

      for (int i_child = 0; i_child < matchedMCParticle->NumberDaughters(); i_child++)
      {
        const simb::MCParticle* child_particle = pi_serv->ParticleList().at(matchedMCParticle->Daughter(i_child));
        int abs_pdg = std::abs(child_particle->PdgCode());
        if (abs_pdg == 11) has_electron = 1;
        else if (abs_pdg == 12) has_nue = 1;
        else if (abs_pdg == 14) has_numu = 1;
      }

      if (has_electron && has_nue && has_numu) 
        fVarHolder.BoolVars["PFPTrueIsMichelDecay"] = 1;
    }
  }
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::FillMVAVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{

  if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
    return;

  art::Ptr<recob::Track> pfp_track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);
  art::FindOneP<anab::MVAPIDResult> findPIDResult(std::vector<art::Ptr<recob::Track> >{pfp_track}, evt, fPIDModuleLabel);
  art::Ptr<anab::MVAPIDResult> pid(findPIDResult.at(0));

  if (pid.isAvailable())
  {
    fVarHolder.FloatVars["PFPTrackEvalRatio"] = pid->evalRatio;
    fVarHolder.FloatVars["PFPTrackConcentration"] = pid->concentration;
    fVarHolder.FloatVars["PFPTrackCoreHaloRatio"] = pid->coreHaloRatio;
    fVarHolder.FloatVars["PFPTrackConicalness"] = pid->conicalness;
    fVarHolder.FloatVars["PFPTrackdEdxStart"] = pid->dEdxStart;
    fVarHolder.FloatVars["PFPTrackdEdxEnd"] = pid->dEdxEnd;
    fVarHolder.FloatVars["PFPTrackdEdxEndRatio"] = pid->dEdxEndRatio;
    std::map<std::string,double> mvaOutMap = pid->mvaOutput;  
    fVarHolder.FloatVars["PFPTrackMuonMVA"] = mvaOutMap["muon"];
    fVarHolder.FloatVars["PFPTrackProtonMVA"] = mvaOutMap["proton"];
    fVarHolder.FloatVars["PFPTrackPionMVA"] = mvaOutMap["pion"];
    fVarHolder.FloatVars["PFPTrackPhotonMVA"] = mvaOutMap["photon"];
    fVarHolder.FloatVars["PFPTrackElectronMVA"] = mvaOutMap["electron"];
    fVarHolder.BoolVars["MVAVarsFilled"] = true;
  }
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::FillMichelElectronVariables(const art::Ptr<recob::PFParticle> mu_pfp, const art::Event& evt)
{
  //Assign the branch value to a new default value to indicate that we have a track, but don't necessarily have a Michel candidate (rather than a global default value)
  fVarHolder.FloatVars["PFPMichelDist"] = -2.;
  fVarHolder.IntVars["PFPMichelNHits"] = -2;
  fVarHolder.FloatVars["PFPMichelElectronMVA"] = -2;
  fVarHolder.FloatVars["PFPMichelRecoEnergyPlane0"] = -2;
  fVarHolder.FloatVars["PFPMichelRecoEnergyPlane1"] = -2;
  fVarHolder.FloatVars["PFPMichelRecoEnergyPlane2"] = -2;

  std::vector<art::Ptr<recob::PFParticle>> child_pfps = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(mu_pfp, evt, fPFParticleModuleLabel);
  
  //Get the shower handle
  art::Handle< std::vector<recob::Shower> > showerListHandle;
  if (!(evt.getByLabel(fShowerModuleLabel, showerListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::Shower> with module label: " << fShowerModuleLabel;
    return;
  }

  //Get the track handle
  art::Handle< std::vector<recob::Track> > trackListHandle;
  if (!(evt.getByLabel(fTrackModuleLabel, trackListHandle))){
    mf::LogWarning("PandizzleAlg") << "Unable to find std::vector<recob::Track> with module label: " << fTrackModuleLabel;
    return;
  }

  art::FindManyP<anab::MVAPIDResult> fmpidt(trackListHandle, evt, fPIDModuleLabel);
  art::FindManyP<anab::MVAPIDResult> fmpids(showerListHandle, evt, fPIDModuleLabel);

  // Find closest particle to end of track
  art::Ptr<recob::Track> mu_track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(mu_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);
  TVector3 mu_end_position(mu_track->End().X(), mu_track->End().Y(), mu_track->End().Z());;
  art::Ptr<recob::PFParticle> closest_child_pfp;
  art::Ptr<anab::MVAPIDResult> closest_child_pfp_mva_pid_result;
  double closest_distance = std::numeric_limits<double>::max();

  for (art::Ptr<recob::PFParticle> child_pfp : child_pfps)
  {
    art::Ptr<anab::MVAPIDResult> child_pfp_mva_pid_result;
    TVector3 child_start_pos;

    if (dune_ana::DUNEAnaPFParticleUtils::IsTrack(child_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
    {
      art::Ptr<recob::Track> child_track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(child_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);
      child_start_pos.SetXYZ(child_track->Start().X(), child_track->Start().Y(), child_track->Start().Z());
      std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(child_track.key());
      child_pfp_mva_pid_result = pids.at(0);
    }
    else if (dune_ana::DUNEAnaPFParticleUtils::IsShower(child_pfp, evt, fPFParticleModuleLabel, fShowerModuleLabel))
    {
      art::Ptr<recob::Shower> child_shower = dune_ana::DUNEAnaPFParticleUtils::GetShower(child_pfp, evt, fPFParticleModuleLabel, fShowerModuleLabel);
      child_start_pos.SetXYZ(child_shower->ShowerStart().X(), child_shower->ShowerStart().Y(), child_shower->ShowerStart().Z());
      std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpids.at(child_shower.key());
      child_pfp_mva_pid_result = pids.at(0);
    }
    else {
      continue;
    }

    const double curr_distance = (child_start_pos - mu_end_position).Mag();
    if (curr_distance < closest_distance)
    {
      closest_distance = curr_distance;
      closest_child_pfp = child_pfp;
      closest_child_pfp_mva_pid_result = child_pfp_mva_pid_result;
    }
  }

  //Check that the PFP passes our Michel requirement
  if (closest_distance > fIPMichelCandidateDistance) 
    return;

  fVarHolder.FloatVars["PFPMichelDist"] = closest_distance;

  std::vector<art::Ptr<recob::Hit> > michel_hits = dune_ana::DUNEAnaPFParticleUtils::GetHits(closest_child_pfp, evt, fClusterModuleLabel);
  fVarHolder.IntVars["PFPMichelNHits"] = std::min((int)(michel_hits.size()), 100);

  simb::MCParticle* matchedMCParticle = nullptr;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, michel_hits, 1); 
  if (g4id > 0)
    matchedMCParticle = pi_serv->ParticleList().at(g4id);

  if (matchedMCParticle)
  {
    fVarHolder.IntVars["PFPMichelTrueID"] = matchedMCParticle->TrackId();
    fVarHolder.IntVars["PFPMichelTrueMotherID"] = matchedMCParticle->Mother();
    fVarHolder.IntVars["PFPMichelTruePDG"] = matchedMCParticle->PdgCode();
  }

  std::map<std::string,double> mvaOutMap = closest_child_pfp_mva_pid_result->mvaOutput;
  fVarHolder.FloatVars["PFPMichelElectronMVA"] = mvaOutMap["electron"];

  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob(clockData);

  fVarHolder.FloatVars["PFPMichelRecoEnergyPlane0"] = fShowerEnergyAlg.ShowerEnergy(clockData, detProp, michel_hits, 0);
  fVarHolder.FloatVars["PFPMichelRecoEnergyPlane1"] = fShowerEnergyAlg.ShowerEnergy(clockData, detProp, michel_hits, 1);
  fVarHolder.FloatVars["PFPMichelRecoEnergyPlane2"] = fShowerEnergyAlg.ShowerEnergy(clockData, detProp, michel_hits, 2);

  return;
}

///////////////////////////////////////////////////////

void FDSelection::PandizzleAlg::FillTrackVariables(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt)
{
  if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
    return;

  art::Ptr<recob::Track> pfp_track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);

  CalculateTrackDeflection(pfp_track);
  CalculateTrackLengthVariables(pfp, evt);

  return;
}

/////////////////////////////////////////////////////// 

void FDSelection::PandizzleAlg::CalculateTrackDeflection(const art::Ptr<recob::Track> track){
  // This follows the MicroBooNE method for MCS

  //Get the number of points
  size_t NPoints = track->NumberTrajectoryPoints();
  //Store the directions between adjacent points on a vector
  std::vector<TVector3> directions;
  for (size_t i_point = 0; i_point < NPoints-1; i_point++){
    TVector3 position_i(track->TrajectoryPoint(i_point).position.X(), track->TrajectoryPoint(i_point).position.Y(), track->TrajectoryPoint(i_point).position.Z());
    TVector3 position_iplus1(track->TrajectoryPoint(i_point+1).position.X(), track->TrajectoryPoint(i_point+1).position.Y(), track->TrajectoryPoint(i_point+1).position.Z());
    TVector3 direction = (position_iplus1-position_i).Unit();
    directions.push_back(direction);
  }

  //Loop through the direction and compare adjacent elements
  std::vector<double> deflection_angles;
  for (size_t i_dir = 0; i_dir < directions.size()-1; i_dir++){
    //Aim: rotate both direction so that the first direction is parallel to the z-axis.  
    //Then take the x-projection of scattered track and calculate the angle between that and the first direction
    TVector3 z_axis(0,0,1);
    TVector3 direction_first = directions[i_dir];
    TVector3 direction_second = directions[i_dir+1];

    //Ignore if either direction is 0 (i.e. not 1)
    if (direction_first.Mag() < 0.999 || direction_second.Mag() < 0.999){
      continue;
    }
    double angle_dir_first_z_axis = direction_first.Angle(z_axis);
    TVector3 orthogonal_vector = direction_first.Cross(z_axis);
    if (orthogonal_vector.Unit().Mag() < 0.999){
        continue;
    }
    direction_first.Rotate(angle_dir_first_z_axis, orthogonal_vector);
    direction_second.Rotate(angle_dir_first_z_axis, orthogonal_vector);

    //Now work out the angle between the vectors in the x-z plane
    direction_first.SetY(0);
    direction_second.SetY(0);
    double dot_product = direction_first.Dot(direction_second);
    dot_product = std::min(std::max(dot_product,-1.),1.);
    double angle = acos(dot_product) * 180/3.142;

    //define +x as a +angle
    if (direction_second.X() < 0) angle*=-1;
    deflection_angles.push_back(angle);
  }

  double angle_mean = 0;
  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++){
    angle_mean += deflection_angles[i_angle];
  }

  if (deflection_angles.size()>0) angle_mean/=deflection_angles.size();
  else angle_mean=-100;

  double angle_var = 0;
  for (size_t i_angle = 0; i_angle < deflection_angles.size(); i_angle++){
    angle_var = (deflection_angles[i_angle] - angle_mean)*(deflection_angles[i_angle] - angle_mean);
  }

  if (deflection_angles.size() > 1) angle_var /= (deflection_angles.size()-1);
  else angle_var = -2.;

  fVarHolder.FloatVars["PFPTrackDeflecAngleMean"] = angle_mean;
  fVarHolder.FloatVars["PFPTrackDeflecAngleVar"] = angle_var;
  fVarHolder.FloatVars["PFPTrackDeflecAngleSD"] = (angle_var > 0.0) ? sqrt(angle_var) : -2.f;
  fVarHolder.IntVars["PFPTrackDeflecNAngles"] = deflection_angles.size();
  return;
}

/////////////////////////////////////////////////////// 

void FDSelection::PandizzleAlg::CalculateTrackLengthVariables(const art::Ptr<recob::PFParticle> track_pfp, const art::Event& evt)
{
  if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(track_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
    return;

  art::Ptr<recob::Track> track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(track_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);

  double average_length = 0;
  int ntracks = 0;
  double candidate_track_length = track->Length(); 
  bool is_longest_track = 1;

  art::Ptr<recob::PFParticle> neutrinoPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> childPFPs = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(neutrinoPFP, evt, fPFParticleModuleLabel);

  //Now grab the primary children of these PFP
  for (art::Ptr<recob::PFParticle> child_pfp : childPFPs)
  {
    //Don't assess the particle we care about in the average
    if (child_pfp->Self() == track_pfp->Self()) continue;

    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(child_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel))
      continue;

    art::Ptr<recob::Track> matched_track = dune_ana::DUNEAnaPFParticleUtils::GetTrack(child_pfp, evt, fPFParticleModuleLabel, fTrackModuleLabel);
    const double matched_track_length = matched_track->Length();

    if (matched_track_length > candidate_track_length)
      is_longest_track = 0;

    average_length += matched_track->Length();
    ntracks++;
  }

  if (ntracks > 0) average_length/=ntracks;
  else average_length = 0;
  fVarHolder.FloatVars["PFPTrackOtherLengths"] = average_length;
  fVarHolder.FloatVars["PFPTrackLength"] = candidate_track_length;
  fVarHolder.FloatVars["PFPTrackLengthRatio"] = average_length/candidate_track_length;
  fVarHolder.IntVars["PFPTrackIsLongest"] = is_longest_track;
}

/////////////////////////

void FDSelection::PandizzleAlg::FillTree()
{
  if (fVarHolder.FloatVars["PFPTrackLength"] < 0)
  {
    mf::LogWarning("PandizzleAlg") << "Found candidate track with length < 0.  Not filling tree";
    return;
  }

  if (std::abs(fVarHolder.IntVars["PFPTruePDG"]) == 13)
  {
      fSignalTrackTree->Fill();
  }
  else 
  {
      fBackgroundTrackTree->Fill();
  }

  return;
}

/////////////////////////////////////////////////////// 

void FDSelection::PandizzleAlg::ResetTreeVariables()
{
  for (std::map<std::string, int>::iterator mapIt = fVarHolder.IntVars.begin(); mapIt != fVarHolder.IntVars.end(); mapIt++){
    mapIt->second = -9999;
  }
  for (std::map<std::string, float>::iterator mapIt = fVarHolder.FloatVars.begin(); mapIt != fVarHolder.FloatVars.end(); mapIt++){
    mapIt->second = -9999.;
  }
  for (std::map<std::string, bool>::iterator mapIt = fVarHolder.BoolVars.begin(); mapIt != fVarHolder.BoolVars.end(); mapIt++){
    mapIt->second = 0;
  }

  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

FDSelection::PandizzleAlg::Record FDSelection::PandizzleAlg::RunPID(const art::Ptr<recob::Track> pTrack, const art::Event& evt) 
{
  art::Ptr<recob::PFParticle> pfp = dune_ana::DUNEAnaTrackUtils::GetPFParticle(pTrack, evt, fTrackModuleLabel);

  fVarHolder.BoolVars["MVAVarsFilled"] = false;

  ResetTreeVariables();
  ProcessPFParticle(pfp, evt);

  if (fVarHolder.BoolVars["MVAVarsFilled"])
  {
    SetVar(kMichelNHits, (float)fVarHolder.IntVars["PFPMichelNHits"]);
    SetVar(kMichelElectronMVA, fVarHolder.FloatVars["PFPMichelElectronMVA"]);
    SetVar(kMichelRecoEnergyPlane2, fVarHolder.FloatVars["PFPMichelRecoEnergyPlane2"]);
    SetVar(kTrackDeflecAngleSD, fVarHolder.FloatVars["PFPTrackDeflecAngleSD"]);
    SetVar(kTrackLength, fVarHolder.FloatVars["PFPTrackLength"]);
    SetVar(kEvalRatio, fVarHolder.FloatVars["PFPTrackEvalRatio"]);
    SetVar(kConcentration, fVarHolder.FloatVars["PFPTrackConcentration"]);
    SetVar(kCoreHaloRatio, fVarHolder.FloatVars["PFPTrackCoreHaloRatio"]);
    SetVar(kConicalness, fVarHolder.FloatVars["PFPTrackConicalness"]);
    SetVar(kdEdxStart, fVarHolder.FloatVars["PFPTrackdEdxStart"]);
    SetVar(kdEdxEnd, fVarHolder.FloatVars["PFPTrackdEdxEnd"]);
    SetVar(kdEdxEndRatio, fVarHolder.FloatVars["PFPTrackdEdxEndRatio"]);

    return Record(fInputsToReader, fPandizzleReader.EvaluateMVA("BDTG"), true);
  }

  return ReturnEmptyRecord();
}

////////////////////////

FDSelection::PandizzleAlg::Record FDSelection::PandizzleAlg::ReturnEmptyRecord()
{
    Reset(fInputsToReader);
    return Record(fInputsToReader, kDefValue, false);
}

////////////////////////////////////////////////////////////////////////////////////////////////

