///////////////////////////////////////////////
// PandoraNuSelection produced module
//
// Add Pandora-based numu/nue selection scores to root file
// I Mawby        May 2023
///////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

// ART
#include "art/Utilities/make_tool.h" 

// LARSOFT
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "dunereco/FDSelections/FDSelectionData/PandoraNuSelection.h"
#include "dunereco/FDSelections/pandizzle/PandizzleAlg.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"
#include "tools/RecoShowerSelector.h"
#include "tools/RecoTrackSelector.h"

constexpr int kDefInt = -9999;
constexpr int kDefDoub = (double)(kDefInt);

class PandoraNuSelection;

class PandoraNuSelection : public art::EDProducer {
public:
  explicit PandoraNuSelection(fhicl::ParameterSet const& pset);
  PandoraNuSelection(PandoraNuSelection const&) = delete;
  PandoraNuSelection(PandoraNuSelection&&) = delete;
  PandoraNuSelection& operator=(PandoraNuSelection const&) = delete;
  PandoraNuSelection& operator=(PandoraNuSelection&&) = delete;
  void produce(art::Event& evt) override;
  void beginJob() override;
  void endJob() override;

  bool GetSelectedTrack(art::Event& evt, art::Ptr<recob::Track>& selTrack);
  double GetPandizzleScore(art::Event& evt, art::Ptr<recob::Track> selTrack);
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt);
  bool GetSelectedShower(art::Event& evt, art::Ptr<recob::Shower>& selShower);
  void SetPandrizzleScores(art::Event& evt, art::Ptr<recob::Shower> selShower, std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection);
  void FillVertexInformation(art::Event const & evt);

private:
  // Tools
  std::unique_ptr<FDSelectionTools::RecoShowerSelector> fRecoShowerSelector;
  std::unique_ptr<FDSelectionTools::RecoTrackSelector> fRecoTrackSelector;

  // Algs
  FDSelection::PandizzleAlg fPandizzleAlg;
  FDSelection::PandrizzleAlg fPandrizzleAlg;
};

///////////////////////////////////////////////////////////////////////////////////////////////////

PandoraNuSelection::PandoraNuSelection(fhicl::ParameterSet const& pset) : 
  EDProducer{pset},
  fRecoShowerSelector{art::make_tool<FDSelectionTools::RecoShowerSelector>(pset.get<fhicl::ParameterSet>("RecoShowerSelectorTool"))},
  fRecoTrackSelector{art::make_tool<FDSelectionTools::RecoTrackSelector>(pset.get<fhicl::ParameterSet>("RecoTrackSelectorTool"))},
  fPandizzleAlg(pset),
  fPandrizzleAlg(pset)
{
  produces<pandoranusel::PandoraNuSelection>();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandoraNuSelection::produce(art::Event& evt)
{
  std::unique_ptr<pandoranusel::PandoraNuSelection> pandoraNuSelection = std::make_unique<pandoranusel::PandoraNuSelection>();

  art::Ptr<recob::Track> selTrack;

  if (GetSelectedTrack(evt, selTrack))
    pandoraNuSelection->selTrackPandizzleScore = GetPandizzleScore(evt, selTrack);
  else
    pandoraNuSelection->selTrackPandizzleScore = kDefDoub;

  art::Ptr<recob::Shower> selShower;

  if (GetSelectedShower(evt, selShower))
  {
    SetPandrizzleScores(evt, selShower, pandoraNuSelection);
  }
  else
  {
    pandoraNuSelection->selShowerBackupPandrizzleScore = kDefDoub;
    pandoraNuSelection->selShowerEnhancedPandrizzleScore = kDefDoub;
  }

  evt.put(std::move(pandoraNuSelection));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

bool PandoraNuSelection::GetSelectedTrack(art::Event& evt, art::Ptr<recob::Track>& selTrack)
{
  selTrack = fRecoTrackSelector->FindSelectedTrack(evt);

  if (!selTrack.isAvailable()) 
  {
    std::cout << "FDSelection::PandoraNuSelection - no track returned from selection" << std::endl; 
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

double PandoraNuSelection::GetPandizzleScore(art::Event& evt, art::Ptr<recob::Track> selTrack)
{
  FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(selTrack, evt));
  return pandizzleRecord.GetMVAScore();
}

///////////////////////////////////////////////////////////////////////////////////////////////////

bool PandoraNuSelection::GetSelectedShower(art::Event& evt, art::Ptr<recob::Shower>& selShower)
{
  selShower = fRecoShowerSelector->FindSelectedShower(evt);

  if (!selShower.isAvailable()) 
  {
    std::cout << "FDSelection::PandoraNuSelection::RunShowerSelection - no shower returned from selection" << std::endl;
    return false;
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandoraNuSelection::SetPandrizzleScores(art::Event& evt, art::Ptr<recob::Shower> selShower, std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection)
{
  FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(selShower, evt));
  float pandrizzleBDTMethod(pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod));
  float pandrizzleScore(pandrizzleRecord.GetMVAScore());

  pandoraNuSelection->selShowerBackupPandrizzleScore = (std::fabs(pandrizzleBDTMethod - 1.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;
  pandoraNuSelection->selShowerEnhancedPandrizzleScore = (std::fabs(pandrizzleBDTMethod - 2.0) < std::numeric_limits<float>::epsilon()) ? pandrizzleScore : -9999.f;

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandoraNuSelection::beginJob()
{
}

void PandoraNuSelection::endJob()
{
}

DEFINE_ART_MODULE(PandoraNuSelection)
