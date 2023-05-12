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
  void SetPandizzleScores(art::Event& evt, art::Ptr<recob::Track> selTrack, std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection);
  void SetNullPandizzleScores(std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection);
  art::Ptr<recob::PFParticle> GetPFParticleMatchedToTrack(art::Ptr<recob::Track> const track, art::Event const & evt);
  bool GetSelectedShower(art::Event& evt, art::Ptr<recob::Shower>& selShower);
  void SetPandrizzleScores(art::Event& evt, art::Ptr<recob::Shower> selShower, std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection);
  void SetNullPandrizzleScores(std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection);
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
      SetPandizzleScores(evt, selTrack, pandoraNuSelection);
  else
      SetNullPandizzleScores(pandoraNuSelection);

  art::Ptr<recob::Shower> selShower;

  if (GetSelectedShower(evt, selShower))
    SetPandrizzleScores(evt, selShower, pandoraNuSelection);
  else
    SetNullPandrizzleScores(pandoraNuSelection);

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

void PandoraNuSelection::SetPandizzleScores(art::Event& evt, art::Ptr<recob::Track> selTrack, std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection)
{
  FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(selTrack, evt));

  pandoraNuSelection->selTrackPandizzleScore = pandizzleRecord.GetMVAScore();
  pandoraNuSelection->selTrackMichelNHits = (double)pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelNHits);
  pandoraNuSelection->selTrackMichelElectronMVA = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelElectronMVA);
  pandoraNuSelection->selTrackMichelRecoEnergyPlane2 = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kMichelRecoEnergyPlane2);
  pandoraNuSelection->selTrackDeflecAngleSD = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackDeflecAngleSD);
  pandoraNuSelection->selTrackLength = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kTrackLength);
  pandoraNuSelection->selTrackEvalRatio = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kEvalRatio);
  pandoraNuSelection->selTrackConcentration = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConcentration);
  pandoraNuSelection->selTrackCoreHaloRatio = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kCoreHaloRatio);
  pandoraNuSelection->selTrackConicalness = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kConicalness);
  pandoraNuSelection->selTrackdEdxStart = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxStart);
  pandoraNuSelection->selTrackdEdxEnd = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEnd);
  pandoraNuSelection->selTrackdEdxEndRatio = pandizzleRecord.GetVar(FDSelection::PandizzleAlg::kdEdxEndRatio);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandoraNuSelection::SetNullPandizzleScores(std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection)
{
  pandoraNuSelection->selTrackPandizzleScore = kDefDoub;
  pandoraNuSelection->selTrackMichelNHits = kDefDoub;
  pandoraNuSelection->selTrackMichelElectronMVA = kDefDoub;
  pandoraNuSelection->selTrackMichelRecoEnergyPlane2 = kDefDoub;
  pandoraNuSelection->selTrackDeflecAngleSD = kDefDoub;
  pandoraNuSelection->selTrackLength = kDefDoub;
  pandoraNuSelection->selTrackEvalRatio = kDefDoub;
  pandoraNuSelection->selTrackConcentration = kDefDoub;
  pandoraNuSelection->selTrackCoreHaloRatio = kDefDoub;
  pandoraNuSelection->selTrackConicalness = kDefDoub;
  pandoraNuSelection->selTrackdEdxStart = kDefDoub;
  pandoraNuSelection->selTrackdEdxEnd = kDefDoub;
  pandoraNuSelection->selTrackdEdxEndRatio = kDefDoub;

  return;
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

  pandoraNuSelection->selShowerEvalRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEvalRatio);
  pandoraNuSelection->selShowerConcentration = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConcentration);
  pandoraNuSelection->selShowerCoreHaloRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kCoreHaloRatio);
  pandoraNuSelection->selShowerConicalness = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kConicalness);
  pandoraNuSelection->selShowerdEdxBestPlane = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kdEdxBestPlane);
  pandoraNuSelection->selShowerDisplacement = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDisplacement);
  pandoraNuSelection->selShowerDCA = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kDCA);
  pandoraNuSelection->selShowerWideness = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kWideness);
  pandoraNuSelection->selShowerEnergyDensity = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kEnergyDensity);
  pandoraNuSelection->selShowerPathwayLengthMin = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kPathwayLengthMin);
  pandoraNuSelection->selShowerMaxShowerStartPathwayScatteringAngle2D = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxShowerStartPathwayScatteringAngle2D);
  pandoraNuSelection->selShowerMaxNPostShowerStartHits = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxNPostShowerStartHits);
  pandoraNuSelection->selShowerMaxPostShowerStartScatterAngle = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartScatterAngle);
  pandoraNuSelection->selShowerMaxPostShowerStartNuVertexEnergyAsymmetry = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyAsymmetry);
  pandoraNuSelection->selShowerMaxPostShowerStartShowerStartEnergyAsymmetry = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartShowerStartEnergyAsymmetry);
  pandoraNuSelection->selShowerMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = 
    pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance);
  pandoraNuSelection->selShowerMinPostShowerStartShowerStartMoliereRadius = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinPostShowerStartShowerStartMoliereRadius);
  pandoraNuSelection->selShowerMaxPostShowerStartOpeningAngle = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxPostShowerStartOpeningAngle);
  pandoraNuSelection->selShowerMaxFoundHitRatio = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxFoundHitRatio);
  pandoraNuSelection->selShowerMaxInitialGapSize = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMaxInitialGapSize);
  pandoraNuSelection->selShowerMinLargestProjectedGapSize = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kMinLargestProjectedGapSize);
  pandoraNuSelection->selShowerNViewsWithAmbiguousHits = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kNViewsWithAmbiguousHits);
  pandoraNuSelection->selShowereAmbiguousHitMaxUnaccountedEnergy = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kAmbiguousHitMaxUnaccountedEnergy);

  return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

void PandoraNuSelection::SetNullPandrizzleScores(std::unique_ptr<pandoranusel::PandoraNuSelection> &pandoraNuSelection)
{
  pandoraNuSelection->selShowerBackupPandrizzleScore = kDefDoub;
  pandoraNuSelection->selShowerEnhancedPandrizzleScore = kDefDoub;
  pandoraNuSelection->selShowerEvalRatio = kDefDoub;
  pandoraNuSelection->selShowerConcentration = kDefDoub;
  pandoraNuSelection->selShowerCoreHaloRatio = kDefDoub;
  pandoraNuSelection->selShowerConicalness = kDefDoub;
  pandoraNuSelection->selShowerdEdxBestPlane = kDefDoub;
  pandoraNuSelection->selShowerDisplacement = kDefDoub;
  pandoraNuSelection->selShowerDCA = kDefDoub;
  pandoraNuSelection->selShowerWideness = kDefDoub;
  pandoraNuSelection->selShowerEnergyDensity = kDefDoub;
  pandoraNuSelection->selShowerPathwayLengthMin = kDefDoub;
  pandoraNuSelection->selShowerMaxShowerStartPathwayScatteringAngle2D = kDefDoub;
  pandoraNuSelection->selShowerMaxNPostShowerStartHits = kDefDoub;
  pandoraNuSelection->selShowerMaxPostShowerStartScatterAngle = kDefDoub;
  pandoraNuSelection->selShowerMaxPostShowerStartNuVertexEnergyAsymmetry = kDefDoub;
  pandoraNuSelection->selShowerMaxPostShowerStartShowerStartEnergyAsymmetry = kDefDoub;
  pandoraNuSelection->selShowerMaxPostShowerStartNuVertexEnergyWeightedMeanRadialDistance = kDefDoub; 
  pandoraNuSelection->selShowerMinPostShowerStartShowerStartMoliereRadius = kDefDoub;
  pandoraNuSelection->selShowerMaxPostShowerStartOpeningAngle = kDefDoub;
  pandoraNuSelection->selShowerMaxFoundHitRatio = kDefDoub;
  pandoraNuSelection->selShowerMaxInitialGapSize = kDefDoub;
  pandoraNuSelection->selShowerMinLargestProjectedGapSize = kDefDoub;
  pandoraNuSelection->selShowerNViewsWithAmbiguousHits = kDefDoub;
  pandoraNuSelection->selShowereAmbiguousHitMaxUnaccountedEnergy = kDefDoub;

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
