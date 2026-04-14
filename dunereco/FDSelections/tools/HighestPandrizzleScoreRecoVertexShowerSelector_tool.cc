#include "HighestPandrizzleScoreRecoVertexShowerSelector.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector::HighestPandrizzleScoreRecoVertexShowerSelector(fhicl::ParameterSet const& ps) :
    fShowerModuleLabel(ps.get< std::string> ("ModuleLabels.ShowerModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fPandrizzleAlg(ps)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Shower> FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector::SelectShower(art::Event const & evt)
{
  art::Ptr<recob::Shower> selShower;
  art::Ptr<recob::Shower> backupSelShower;

  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return selShower;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel);

  double highestEnhancedPandrizzleScore = std::numeric_limits<double>::lowest();
  double highestBackupPandrizzleScore = std::numeric_limits<double>::lowest();
  bool found = false;

  for (art::Ptr<recob::PFParticle> childPFP : nuChildren) 
  {
    if (!dune_ana::DUNEAnaPFParticleUtils::HasShower(childPFP, evt, fPFParticleModuleLabel, fShowerModuleLabel))
      continue;

    art::Ptr<recob::Shower> childShower = dune_ana::DUNEAnaPFParticleUtils::GetShower(childPFP, evt, fPFParticleModuleLabel, fShowerModuleLabel);

    FDSelection::PandrizzleAlg::Record pandrizzleRecord(fPandrizzleAlg.RunPID(childShower, evt));
    float pandrizzleBDTMethod = pandrizzleRecord.GetVar(FDSelection::PandrizzleAlg::kBDTMethod);
    double pandrizzleScore = pandrizzleRecord.GetMVAScore();

    if ((std::fabs(pandrizzleBDTMethod - 2.0) < std::numeric_limits<float>::epsilon()) && (pandrizzleScore > highestEnhancedPandrizzleScore))
    {
        highestEnhancedPandrizzleScore = pandrizzleScore;
        selShower = childShower;
        found = true;
    }

    if ((std::fabs(pandrizzleBDTMethod - 1.0) < std::numeric_limits<float>::epsilon()) && (pandrizzleScore > highestBackupPandrizzleScore))
    {
        highestBackupPandrizzleScore = pandrizzleScore;
        backupSelShower = childShower;
    }
  }

  return (found ? selShower : backupSelShower);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

