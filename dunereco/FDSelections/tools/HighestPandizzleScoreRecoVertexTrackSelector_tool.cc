#include "HighestPandizzleScoreRecoVertexTrackSelector.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

FDSelectionTools::HighestPandizzleScoreRecoVertexTrackSelector::HighestPandizzleScoreRecoVertexTrackSelector(fhicl::ParameterSet const& ps) :
    fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
    fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")),
    fPandizzleAlg(ps)
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Track> FDSelectionTools::HighestPandizzleScoreRecoVertexTrackSelector::SelectTrack(art::Event const & evt)
{
  art::Ptr<recob::Track> selTrack;

  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return selTrack;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel);

  // Loop over primaries, get the associated track and then find that with the highest pandizzle score
  double highestPandizzleScore(std::numeric_limits<double>::lowest());

  for (art::Ptr<recob::PFParticle> childPFP : nuChildren) 
  {
    if (!dune_ana::DUNEAnaPFParticleUtils::IsTrack(childPFP, evt, fPFParticleModuleLabel, fTrackModuleLabel))
      continue;

    art::Ptr<recob::Track> childTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(childPFP, evt, fPFParticleModuleLabel, fTrackModuleLabel);

    FDSelection::PandizzleAlg::Record pandizzleRecord(fPandizzleAlg.RunPID(childTrack, evt));
    const double pandizzleScore = pandizzleRecord.GetMVAScore();

    if (pandizzleScore > highestPandizzleScore)
    {
        highestPandizzleScore = pandizzleScore;
        selTrack = childTrack;
    }
  }

  return selTrack;
}
