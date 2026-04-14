#include "LongestRecoVertexTrackSelector.h"

#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"

FDSelectionTools::LongestRecoVertexTrackSelector::LongestRecoVertexTrackSelector(fhicl::ParameterSet const& ps) :
  fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
  fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel"))
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

art::Ptr<recob::Track> FDSelectionTools::LongestRecoVertexTrackSelector::SelectTrack(art::Event const & evt)
{
  art::Ptr<recob::Track> selTrack;

  if (!dune_ana::DUNEAnaEventUtils::HasNeutrino(evt, fPFParticleModuleLabel))
    return selTrack;

  art::Ptr<recob::PFParticle> nuPFP = dune_ana::DUNEAnaEventUtils::GetNeutrino(evt, fPFParticleModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> nuChildren = dune_ana::DUNEAnaPFParticleUtils::GetChildParticles(nuPFP, evt, fPFParticleModuleLabel);

  double longestLength = -999.0;

  for (art::Ptr<recob::PFParticle> childPFP : nuChildren) 
  {
    if (!dune_ana::DUNEAnaPFParticleUtils::HasTrack(childPFP, evt, fPFParticleModuleLabel, fTrackModuleLabel))
      continue;

    art::Ptr<recob::Track> childTrack = dune_ana::DUNEAnaPFParticleUtils::GetTrack(childPFP, evt, fPFParticleModuleLabel, fTrackModuleLabel);

    double childTrackLength = childTrack->Length();

    if (childTrackLength > longestLength)
    {
      longestLength = childTrackLength;
      selTrack = childTrack;
    }
  }

  return selTrack;
}
