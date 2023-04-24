#ifndef HIGHESTPANDIZZLESCOREVERTEXTRACKSELECTOR_H_SEEN
#define HIGHESTPANDIZZLESCOREVERTEXTRACKSELECTOR_H_SEEN

//STL
#include <iostream>
#include <limits>

//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"

//LArSoft
#include "lardataobj/RecoBase/PFParticle.h"

//DUNE
#include "dunereco/FDSelections/pandizzle/PandizzleAlg.h"

//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class HighestPandizzleScoreRecoVertexTrackSelector : public RecoTrackSelector{
    public:
      explicit HighestPandizzleScoreRecoVertexTrackSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;

      std::string fTrackModuleLabel;
      std::string fPFParticleModuleLabel;

      FDSelection::PandizzleAlg fPandizzleAlg;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestPandizzleScoreRecoVertexTrackSelector)
#endif
