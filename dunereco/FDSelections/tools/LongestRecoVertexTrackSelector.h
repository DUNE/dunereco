#ifndef LONGESTVERTEXTRACKSELECTOR_H_SEEN
#define LONGESTVERTEXTRACKSELECTOR_H_SEEN

//STL
#include <iostream>

//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"

//LArSoft
#include "lardataobj/RecoBase/PFParticle.h"

//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class LongestRecoVertexTrackSelector : public RecoTrackSelector{
    public:
    explicit LongestRecoVertexTrackSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;
      std::string fTrackModuleLabel;
      std::string fPFParticleModuleLabel;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::LongestRecoVertexTrackSelector)
#endif
