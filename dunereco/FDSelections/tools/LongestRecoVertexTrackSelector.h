#ifndef LONGESTVERTEXTRACKSELECTOR_H_SEEN
#define LONGESTVERTEXTRACKSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"


//LArSoft
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
//CUSTOM
#include "RecoTrackSelector.h"

namespace FDSelectionTools{
  class LongestRecoVertexTrackSelector : public RecoTrackSelector{
    public:
      explicit LongestRecoVertexTrackSelector(fhicl::ParameterSet const& ps) 
        :
        fTrackModuleLabel(ps.get< std::string> ("ModuleLabels.TrackModuleLabel")),
        fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")) {};


    private:
      art::Ptr<recob::Track> SelectTrack(art::Event const & evt) override;
      std::string fTrackModuleLabel;
      std::string fPFParticleModuleLabel;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::LongestRecoVertexTrackSelector)
#endif
