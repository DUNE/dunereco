#ifndef HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTOR_H_SEEN
#define HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTOR_H_SEEN

//STL
#include <iostream>

//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"

//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class HighestPandrizzleScoreRecoVertexShowerSelector : public RecoShowerSelector{
    public:
      explicit HighestPandrizzleScoreRecoVertexShowerSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;

      std::string fShowerModuleLabel;
      std::string fPFParticleModuleLabel;
      FDSelection::PandrizzleAlg fPandrizzleAlg;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector)
#endif
