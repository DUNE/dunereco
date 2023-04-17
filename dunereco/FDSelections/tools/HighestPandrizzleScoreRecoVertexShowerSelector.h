#ifndef HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTOR_H_SEEN
#define HIGHESTPANDRIZZLESCORERECOVERTEXSHOWERSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"

//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class HighestPandrizzleScoreRecoVertexShowerSelector : public RecoShowerSelector{
    public:
      explicit HighestPandrizzleScoreRecoVertexShowerSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;
      std::vector<art::Ptr<recob::Hit> > GetPFPHits(const art::Ptr<recob::PFParticle> pfp, const art::Event& evt);

      std::string fNuGenModuleLabel;
      std::string fClusterModuleLabel;
      std::string fShowerModuleLabel;
      std::string fCheatShowerModuleLabel;
      std::string fPFParticleModuleLabel;
      bool fCheatCharacterisation;
      unsigned int fPFParticleHitCut;
      bool fDemandBDTScore;
      FDSelection::PandrizzleAlg fPandrizzleAlg;
      std::vector<int> fShowerPDGToCheat;

      // Pandrizzle Stuff
      double fRecoNuVtxX;
      double fRecoNuVtxY;
      double fRecoNuVtxZ;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestPandrizzleScoreRecoVertexShowerSelector)
#endif
