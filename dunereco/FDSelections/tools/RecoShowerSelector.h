#ifndef SHOWERSELECTOR_H_SEEN
#define SHOWERSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
//LARSOFT
#include "lardataobj/RecoBase/Shower.h"

//Base class for selecting reconstructed tracks in the DUNE far detector

namespace FDSelectionTools{
  class RecoShowerSelector {
    public:
      virtual ~RecoShowerSelector() noexcept = default;
      art::Ptr<recob::Shower> FindSelectedShower(art::Event const & evt) { return SelectShower(evt); };
    private:
      virtual art::Ptr<recob::Shower> SelectShower(art::Event const & evt) = 0;


  };
}
#endif
