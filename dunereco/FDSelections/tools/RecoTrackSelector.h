#ifndef TRACKSELECTOR_H_SEEN
#define TRACKSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
//LARSOFT
#include "lardataobj/RecoBase/Track.h"

//Base class for selecting reconstructed tracks in the DUNE far detector

namespace FDSelectionTools{
  class RecoTrackSelector {
    public:
      virtual ~RecoTrackSelector() noexcept = default;
      art::Ptr<recob::Track> FindSelectedTrack(art::Event const & evt) { return SelectTrack(evt); };
    private:
      virtual art::Ptr<recob::Track> SelectTrack(art::Event const & evt) = 0;


  };
}
#endif
