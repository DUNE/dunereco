#ifndef SNSLICE_H
#define SNSLICE_H

#include "canvas/Persistency/Common/Ptr.h"
#include "lardataobj/RecoBase/Hit.h"

namespace sn
{
  class SNSlice
  {
  public:
    float meanT;
    float meanZ;
    float totQ;

    std::vector<art::Ptr<recob::Hit>> hits;
  };
}

#endif
