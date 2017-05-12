#include "dune/CVN/func/TrainingData.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/Result.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Wrapper.h"

template class std::vector<float>;
template class std::vector<cvn::PixelMap>;
template class art::Ptr<cvn::PixelMap>;
template class std::vector<art::Ptr<cvn::PixelMap> >;

template class art::Wrapper< std::vector<cvn::PixelMap> >;

template class std::vector<cvn::Result>;
template class art::Ptr<cvn::Result>;
template class std::vector<art::Ptr<cvn::Result> >;

template class art::Wrapper< std::vector<cvn::Result> >;


