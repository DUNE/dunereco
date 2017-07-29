#include "dunetpc/dune/SNSlicer/SNSlice.h"

#include "canvas/Persistency/Common/Wrapper.h"

#include <vector>

template class std::vector<sn::SNSlice>;

template class art::Wrapper<sn::SNSlice>;
template class art::Wrapper<std::vector<sn::SNSlice>>;
