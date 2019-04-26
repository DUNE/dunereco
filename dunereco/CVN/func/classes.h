#include "dune/CVN/func/TrainingData.h"
#include "dune/CVN/func/PixelMap.h"
#include "dune/CVN/func/GCNGraph.h"
#include "dune/CVN/func/GCNGraphNode.h"
#include "dune/CVN/func/Result.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Wrapper.h"

template class std::vector<cvn::GCNGraph>;
template class art::Ptr<cvn::GCNGraph>;
template class std::vector<art::Ptr<cvn::GCNGraph> >;
template class art::Wrapper< std::vector<cvn::GCNGraph> >;

template class std::vector<cvn::GCNGraphNode>;
template class art::Ptr<cvn::GCNGraphNode>;
template class std::vector<art::Ptr<cvn::GCNGraphNode> >;
template class art::Wrapper< std::vector<cvn::GCNGraphNode> >;

template class std::vector<cvn::Result>;
template class art::Ptr<cvn::Result>;
template class std::vector<art::Ptr<cvn::Result> >;



