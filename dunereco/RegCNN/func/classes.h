#include "dune/RegCNN/func/RegPixelMap.h"
#include "dune/RegCNN/func/RegCNNResult.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Wrapper.h"

template class std::vector<float>;
template class std::vector<cnn::RegPixelMap>;
template class art::Ptr<cnn::RegPixelMap>;
template class std::vector<art::Ptr<cnn::RegPixelMap> >;
template class art::Wrapper< std::vector<cnn::RegPixelMap> >;


template class std::vector<cnn::RegCNNResult>;
template class art::Ptr<cnn::RegCNNResult>;
template class std::vector<art::Ptr<cnn::RegCNNResult> >;

template class art::Wrapper< std::vector<cnn::RegCNNResult> >;


