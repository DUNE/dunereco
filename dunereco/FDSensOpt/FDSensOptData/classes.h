#include "dune/FDSensOpt/FDSensOptData/MVASelectPID.h"
#include "dune/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Wrapper.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

template class art::Ptr<dunemva::MVASelectPID>;
template class art::Ptr<dune::EnergyRecoOutput>;


