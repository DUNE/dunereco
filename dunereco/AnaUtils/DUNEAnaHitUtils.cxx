/**
*
* @file dune/AnaUtils/DUNEAnaHitUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about hits
*/

//STL
//ROOT
//ART
#include "art/Framework/Services/Registry/ServiceHandle.h"
//LARSOFT
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
//DUNE
#include "dune/AnaUtils/DUNEAnaHitUtils.h"

namespace dune_ana
{

double DUNEAnaHitUtils::LifetimeCorrection(const art::Ptr<recob::Hit> &pHit)
{
    auto clocks(art::ServiceHandle<detinfo::DetectorClocksService const>()->provider());
    return DUNEAnaHitUtils::LifetimeCorrection(pHit->PeakTime(), clocks->TriggerTime());
}

double DUNEAnaHitUtils::LifetimeCorrection(const double timeInTicks, const double t0InMicroS)
{
    auto clocks(art::ServiceHandle<detinfo::DetectorClocksService const>()->provider());
    auto detProp(art::ServiceHandle<detinfo::DetectorPropertiesService const>()->provider());

    const double tpcSamplingRateInMicroS(clocks->TPCClock().TickPeriod());
    const double tpcTriggerOffsetInTicks(detProp->TriggerOffset());

    const double timeCorrectedForT0((timeInTicks - tpcTriggerOffsetInTicks)*tpcSamplingRateInMicroS - t0InMicroS);
    const double tauLifetime(detProp->ElectronLifetime());

    return std::exp(timeCorrectedForT0/tauLifetime);
}


} // namespace dune_ana


