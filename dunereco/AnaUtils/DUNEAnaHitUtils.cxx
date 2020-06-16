/**
*
* @file dune/AnaUtils/DUNEAnaHitUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about hits
*/

#include "dune/AnaUtils/DUNEAnaHitUtils.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocks.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include <iostream>
namespace dune_ana
{

double DUNEAnaHitUtils::GetLifetimeCorrection(const art::Ptr<recob::Hit> &pHit)
{
    auto clocks(art::ServiceHandle<detinfo::DetectorClocksService const>()->provider());
    std::cout<<"trig time: " << clocks->TriggerTime() << "   beam time: " << clocks->BeamGateTime() << std::endl;

    return DUNEAnaHitUtils::GetLifetimeCorrection(pHit->PeakTime(), clocks->TriggerTime());
}

double DUNEAnaHitUtils::GetLifetimeCorrection(const double timeInTicks, const double t0InMicroS)
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


