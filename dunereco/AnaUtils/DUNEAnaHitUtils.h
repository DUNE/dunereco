/**
 *
 * @file dune/AnaUtils/DUNEAnaHitUtils.h
 *
 * @brief Utility containing helpful functions for end users to access information about Hits
*/

#ifndef DUNE_ANA_HIT_UTILS_H
#define DUNE_ANA_HIT_UTILS_H

//STL
#include <cmath>
#include <string>
#include <vector>
//ROOT
//ART
#include "art/Framework/Principal/Event.h"
//LARSOFT
#include "lardataobj/RecoBase/Hit.h"
//DUNE
#include "dune/AnaUtils/DUNEAnaUtilsBase.h"

namespace dune_ana
{
/**
 *
 * @brief DUNEAnaHitUtils class
 *
*/
class DUNEAnaHitUtils:DUNEAnaUtilsBase
{
public:
    /**
    * @brief  get the lifetime correction for a hit, assumes the detector properties GetTriggerOffset is T0
    *
    * @param  phit is the hit 
    * 
    * @return the charge normalisation correction
    */
    static double GetLifetimeCorrection(const art::Ptr<recob::Hit> &pHit);

    /**
    * @brief  get the lifetime correction for a particular time
    *
    * @param  timeinTicks the time in ticks
    * @param  t0InMicroS the t0 time in micro seconds
    * 
    * @return the charge normalisation correction
    */
    static double GetLifetimeCorrection(const double timeInTicks, const double t0InMicroS);
};

} // namespace dune_ana


#endif // DUNE_ANA_HIT_UTILS_H

