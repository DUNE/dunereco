#include "utils.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "dune/AnaUtils/DUNEAnaHitUtils.h"

namespace VLN {

std::pair<double, double> calcHitsChargeCalE(
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const art::Event                        &evt,
    calo::CalorimetryAlg                    &algCalorimetry,
    unsigned int                            plane
)
{
    const auto clockData =
        art::ServiceHandle<const detinfo::DetectorClocksService>()
            ->DataFor(evt);

    const auto detProp =
        art::ServiceHandle<const detinfo::DetectorPropertiesService>()
            ->DataFor(evt, clockData);

    const double charge =
        dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(
            clockData, detProp, hits
        );

    const double calE = algCalorimetry.ElectronsFromADCArea(charge, plane);

    return std::make_pair(charge, calE);
}

}

