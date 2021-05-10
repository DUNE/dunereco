#pragma once

#include <string>
#include <vector>

#include "art/Framework/Principal/Event.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace VLN {

std::pair<double, double> calcHitsChargeCalE(
    const std::vector<art::Ptr<recob::Hit>> &hits,
    const art::Event                        &evt,
    calo::CalorimetryAlg                    &algCalorimetry,
    unsigned int plane = 2
);

}

