/**
 *  @file  dunetpc/DUNEPandora/DUNE4APAGeometryHelper.cxx
 *
 *  @brief helper function for DUNE 4APA geometry
 *
 */

#include "dune/DUNEPandora/DUNE4APAGeometryHelper.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"

#include "cetlib/exception.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

DUNE4APAGeometryHelper::DUNE4APAVolume DUNE4APAGeometryHelper::GetVolumeID(const unsigned int cstat, const unsigned int tpc)
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc,cstat));

    // Left drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        return DUNE4APAGeometryHelper::kLeftVolume;
    }

    // Right drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        return DUNE4APAGeometryHelper::kRightVolume;
    }

    throw cet::exception("LArPandora") << " DUNE4APAGeometryHelper::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
