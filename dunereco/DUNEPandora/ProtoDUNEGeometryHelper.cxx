/**
 *  @file  dunetpc/DUNEPandora/ProtoDUNEGeometryHelper.cxx
 *
 *  @brief helper function for ProtoDUNE geometry
 *
 */

#include "dune/DUNEPandora/ProtoDUNEGeometryHelper.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"

#include "cetlib/exception.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

ProtoDUNEGeometryHelper::ProtoDUNEVolume ProtoDUNEGeometryHelper::GetVolumeID(const unsigned int cstat, const unsigned int tpc)
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc,cstat));

    // Left drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        return ProtoDUNEGeometryHelper::kLeftVolume;
    }

    // Right drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        return ProtoDUNEGeometryHelper::kRightVolume;
    }

    throw cet::exception("LArPandora") << " ProtoDUNEGeometryHelper::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
