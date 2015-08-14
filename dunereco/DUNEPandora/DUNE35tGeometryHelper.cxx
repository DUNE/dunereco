/**
 *  @file  dunetpc/DUNEPandora/DUNE35tGeometryHelper.cxx
 *
 *  @brief helper function for DUNE 35t geometry
 *
 */

#include "dune/DUNEPandora/DUNE35tGeometryHelper.h"

#include "Geometry/TPCGeo.h"
#include "Geometry/Geometry.h"

#include "cetlib/exception.h"

#include <limits>
#include <iostream>

namespace lar_pandora {

DUNE35tGeometryHelper::DUNE35tVolume DUNE35tGeometryHelper::GetVolumeID(const unsigned int cstat, const unsigned int tpc)
{
    art::ServiceHandle<geo::Geometry> theGeometry;

    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc,cstat));

    // Long drift volume: negative drift direction, odd TPC numbers (1 == tpc%2)
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        return DUNE35tGeometryHelper::kLongVolume;
    }

    // Short drift volume: positive drift direction, even TPC numbers (0 == tpc%2)
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        return DUNE35tGeometryHelper::kShortVolume;
    }

    throw cet::exception("LArPandora") << " DUNE35tGeometryHelper::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
