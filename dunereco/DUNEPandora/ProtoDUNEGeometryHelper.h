/**
 *  @file  dunetpc/DUNEPandora/ProtoDUNEGeometryHelper.h
 *
 *  @brief helper function for ProtoDUNE geometry
 *
 */
#ifndef PROTO_DUNE_PANDORA_HELPER_H
#define PROTO_DUNE_PANDORA_HELPER_H

namespace lar_pandora 
{

class ProtoDUNEGeometryHelper 
{
public:

    enum ProtoDUNEVolume
    {
        kLeftVolume = 0,
        kRightVolume = 1,
        kUnknownVolume = 2
    };

    /**
     *  @brief Assign a drift volume ID based on cryostate and TPC
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc 
     */
     static ProtoDUNEGeometryHelper::ProtoDUNEVolume GetVolumeID(const unsigned int cstat, const unsigned int tpc);
};

} // namespace lar_pandora

#endif //  PROTO_DUNE_PANDORA_HELPER_H
