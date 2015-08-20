/**
 *  @file  dunetpc/DUNEPandora/DUNE35tGeometryHelper.h
 *
 *  @brief helper function for DUNE 35t geometry
 *
 */
#ifndef DUNE_35T_PANDORA_HELPER_H
#define DUNE_35T_PANDORA_HELPER_H

namespace lar_pandora 
{

class DUNE35tGeometryHelper 
{
public:

    enum DUNE35tVolume
    {
        kShortVolume = 0,
        kLongVolume = 1,
        kUnknownVolume = 2
    };

    /**
     *  @brief Assign a drift volume ID based on cryostate and TPC
     *
     *  @param cstat the cryostat
     *  @param tpc the tpc 
     */
     static DUNE35tGeometryHelper::DUNE35tVolume GetVolumeID(const unsigned int cstat, const unsigned int tpc);
};

} // namespace lar_pandora

#endif //  DUNE_35T_PANDORA_HELPER_H
