/**
 *  @file  dunetpc/DUNEPandora/DUNE4APAGeometryHelper.h
 *
 *  @brief helper function for DUNE 4APA geometry
 *
 */
#ifndef DUNE_4APA_PANDORA_HELPER_H
#define DUNE_4APA_PANDORA_HELPER_H

namespace lar_pandora 
{

class DUNE4APAGeometryHelper 
{
public:

    enum DUNE4APAVolume
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
     static DUNE4APAGeometryHelper::DUNE4APAVolume GetVolumeID(const unsigned int cstat, const unsigned int tpc);
};

} // namespace lar_pandora

#endif //  LAR_PANDORA_DUNE_4APA_PANDORA_HELPER_H
