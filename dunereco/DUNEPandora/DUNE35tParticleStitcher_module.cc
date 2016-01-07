/**
 *  @file   dunetpc/DUNEPandora/DUNE35tParticleStitcher_module.cc
 *
 *  @brief  Stitching module for DUNE 35t detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "larpandora/LArPandoraInterface/PFParticleStitcher.h"

#include "cetlib/exception.h"

#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  DUNE35tParticleStitcher class
 */
class DUNE35tParticleStitcher : public PFParticleStitcher
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    DUNE35tParticleStitcher(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~DUNE35tParticleStitcher();

private:
    unsigned int GetVolumeID(const unsigned int cstat, const unsigned int tpc) const;

};

DEFINE_ART_MODULE(DUNE35tParticleStitcher)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "dune/DUNEPandora/DUNE35tGeometryHelper.h"

namespace lar_pandora {

DUNE35tParticleStitcher::DUNE35tParticleStitcher(fhicl::ParameterSet const &pset) : PFParticleStitcher(pset)
{
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

DUNE35tParticleStitcher::~DUNE35tParticleStitcher()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DUNE35tParticleStitcher::GetVolumeID(const unsigned int cstat, const unsigned int tpc) const
{
    DUNE35tGeometryHelper::DUNE35tVolume volumeID(DUNE35tGeometryHelper::GetVolumeID(cstat, tpc));

    if (DUNE35tGeometryHelper::kShortVolume == volumeID) 
        return 0;

    if (DUNE35tGeometryHelper::kLongVolume == volumeID) 
        return 1;

    throw cet::exception("LArPandora") << " DUNE35tParticleStitcher::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
