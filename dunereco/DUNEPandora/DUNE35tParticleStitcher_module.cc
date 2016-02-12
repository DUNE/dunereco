/**
 *  @file   dunetpc/DUNEPandora/DUNE35tParticleStitcher_module.cc
 *
 *  @brief  Stitching module for DUNE 35t detector (soon to be replaced with stitching by algorithms running within a Pandora instance)
 *
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "TVector3.h" // ATTN Should be included in PFParticleStitcher.h
#include "LArPandoraInterface/PFParticleStitcher.h"

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
     *  @param  pset the parameter set
     */
    DUNE35tParticleStitcher(fhicl::ParameterSet const &pset);

private:
    unsigned int GetVolumeID(const unsigned int cstat, const unsigned int tpc) const;
};

DEFINE_ART_MODULE(DUNE35tParticleStitcher)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib/exception.h"

#include "Geometry/Geometry.h"

namespace lar_pandora
{

DUNE35tParticleStitcher::DUNE35tParticleStitcher(fhicl::ParameterSet const &pset) :
    PFParticleStitcher(pset)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DUNE35tParticleStitcher::GetVolumeID(const unsigned int cstat, const unsigned int tpc) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;
    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc, cstat));

    // Long drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
        return 1;

    // Short drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
        return 0;

    throw cet::exception("LArPandora") << " DUNE35tParticleStitcher::GetVolumeID --- No valid ID for cstat=" << cstat << ", tpc=" << tpc;
}

} // namespace lar_pandora
