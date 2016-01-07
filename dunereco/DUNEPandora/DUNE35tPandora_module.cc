/**
 *  @file   dunetpc/DUNEPandora/DUNE35tPandora_module.cc
 *
 *  @brief  Producer module for DUNE 35t detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "larpandora/LArPandoraInterface/LArPandoraParticleCreator.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  DUNE35tPandora class
 */
class DUNE35tPandora : public LArPandoraParticleCreator
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    DUNE35tPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~DUNE35tPandora();

private:

    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;

    bool            m_useShortVolume;      ///<
    bool            m_useLongVolume;       ///<
};

DEFINE_ART_MODULE(DUNE35tPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "larcore/Geometry/Geometry.h"

// Local includes (LArContent) 
#include "LArContent.h"

// Local includes (LArPandora)
#include "dune/DUNEPandora/DUNE35tPseudoLayerPlugin.h"
#include "dune/DUNEPandora/DUNE35tTransformationPlugin.h"
#include "dune/DUNEPandora/DUNE35tGeometryHelper.h"

namespace lar_pandora {

DUNE35tPandora::DUNE35tPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{
    m_useShortVolume = pset.get<bool>("UseShortVolume",true);
    m_useLongVolume = pset.get<bool>("UseLongVolume",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DUNE35tPandora::~DUNE35tPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE35tPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** DUNE35tPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("dune35t"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " DUNE35tPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const unsigned int      volumeID = pIter->first;
        const pandora::Pandora *pPandora = pIter->second;

        const bool isForward((0 == volumeID) ? true : false); // ATTN: Sign of rotation matrix is taken from Volume ID
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new DUNE35tPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new DUNE35tTransformationPlugin(isForward)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DUNE35tPandora::GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const
{    
    const DUNE35tGeometryHelper::DUNE35tVolume volumeID(DUNE35tGeometryHelper::GetVolumeID(cstat, tpc));

    if (DUNE35tGeometryHelper::kShortVolume == volumeID) 
    {
        if (m_useShortVolume) 
            return 0;
    }

    if (DUNE35tGeometryHelper::kLongVolume == volumeID) 
    {
        if (m_useLongVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_pandora
