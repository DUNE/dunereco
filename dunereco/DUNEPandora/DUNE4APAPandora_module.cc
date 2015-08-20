/**
 *  @file   dunetpc/DUNEPandora/DUNE4APAPandora_module.cc
 *
 *  @brief  Producer module for DUNE 4APA detector.
 *
 */

// Framework Includes
#include "art/Framework/Core/ModuleMacros.h"

// Local includes
#include "LArPandoraInterface/LArPandoraParticleCreator.h"

// std includes
#include <string>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora
{

/**
 *  @brief  DUNE4APAPandora class
 */
class DUNE4APAPandora : public LArPandoraParticleCreator
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    DUNE4APAPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~DUNE4APAPandora();

private:

    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;

    bool            m_useLeftVolume;      ///<
    bool            m_useRightVolume;       ///<
};

DEFINE_ART_MODULE(DUNE4APAPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "Geometry/Geometry.h"

// Local includes (LArContent) 
#include "LArContent.h"

// Local includes (LArPandora)
#include "dune/DUNEPandora/DUNE4APAPseudoLayerPlugin.h"
#include "dune/DUNEPandora/DUNE4APATransformationPlugin.h"
#include "dune/DUNEPandora/DUNE4APAGeometryHelper.h"

namespace lar_pandora {

DUNE4APAPandora::DUNE4APAPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{
    m_useLeftVolume = pset.get<bool>("UseLeftVolume",true);
    m_useRightVolume = pset.get<bool>("UseRightVolume",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

DUNE4APAPandora::~DUNE4APAPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE4APAPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** DUNE4APAPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("dune10kt"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " DUNE4APAPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const unsigned int      volumeID = pIter->first;
        const pandora::Pandora *pPandora = pIter->second;

        const bool isForward((0 == volumeID) ? true : false); // ATTN: Sign of rotation matrix is taken from Volume ID
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new DUNE4APAPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new DUNE4APATransformationPlugin(isForward)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int DUNE4APAPandora::GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const
{    
    const DUNE4APAGeometryHelper::DUNE4APAVolume volumeID(DUNE4APAGeometryHelper::GetVolumeID(cstat, tpc));

    if (DUNE4APAGeometryHelper::kLeftVolume == volumeID) 
    {
        if (m_useLeftVolume) 
            return 0;
    }

    if (DUNE4APAGeometryHelper::kRightVolume == volumeID) 
    {
        if (m_useRightVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_pandora
