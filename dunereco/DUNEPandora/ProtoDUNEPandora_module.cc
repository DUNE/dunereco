/**
 *  @file   dunetpc/DUNEPandora/ProtoDUNEPandora_module.cc
 *
 *  @brief  Producer module for ProtoDUNE detector.
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
 *  @brief  ProtoDUNEPandora class
 */
class ProtoDUNEPandora : public LArPandoraParticleCreator
{
public: 

    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    ProtoDUNEPandora(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~ProtoDUNEPandora();

private:

    unsigned int GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const;
    void ConfigurePandoraGeometry() const;

    bool            m_useLeftVolume;      ///<
    bool            m_useRightVolume;       ///<
};

DEFINE_ART_MODULE(ProtoDUNEPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

// LArSoft includes
#include "larcore/Geometry/Geometry.h"

// Local includes (LArContent) 
#include "LArContent.h"

// Local includes (LArPandora)
#include "dune/DUNEPandora/ProtoDUNEPseudoLayerPlugin.h"
#include "dune/DUNEPandora/ProtoDUNETransformationPlugin.h"
#include "dune/DUNEPandora/ProtoDUNEGeometryHelper.h"

namespace lar_pandora {

ProtoDUNEPandora::ProtoDUNEPandora(fhicl::ParameterSet const &pset) : LArPandoraParticleCreator(pset)
{
    m_useLeftVolume = pset.get<bool>("UseLeftVolume",true);
    m_useRightVolume = pset.get<bool>("UseRightVolume",true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ProtoDUNEPandora::~ProtoDUNEPandora()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProtoDUNEPandora::ConfigurePandoraGeometry() const
{
    mf::LogDebug("LArPandora") << " *** ProtoDUNEPandora::ConfigurePandoraGeometry(...) *** " << std::endl;

    // Identify the Geometry and load the plugins
    art::ServiceHandle<geo::Geometry> theGeometry;

    if (std::string::npos == theGeometry->DetectorName().find("protodune"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " ProtoDUNEPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    for (PandoraInstanceMap::const_iterator pIter = m_pandoraInstanceMap.begin(), pIterEnd = m_pandoraInstanceMap.end(); 
        pIter != pIterEnd; ++pIter)
    {
        const unsigned int      volumeID = pIter->first;
        const pandora::Pandora *pPandora = pIter->second;

        const bool isForward((0 == volumeID) ? true : false); // ATTN: Sign of rotation matrix is taken from Volume ID
    
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, 
            new ProtoDUNEPseudoLayerPlugin));

        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, 
            new ProtoDUNETransformationPlugin(isForward)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int ProtoDUNEPandora::GetPandoraVolumeID(const unsigned int cstat, const unsigned int tpc) const
{    
    const ProtoDUNEGeometryHelper::ProtoDUNEVolume volumeID(ProtoDUNEGeometryHelper::GetVolumeID(cstat, tpc));

    if (ProtoDUNEGeometryHelper::kLeftVolume == volumeID) 
    {
        if (m_useLeftVolume) 
            return 0;
    }

    if (ProtoDUNEGeometryHelper::kRightVolume == volumeID) 
    {
        if (m_useRightVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

} // namespace lar_pandora
