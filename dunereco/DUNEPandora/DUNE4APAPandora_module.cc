/**
 *  @file   dunetpc/DUNEPandora/DUNE4APAPandora_module.cc
 *
 *  @brief  LArPandora producer module for DUNE4APA detector.
 *
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "Geometry/Geometry.h"

#include "LArStitching/MultiPandoraApi.h"

#include "LArPandoraInterface/LArPandora.h"

#include <set>
#include <string>

namespace lar_pandora
{

/**
 *  @brief  DUNE4APAPandora class
 */
class DUNE4APAPandora : public LArPandora
{
public: 
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    DUNE4APAPandora(fhicl::ParameterSet const &pset);

    int GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const;

private:
    void CreatePandoraInstances();

    /**
     *  @brief  Create primary pandora instance
     *
     *  @param  stitchingConfigFileName the pandora settings stitching config file name
     */
    void CreatePrimaryPandoraInstance(const std::string &stitchingConfigFileName);

    /**
     *  @brief  Create daughter pandora instances
     *
     *  @param  theGeometry the geometry handle
     *  @param  configFileName the pandora settings config file name
     */
    void CreateDaughterPandoraInstances(const art::ServiceHandle<geo::Geometry> &theGeometry, const std::string &configFileName);

    typedef std::set<int> IntSet;   ///<

    bool    m_useLeftVolume;        ///<
    bool    m_useRightVolume;       ///<
};

DEFINE_ART_MODULE(DUNE4APAPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib/exception.h"

#include "Api/PandoraApi.h"

#include "LArContent.h"

#include "dune/DUNEPandora/DUNE4APAPseudoLayerPlugin.h"
#include "dune/DUNEPandora/DUNE4APATransformationPlugin.h"

namespace lar_pandora
{

DUNE4APAPandora::DUNE4APAPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
    m_useLeftVolume = pset.get<bool>("UseLeftVolume", true);
    m_useRightVolume = pset.get<bool>("UseRightVolume", true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DUNE4APAPandora::GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;
    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc, cryostat));

    // Left drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        if (m_useLeftVolume) 
            return 0;
    }

    // Right drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        if (m_useRightVolume) 
            return 1;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE4APAPandora::CreatePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** DUNE4APAPandora::CreatePandoraInstances(...) *** " << std::endl;

    art::ServiceHandle<geo::Geometry> theGeometry;
    if (std::string::npos == theGeometry->DetectorName().find("dune10kt"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " DUNE4APAPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
    }

    cet::search_path sp("FW_SEARCH_PATH");
    std::string stitchingConfigFileName, configFileName;
    if (!sp.find_file(m_stitchingConfigFile, stitchingConfigFileName) || !sp.find_file(m_configFile, configFileName))
    {                                               
        mf::LogError("LArPandora") << "   Failed to find one of: " << m_stitchingConfigFile << ", " << m_configFile << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);              
    }

    this->CreatePrimaryPandoraInstance(stitchingConfigFileName);
    this->CreateDaughterPandoraInstances(theGeometry, configFileName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE4APAPandora::CreatePrimaryPandoraInstance(const std::string &stitchingConfigFileName)
{
    m_pPrimaryPandora = this->CreateNewPandora();
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, stitchingConfigFileName));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora, new lar_pandora::DUNE4APAPseudoLayerPlugin));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE4APAPandora::CreateDaughterPandoraInstances(const art::ServiceHandle<geo::Geometry> &theGeometry, const std::string &configFileName)
{
    if (!m_pPrimaryPandora)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    IntSet volumeIdNumbers;

    for (unsigned int icstat = 0; icstat < theGeometry->Ncryostats(); ++icstat)
    {
        for (unsigned int itpc = 0; itpc < theGeometry->NTPC(icstat); ++itpc)
        {
            try
            {
                const int volumeIdNumber(this->GetVolumeIdNumber(icstat, itpc));

                if (!volumeIdNumbers.insert(volumeIdNumber).second)
                    continue;

                const bool isForward(0 == volumeIdNumber); // ATTN: Sign of rotation matrix is taken from Volume ID
                const bool isPositiveDrift(1 == volumeIdNumber);

                std::ostringstream volumeIdString("dune10kt_");
                volumeIdString << volumeIdNumber;

                const pandora::Pandora *const pPandora = this->CreateNewPandora();
                MultiPandoraApi::AddDaughterPandoraInstance(m_pPrimaryPandora, pPandora);
                MultiPandoraApi::SetVolumeInfo(pPandora, new VolumeInfo(volumeIdNumber, volumeIdString.str(), pandora::CartesianVector(0.f, 0.f, 0.f), isPositiveDrift));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, new lar_pandora::DUNE4APAPseudoLayerPlugin));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, new lar_pandora::DUNE4APATransformationPlugin(isForward)));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, configFileName));
            }
            catch (pandora::StatusCodeException &)
            {
                mf::LogDebug("DUNE4APAPandora") << "    No volume ID for this TPC..." << std::endl;
            }
        }
    }
}

} // namespace lar_pandora
