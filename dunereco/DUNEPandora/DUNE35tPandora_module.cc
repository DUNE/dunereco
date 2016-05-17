/**
 *  @file   dunetpc/DUNEPandora/DUNE35tPandora_module.cc
 *
 *  @brief  LArPandora producer module for DUNE35t detector.
 *
 */

#include "art/Framework/Core/ModuleMacros.h"

#include "larcore/Geometry/Geometry.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include "larpandora/LArPandoraInterface/LArPandora.h"

#include <set>
#include <string>

namespace lar_pandora
{

/**
 *  @brief  DUNE35tPandora class
 */
class DUNE35tPandora : public LArPandora
{
public: 
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    DUNE35tPandora(fhicl::ParameterSet const &pset);

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

    bool    m_useShortVolume;       ///<
    bool    m_useLongVolume;        ///<
};

DEFINE_ART_MODULE(DUNE35tPandora)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "cetlib/exception.h"

#include "Api/PandoraApi.h"

#include "larpandoracontent/LArContent.h"

#include "dune/DUNEPandora/DUNE35tPseudoLayerPlugin.h"
#include "dune/DUNEPandora/DUNE35tTransformationPlugin.h"

namespace lar_pandora
{

DUNE35tPandora::DUNE35tPandora(fhicl::ParameterSet const &pset) :
    LArPandora(pset)
{
    m_useShortVolume = pset.get<bool>("UseShortVolume", true);
    m_useLongVolume = pset.get<bool>("UseLongVolume", true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DUNE35tPandora::GetVolumeIdNumber(const unsigned int cryostat, const unsigned int tpc) const
{
    art::ServiceHandle<geo::Geometry> theGeometry;
    const geo::TPCGeo &theTpcGeo(theGeometry->TPC(tpc, cryostat));

    // Long drift volume: negative drift direction
    if (theTpcGeo.DriftDirection() == geo::kNegX)
    {
        if (m_useLongVolume) 
            return 1;
    }

    // Short drift volume: positive drift direction
    if (theTpcGeo.DriftDirection() == geo::kPosX)
    {
        if (m_useShortVolume) 
            return 0;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE35tPandora::CreatePandoraInstances()
{
    mf::LogDebug("LArPandora") << " *** DUNE35tPandora::CreatePandoraInstances(...) *** " << std::endl;

    art::ServiceHandle<geo::Geometry> theGeometry;
    if (std::string::npos == theGeometry->DetectorName().find("dune35t"))
    {
        mf::LogError("LArPandora") << " Error! Using invalid geometry: " << theGeometry->DetectorName() << std::endl;
        throw cet::exception("LArPandora") << " DUNE35tPandora::ConfigurePandoraGeometry --- Invalid Geometry: " << theGeometry->DetectorName();
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

void DUNE35tPandora::CreatePrimaryPandoraInstance(const std::string &stitchingConfigFileName)
{
    m_pPrimaryPandora = this->CreateNewPandora();
    MultiPandoraApi::AddPrimaryPandoraInstance(m_pPrimaryPandora);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPrimaryPandora, stitchingConfigFileName));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*m_pPrimaryPandora, new lar_pandora::DUNE35tPseudoLayerPlugin));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DUNE35tPandora::CreateDaughterPandoraInstances(const art::ServiceHandle<geo::Geometry> &theGeometry, const std::string &configFileName)
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
                const bool isPositiveDrift(0 == volumeIdNumber);

                std::ostringstream volumeIdString("dune35t_");
                volumeIdString << volumeIdNumber;

                const pandora::Pandora *const pPandora = this->CreateNewPandora();
                MultiPandoraApi::AddDaughterPandoraInstance(m_pPrimaryPandora, pPandora);
                MultiPandoraApi::SetVolumeInfo(pPandora, new VolumeInfo(volumeIdNumber, volumeIdString.str(), pandora::CartesianVector(0.f, 0.f, 0.f), isPositiveDrift));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArPseudoLayerPlugin(*pPandora, new lar_pandora::DUNE35tPseudoLayerPlugin));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LArContent::SetLArTransformationPlugin(*pPandora, new lar_pandora::DUNE35tTransformationPlugin(isForward)));
                PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, configFileName));
            }
            catch (pandora::StatusCodeException &)
            {
                mf::LogDebug("DUNE35tPandora") << "    No volume ID for this TPC..." << std::endl;
            }
        }
    }
}

} // namespace lar_pandora
