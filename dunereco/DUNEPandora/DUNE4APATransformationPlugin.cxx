/**
 *  @file   dunetpc/DUNEPandora/DUNE4APATransformationPlugin.cxx
 *
 *  @brief  Implementation of the DUNE4APA transformation plugin class.
 *
 *  $Log: $
 */

#include "dune/DUNEPandora/DUNE4APATransformationPlugin.h"

#include <cmath>

namespace lar_pandora
{

DUNE4APATransformationPlugin::DUNE4APATransformationPlugin(const bool isForward) :
    lar_content::LArRotationalTransformationPlugin( (isForward ? M_PI * (36.f / 180.f) : -M_PI * (36.f/ 180.f)),
        (isForward ? M_PI * (36.f/ 180.f) : -M_PI * (36.f / 180.f)), 1.33f)
{
}

} // namespace lar_pandora
