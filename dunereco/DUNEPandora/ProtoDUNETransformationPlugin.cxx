/**
 *  @file   dunetpc/DUNEPandora/ProtoDUNETransformationPlugin.cxx
 *
 *  @brief  Implementation of the ProtoDUNE transformation plugin class.
 *
 *  $Log: $
 */

#include "dune/DUNEPandora/ProtoDUNETransformationPlugin.h"

#include <cmath>

namespace lar_pandora
{

ProtoDUNETransformationPlugin::ProtoDUNETransformationPlugin(const bool isForward) :
    lar_content::LArRotationalTransformationPlugin( (isForward ? M_PI * (35.71f / 180.f) : -M_PI * (35.71f/ 180.f)),
        (isForward ? M_PI * (35.71f/ 180.f) : -M_PI * (35.71f / 180.f)), 1.33f)
{
}

} // namespace lar_pandora
