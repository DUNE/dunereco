/**
 *  @file   dunetpc/DUNEPandora/DUNE35tTransformationPlugin.cxx
 *
 *  @brief  Implementation of the DUNE35t transformation plugin class.
 *
 *  $Log: $
 */

#include "dune/DUNEPandora/DUNE35tTransformationPlugin.h"

#include <cmath>

namespace lar_pandora
{

DUNE35tTransformationPlugin::DUNE35tTransformationPlugin(const bool isForward) :
    lar_content::LArRotationalTransformationPlugin( (isForward ? M_PI * (45.707f / 180.f) : -M_PI * (45.707f/ 180.f)),
        (isForward ? M_PI * (44.275f/ 180.f) : -M_PI * (44.275f / 180.f)), 1.33f)
{
}

} // namespace lar_pandora
