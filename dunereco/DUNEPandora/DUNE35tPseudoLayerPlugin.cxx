/**
 *  @file   dunetpc/DUNEPandora/DUNE35tPseudoLayerPlugin.cxx
 *
 *  @brief  Implementation of the DUNE35t pseudo layer Plugin class.
 *
 *  $Log: $
 */

#include "dune/DUNEPandora/DUNE35tPseudoLayerPlugin.h"

namespace lar_pandora
{

DUNE35tPseudoLayerPlugin::DUNE35tPseudoLayerPlugin() : LArPseudoLayerPlugin(0.488f, 0.500f, 0.449f)
{
}

} // namespace lar_pandora
