/**
 *  @file   dunetpc/DUNEPandora/DUNE35tPseudoLayerPlugin.h
 *
 *  @brief  Header file for the DUNE35t pseudo layer plugin class.
 *
 *  $Log: $
 */
#ifndef DUNE_35T_PSEUDO_LAYER_PLUGIN_H
#define DUNE_35T_PSEUDO_LAYER_PLUGIN_H 1

#include "LArPlugins/LArPseudoLayerPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  DUNE35tPseudoLayerPlugin class
 */
class DUNE35tPseudoLayerPlugin : public lar_content::LArPseudoLayerPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    DUNE35tPseudoLayerPlugin();
};

} // namespace lar

#endif // #ifndef DUNE_35T_PSEUDO_LAYER_PLUGIN_H
