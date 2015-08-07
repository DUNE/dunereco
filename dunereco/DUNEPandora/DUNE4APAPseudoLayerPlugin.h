/**
 *  @file   dunetpc/DUNEPandora/DUNE4APAPseudoLayerPlugin.h
 *
 *  @brief  Header file for the DUNE4APA pseudo layer plugin class.
 *
 *  $Log: $
 */
#ifndef DUNE_4APA_PSEUDO_LAYER_PLUGIN_H
#define DUNE_4APA_PSEUDO_LAYER_PLUGIN_H 1

#include "LArPlugins/LArPseudoLayerPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  DUNE4APAPseudoLayerPlugin class
 */
class DUNE4APAPseudoLayerPlugin : public lar_content::LArPseudoLayerPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    DUNE4APAPseudoLayerPlugin();
};

} // namespace lar

#endif // #ifndef DUNE_4APA_PSEUDO_LAYER_PLUGIN_H
