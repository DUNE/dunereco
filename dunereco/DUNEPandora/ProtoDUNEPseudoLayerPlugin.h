/**
 *  @file   dunetpc/DUNEPandora/ProtoDUNEPseudoLayerPlugin.h
 *
 *  @brief  Header file for the ProtoDUNE pseudo layer plugin class.
 *
 *  $Log: $
 */
#ifndef PROTO_DUNE_PSEUDO_LAYER_PLUGIN_H
#define PROTO_DUNE_PSEUDO_LAYER_PLUGIN_H 1

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  ProtoDUNEPseudoLayerPlugin class
 */
class ProtoDUNEPseudoLayerPlugin : public lar_content::LArPseudoLayerPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    ProtoDUNEPseudoLayerPlugin();
};

} // namespace lar

#endif // #ifndef PROTO_DUNE_PSEUDO_LAYER_PLUGIN_H
