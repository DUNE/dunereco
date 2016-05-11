/**
 *  @file   dunetpc/DUNEPandora/ProtoDUNETransformationPlugin.h
 *
 *  @brief  Header file for the ProtoDUNE transformation plugin class.
 *
 *  $Log: $
 */
#ifndef PROTO_DUNE_TRANSFORMATION_PLUGIN_H
#define PROTO_DUNE_TRANSFORMATION_PLUGIN_H 1

#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  ProtoDUNETransformationPlugin class
 */
class ProtoDUNETransformationPlugin : public lar_content::LArRotationalTransformationPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    ProtoDUNETransformationPlugin(const bool isForward);
};

} // namespace lar_pandora

#endif // #ifndef PROTO_DUNE_TRANSFORMATION_PLUGIN_H
