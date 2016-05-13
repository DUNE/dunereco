/**
 *  @file   dunetpc/DUNEPandora/DUNE35tTransformationPlugin.h
 *
 *  @brief  Header file for the DUNE35t transformation plugin class.
 *
 *  $Log: $
 */
#ifndef DUNE_35T_TRANSFORMATION_PLUGIN_H
#define DUNE_35T_TRANSFORMATION_PLUGIN_H 1

#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  DUNE35tTransformationPlugin class
 */
class DUNE35tTransformationPlugin : public lar_content::LArRotationalTransformationPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    DUNE35tTransformationPlugin(const bool isForward);
};

} // namespace lar_pandora

#endif // #ifndef DUNE_35T_TRANSFORMATION_PLUGIN_H
