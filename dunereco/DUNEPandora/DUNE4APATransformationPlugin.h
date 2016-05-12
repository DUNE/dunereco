/**
 *  @file   dunetpc/DUNEPandora/DUNE4APATransformationPlugin.h
 *
 *  @brief  Header file for the DUNE4APA transformation plugin class.
 *
 *  $Log: $
 */
#ifndef DUNE_4APA_TRANSFORMATION_PLUGIN_H
#define DUNE_4APA_TRANSFORMATION_PLUGIN_H 1

#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  DUNE4APATransformationPlugin class
 */
class DUNE4APATransformationPlugin : public lar_content::LArRotationalTransformationPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    DUNE4APATransformationPlugin(const bool isForward);
};

} // namespace lar_pandora

#endif // #ifndef DUNE_4APA_TRANSFORMATION_PLUGIN_H
