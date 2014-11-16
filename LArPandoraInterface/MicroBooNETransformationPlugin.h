/**
 *  @file   larpandora/LArPandoraInterface/MicroBooNETransformationPlugin.h
 *
 *  @brief  Header file for the MicroBooNE transformation plugin class.
 *
 *  $Log: $
 */
#ifndef LAR_PANDORA_MICROBOONE_TRANSFORMATION_PLUGIN_H
#define LAR_PANDORA_MICROBOONE_TRANSFORMATION_PLUGIN_H 1

#include "LArPandoraTransformationPlugin.h"

namespace lar_pandora
{

/**
 *  @brief  MicroBooNETransformationPlugin class
 */
class MicroBooNETransformationPlugin : public LArPandoraTransformationPlugin
{
public:
    /**
     *  @brief  Default constructor
     */
    MicroBooNETransformationPlugin();
};

} // namespace lar_pandora

#endif // #ifndef LAR_PANDORA_MICROBOONE_TRANSFORMATION_PLUGIN_H
