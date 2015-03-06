/**
 *  @file   LArContent/include/LArContent.h
 * 
 *  @brief  Header file detailing content for use with particle flow reconstruction at liquid argon time projection chambers
 * 
 *  $Log: $
 */
#ifndef LAR_CONTENT_H
#define LAR_CONTENT_H 1

#include "LArCheating/CheatingClusterCharacterisationAlgorithm.h"
#include "LArCheating/CheatingClusterCreationAlgorithm.h"
#include "LArCheating/CheatingCosmicRayIdentificationAlg.h"
#include "LArCheating/CheatingCosmicRayShowerMatchingAlg.h"
#include "LArCheating/CheatingPfoCreationAlgorithm.h"
#include "LArCheating/CheatingVertexCreationAlgorithm.h"

#include "LArPandoraAlgorithms/LArHelpers/LArGeometryHelper.h"

#include "LArMonitoring/EventDisplayAlgorithm.h"
#include "LArMonitoring/ParticleAnalysisAlgorithm.h"
#include "LArMonitoring/ParticleMonitoringAlgorithm.h"
#include "LArMonitoring/VisualMonitoringAlgorithm.h"

#include "LArPandoraAlgorithms/LArPlugins/LArParticleIdPlugins.h"
#include "LArPandoraAlgorithms/LArPlugins/LArPseudoLayerPlugin.h"
#include "LArPandoraAlgorithms/LArPlugins/LArTransformationPlugin.h"

#include "LArPandoraAlgorithms/LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArCosmicRay/DeltaRayMatchingAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArEventBuilding/NeutrinoEventBuildingAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArEventBuilding/NeutrinoEventCreationAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArEventBuilding/NeutrinoVertexBuildingAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArEventBuilding/NeutrinoVertexCreationAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArHitCreation/TwoViewShowerHitsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArLongitudinalTrackMatching/ClearLongitudinalTracksTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArPfoMopUp/ParticleRecoveryAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerMatching/ClearShowersTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerMatching/ShowerTensorVisualizationTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerMatching/SimpleShowersTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArShowerMatching/SplitShowersTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/OvershootTracksTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/TrackSplittingTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/TransverseTensorVisualizationTool.h"
#include "LArPandoraAlgorithms/LArThreeDReco/LArTransverseTrackMatching/UndershootTracksTool.h"

#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterAssociation/TransverseAssociationAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterAssociation/TransverseExtensionAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterMopUp/IsolatedHitMergingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterMopUp/ProximityBasedMergingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArCosmicRay/CosmicRayExtensionAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArCosmicRay/DeltaRayGrowingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/KinkSplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/LayerSplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/LArSeedFinding/ShowerGrowingAlgorithm.h"
#include "LArPandoraAlgorithms/LArTwoDReco/TwoDParticleCreationAlgorithm.h"

#include "LArUtility/ListChangingAlgorithm.h"
#include "LArUtility/ListDissolutionAlgorithm.h"
#include "LArUtility/ListMergingAlgorithm.h"
#include "LArUtility/ListPreparationAlgorithm.h"

#include "LArVertex/CandidateVertexCreationAlgorithm.h"
#include "LArVertex/VertexSelectionAlgorithm.h"

/**
 *  @brief  LArContent class
 */
class LArContent
{
public:
    #define LAR_ALGORITHM_LIST(d)                                                                                               \
        d("LArEventDisplay",                        lar_content::EventDisplayAlgorithm::Factory)                                \
        d("LArParticleAnalysis",                    lar_content::ParticleAnalysisAlgorithm::Factory)                            \
        d("LArParticleMonitoring",                  lar_content::ParticleMonitoringAlgorithm::Factory)                          \
        d("LArVisualMonitoring",                    lar_content::VisualMonitoringAlgorithm::Factory)                            \
        d("LArCheatingClusterCharacterisation",     lar_content::CheatingClusterCharacterisationAlgorithm::Factory)             \
        d("LArCheatingClusterCreation",             lar_content::CheatingClusterCreationAlgorithm::Factory)                     \
        d("LArCheatingCosmicRayIdentification",     lar_content::CheatingCosmicRayIdentificationAlg::Factory)                   \
        d("LArCheatingCosmicRayShowerMatching",     lar_content::CheatingCosmicRayShowerMatchingAlg::Factory)                   \
        d("LArCheatingPfoCreation",                 lar_content::CheatingPfoCreationAlgorithm::Factory)                         \
        d("LArCheatingVertexCreation",              lar_content::CheatingVertexCreationAlgorithm::Factory)                      \
        d("LArCosmicRayIdentification",             lar_content::CosmicRayIdentificationAlgorithm::Factory)                     \
        d("LArCosmicRayTrackMatching",              lar_content::CosmicRayTrackMatchingAlgorithm::Factory)                      \
        d("LArCosmicRayVertexBuilding",             lar_content::CosmicRayVertexBuildingAlgorithm::Factory)                     \
        d("LArNeutrinoEventBuilding",               lar_content::NeutrinoEventBuildingAlgorithm::Factory)                       \
        d("LArNeutrinoEventCreation",               lar_content::NeutrinoEventCreationAlgorithm::Factory)                       \
        d("LArNeutrinoVertexBuilding",              lar_content::NeutrinoVertexBuildingAlgorithm::Factory)                      \
        d("LArNeutrinoVertexCreation",              lar_content::NeutrinoVertexCreationAlgorithm::Factory)                      \
        d("LArDeltaRayIdentification",              lar_content::DeltaRayIdentificationAlgorithm::Factory)                      \
        d("LArDeltaRayMatching",                    lar_content::DeltaRayMatchingAlgorithm::Factory)                            \
        d("LArThreeDHitCreation",                   lar_content::ThreeDHitCreationAlgorithm::Factory)                           \
        d("LArThreeDLongitudinalTracks",            lar_content::ThreeDLongitudinalTracksAlgorithm::Factory)                    \
        d("LArParticleRecovery",                    lar_content::ParticleRecoveryAlgorithm::Factory)                            \
        d("LArVertexBasedPfoMerging",               lar_content::VertexBasedPfoMergingAlgorithm::Factory)                       \
        d("LArThreeDRemnants",                      lar_content::ThreeDRemnantsAlgorithm::Factory)                              \
        d("LArThreeDShowers",                       lar_content::ThreeDShowersAlgorithm::Factory)                               \
        d("LArThreeDTrackFragments",                lar_content::ThreeDTrackFragmentsAlgorithm::Factory)                        \
        d("LArThreeDTransverseTracks",              lar_content::ThreeDTransverseTracksAlgorithm::Factory)                      \
        d("LArLongitudinalAssociation",             lar_content::LongitudinalAssociationAlgorithm::Factory)                     \
        d("LArLongitudinalExtension",               lar_content::LongitudinalExtensionAlgorithm::Factory)                       \
        d("LArSimpleClusterMerging",                lar_content::SimpleClusterMergingAlgorithm::Factory)                        \
        d("LArTransverseAssociation",               lar_content::TransverseAssociationAlgorithm::Factory)                       \
        d("LArTransverseExtension",                 lar_content::TransverseExtensionAlgorithm::Factory)                         \
        d("LArSimpleClusterCreation",               lar_content::SimpleClusterCreationAlgorithm::Factory)                       \
        d("LArTrackClusterCreation",                lar_content::TrackClusterCreationAlgorithm::Factory)                        \
        d("LArClusteringParent",                    lar_content::ClusteringParentAlgorithm::Factory)                            \
        d("LArBoundedClusterMerging",               lar_content::BoundedClusterMergingAlgorithm::Factory)                       \
        d("LArConeBasedMerging",                    lar_content::ConeBasedMergingAlgorithm::Factory)                            \
        d("LArIsolatedHitMerging",                  lar_content::IsolatedHitMergingAlgorithm::Factory)                          \
        d("LArProximityBasedMerging",               lar_content::ProximityBasedMergingAlgorithm::Factory)                       \
        d("LArCosmicRayExtension",                  lar_content::CosmicRayExtensionAlgorithm::Factory)                          \
        d("LArCosmicRaySplitting",                  lar_content::CosmicRaySplittingAlgorithm::Factory)                          \
        d("LArDeltaRayExtension",                   lar_content::DeltaRayExtensionAlgorithm::Factory)                           \
        d("LArDeltaRayGrowing",                     lar_content::DeltaRayGrowingAlgorithm::Factory)                             \
        d("LArBranchSplitting",                     lar_content::BranchSplittingAlgorithm::Factory)                             \
        d("LArCrossedTrackSplitting",               lar_content::CrossedTrackSplittingAlgorithm::Factory)                       \
        d("LArDeltaRaySplitting",                   lar_content::DeltaRaySplittingAlgorithm::Factory)                           \
        d("LArKinkSplitting",                       lar_content::KinkSplittingAlgorithm::Factory)                               \
        d("LArLayerSplitting",                      lar_content::LayerSplittingAlgorithm::Factory)                              \
        d("LArTrackConsolidation",                  lar_content::TrackConsolidationAlgorithm::Factory)                          \
        d("LArVertexSplitting",                     lar_content::VertexSplittingAlgorithm::Factory)                             \
        d("LArClusterCharacterisation",             lar_content::ClusterCharacterisationAlgorithm::Factory)                     \
        d("LArShowerGrowing",                       lar_content::ShowerGrowingAlgorithm::Factory)                               \
        d("LArTwoDParticleCreationAlgorithm",       lar_content::TwoDParticleCreationAlgorithm::Factory)                        \
        d("LArListChanging",                        lar_content::ListChangingAlgorithm::Factory)                                \
        d("LArListDissolution",                     lar_content::ListDissolutionAlgorithm::Factory)                             \
        d("LArListMerging",                         lar_content::ListMergingAlgorithm::Factory)                                 \
        d("LArListPreparation",                     lar_content::ListPreparationAlgorithm::Factory)                             \
        d("LArCandidateVertexCreation",             lar_content::CandidateVertexCreationAlgorithm::Factory)                     \
        d("LArVertexSelection",                     lar_content::VertexSelectionAlgorithm::Factory)

    #define LAR_ALGORITHM_TOOL_LIST(d)                                                                                          \
        d("LArClearShowers",                        lar_content::ClearShowersTool::Factory)                                     \
        d("LArShowerTensorVisualization",           lar_content::ShowerTensorVisualizationTool::Factory)                        \
        d("LArSimpleShowers",                       lar_content::SimpleShowersTool::Factory)                                    \
        d("LArSplitShowers",                        lar_content::SplitShowersTool::Factory)                                     \
        d("LArClearTrackFragments",                 lar_content::ClearTrackFragmentsTool::Factory)                              \
        d("LArClearLongitudinalTrackHits",          lar_content::ClearLongitudinalTrackHitsTool::Factory)                       \
        d("LArClearTransverseTrackHits",            lar_content::ClearTransverseTrackHitsTool::Factory)                         \
        d("LArDeltaRayShowerHits",                  lar_content::DeltaRayShowerHitsTool::Factory)                               \
        d("LArMultiValuedLongitudinalTrackHits",    lar_content::MultiValuedLongitudinalTrackHitsTool::Factory)                 \
        d("LArMultiValuedTransverseTrackHits",      lar_content::MultiValuedTransverseTrackHitsTool::Factory)                   \
        d("LArThreeViewShowerHits",                 lar_content::ThreeViewShowerHitsTool::Factory)                              \
        d("LArTwoViewShowerHits",                   lar_content::TwoViewShowerHitsTool::Factory)                                \
        d("LArClearLongitudinalTracks",             lar_content::ClearLongitudinalTracksTool::Factory)                          \
        d("LArMatchedEndPoints",                    lar_content::MatchedEndPointsTool::Factory)                                 \
        d("LArClearRemnants",                       lar_content::ClearRemnantsTool::Factory)                                    \
        d("LArClearTracks",                         lar_content::ClearTracksTool::Factory)                                      \
        d("LArLongTracks",                          lar_content::LongTracksTool::Factory)                                       \
        d("LArMissingTrack",                        lar_content::MissingTrackTool::Factory)                                     \
        d("LArMissingTrackSegment",                 lar_content::MissingTrackSegmentTool::Factory)                              \
        d("LArOvershootTracks",                     lar_content::OvershootTracksTool::Factory)                                  \
        d("LArTrackSplitting",                      lar_content::TrackSplittingTool::Factory)                                   \
        d("LArTransverseTensorVisualization",       lar_content::TransverseTensorVisualizationTool::Factory)                    \
        d("LArUndershootTracks",                    lar_content::UndershootTracksTool::Factory)

    #define LAR_PARTICLE_ID_LIST(d)                                                                                             \
        d("LArMuonId",                              lar_content::LArParticleIdPlugins::LArMuonId)

    /**
     *  @brief  Register all the lar content algorithms and tools with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);

    /**
     *  @brief  Register the basic lar content plugins with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterBasicPlugins(const pandora::Pandora &pandora);

    /**
     *  @brief  Register pseudo layer plugin with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  pPseudoLayerPlugin the address of the pseudo layer plugin
     */
    static pandora::StatusCode SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, lar_content::LArPseudoLayerPlugin *pLArPseudoLayerPlugin);

    /**
     *  @brief  Register lar coordinate transformation plugin with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     *  @param  pLArTransformationPlugin the address of the lar transformation plugin
     */
    static pandora::StatusCode SetLArTransformationPlugin(const pandora::Pandora &pandora, lar_content::LArTransformationPlugin *pLArTransformationPlugin);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ALGORITHM_LIST(PANDORA_REGISTER_ALGORITHM);
    LAR_ALGORITHM_TOOL_LIST(PANDORA_REGISTER_ALGORITHM_TOOL);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::RegisterBasicPlugins(const pandora::Pandora &pandora)
{
    LAR_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArPseudoLayerPlugin(const pandora::Pandora &pandora, lar_content::LArPseudoLayerPlugin *pLArPseudoLayerPlugin)
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, pLArPseudoLayerPlugin));

    return lar_content::LArGeometryHelper::SetLArPseudoLayerPlugin(pandora, pLArPseudoLayerPlugin);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArContent::SetLArTransformationPlugin(const pandora::Pandora &pandora, lar_content::LArTransformationPlugin *pLArTransformationPlugin)
{
    return lar_content::LArGeometryHelper::SetLArTransformationPlugin(pandora, pLArTransformationPlugin);
}

#endif // #ifndef LAR_CONTENT_H
