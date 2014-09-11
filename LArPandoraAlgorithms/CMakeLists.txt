cet_add_compiler_flags( CXX -DMONITORING=1 )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{PANDORA_FQ_DIR}/PandoraSDK/include )
include_directories( $ENV{PANDORA_FQ_DIR}/PandoraMonitoring/include )

art_make( BASENAME_ONLY
          LIBRARY_NAME  LArPandoraAlgorithms
          SUBDIRS LArCheating
                  LArHelpers
                  LArMonitoring
                  LArObjects
                  LArPlugins
                  LArThreeDReco
                  LArThreeDReco/LArClusterSplitting  
                  LArThreeDReco/LArCosmicRay  
                  LArThreeDReco/LArEventBuilding
                  LArThreeDReco/LArHitCreation
                  LArThreeDReco/LArLongitudinalTrackMatching
                  LArThreeDReco/LArShowerFragments
                  LArThreeDReco/LArShowerMatching
                  LArThreeDReco/LArThreeDBase
                  LArThreeDReco/LArTrackFragments
                  LArThreeDReco/LArTransverseTrackMatching
                  LArTwoDReco
                  LArTwoDReco/LArClusterAssociation
                  LArTwoDReco/LArClusterCreation
                  LArTwoDReco/LArClusterMopUp
                  LArTwoDReco/LArClusterSplitting
                  LArTwoDReco/LArCosmicRay
                  LArTwoDReco/LArSeedFinding
                  LArUtility
          LIB_LIBRARIES ${PANDORASDK}
                        ${PANDORAMONITORING}
          )

install_headers( SUBDIRS LArCheating
                         LArHelpers
                         LArMonitoring
                         LArObjects
                         LArPlugins
                         LArThreeDReco
                         LArThreeDReco/LArClusterSplitting
                         LArThreeDReco/LArCosmicRay  
                         LArThreeDReco/LArEventBuilding
                         LArThreeDReco/LArHitCreation
                         LArThreeDReco/LArLongitudinalTrackMatching
                         LArThreeDReco/LArShowerFragments
                         LArThreeDReco/LArShowerMatching
                         LArThreeDReco/LArThreeDBase
                         LArThreeDReco/LArTrackFragments
                         LArThreeDReco/LArTransverseTrackMatching
                         LArTwoDReco
                         LArTwoDReco/LArClusterAssociation
                         LArTwoDReco/LArClusterCreation
                         LArTwoDReco/LArClusterMopUp
                         LArTwoDReco/LArClusterSplitting
                         LArTwoDReco/LArCosmicRay
                         LArTwoDReco/LArSeedFinding
                         LArUtility )
install_fhicl( SUBDIRS LArCheating
                       LArHelpers
                       LArMonitoring
                       LArObjects
                       LArPlugins
                       LArThreeDReco
                       LArThreeDReco/LArClusterSplitting
                       LArThreeDReco/LArCosmicRay  
                       LArThreeDReco/LArEventBuilding
                       LArThreeDReco/LArHitCreation
                       LArThreeDReco/LArLongitudinalTrackMatching
                       LArThreeDReco/LArShowerFragments
                       LArThreeDReco/LArShowerMatching
                       LArThreeDReco/LArThreeDBase
                       LArThreeDReco/LArTrackFragments
                       LArThreeDReco/LArTransverseTrackMatching
                       LArTwoDReco
                       LArTwoDReco/LArClusterAssociation
                       LArTwoDReco/LArClusterCreation
                       LArTwoDReco/LArClusterMopUp
                       LArTwoDReco/LArClusterSplitting
                       LArTwoDReco/LArCosmicRay
                       LArTwoDReco/LArSeedFinding
                       LArUtility )
install_source( SUBDIRS LArCheating
                        LArHelpers
                        LArMonitoring
                        LArObjects
                        LArPlugins
                        LArThreeDReco
                        LArThreeDReco/LArClusterSplitting
                        LArThreeDReco/LArCosmicRay  
                        LArThreeDReco/LArHitCreation
                        LArThreeDReco/LArShowerMatching
                        LArThreeDReco/LArShowerFragments
                        LArThreeDReco/LArTransverseTrackMatching
                        LArThreeDReco/LArLongitudinalTrackMatching
                        LArThreeDReco/LArTrackFragments
                        LArThreeDReco/LArThreeDBase
                        LArTwoDReco
                        LArTwoDReco/LArClusterAssociation
                        LArTwoDReco/LArClusterCreation
                        LArTwoDReco/LArClusterMopUp
                        LArTwoDReco/LArClusterSplitting
                        LArTwoDReco/LArCosmicRay
                        LArTwoDReco/LArSeedFinding
                        LArUtility )

