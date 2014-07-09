cet_add_compiler_flags( CXX -DMONITORING=1 )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{PANDORA_FQ_DIR}/PandoraSDK/include )
include_directories( $ENV{PANDORA_FQ_DIR}/PandoraMonitoring/include )

art_make( BASENAME_ONLY
          LIBRARY_NAME  LArPandoraAlgorithms
          SUBDIRS LArCalculators
                  LArCheating
                  LArHelpers
                  LArMonitoring
                  LArObjects
                  LArThreeDReco
                  LArThreeDReco/LArClusterSplitting  
                  LArThreeDReco/LArCosmicRay  
                  LArThreeDReco/LArHitCreation
                  LArThreeDReco/LArShowerMatching
                  LArThreeDReco/LArShowerFragments
                  LArThreeDReco/LArLongitudinalTrackMatching
                  LArThreeDReco/LArTransverseTrackMatching
                  LArThreeDReco/LArTrackFragments
                  LArTwoDReco
                  LArTwoDReco/LArClusterAssociation
                  LArTwoDReco/LArClusterCreation
                  LArTwoDReco/LArClusterMopUp
                  LArTwoDReco/LArClusterSplitting
                  LArTwoDReco/LArCosmicRay
                  LArTwoDReco/LArSeedFinding
                  LArUtility
                  LArVertex
          LIB_LIBRARIES ${PANDORASDK}
                        ${PANDORAMONITORING}
          )

install_headers( SUBDIRS LArCalculators
                         LArCheating
                         LArHelpers
                         LArMonitoring
                         LArObjects
                         LArThreeDReco
                         LArThreeDReco/LArClusterSplitting
                         LArThreeDReco/LArCosmicRay  
                         LArThreeDReco/LArHitCreation
                         LArThreeDReco/LArShowerMatching
                         LArThreeDReco/LArShowerFragments
                         LArThreeDReco/LArLongitudinalTrackMatching
                         LArThreeDReco/LArTransverseTrackMatching
                         LArThreeDReco/LArTrackFragments
                         LArThreeDReco/LArThreeDBase
                         LArTwoDReco
                         LArTwoDReco/LArClusterAssociation
                         LArTwoDReco/LArClusterCreation
                         LArTwoDReco/LArClusterMopUp
                         LArTwoDReco/LArClusterSplitting
                         LArTwoDReco/LArCosmicRay
                         LArTwoDReco/LArSeedFinding
                         LArUtility
                         LArVertex )
install_fhicl( SUBDIRS LArCalculators
                       LArCheating
                       LArHelpers
                       LArMonitoring
                       LArObjects
                       LArThreeDReco
                       LArThreeDReco/LArClusterSplitting
                       LArThreeDReco/LArCosmicRay  
                       LArThreeDReco/LArHitCreation
                       LArThreeDReco/LArShowerMatching
                       LArThreeDReco/LArShowerFragments
                       LArThreeDReco/LArLongitudinalTrackMatching
                       LArThreeDReco/LArTransverseTrackMatching
                       LArThreeDReco/LArTrackFragments
                       LArThreeDReco/LArThreeDBase
                       LArTwoDReco
                       LArTwoDReco/LArClusterAssociation
                       LArTwoDReco/LArClusterCreation
                       LArTwoDReco/LArClusterMopUp
                       LArTwoDReco/LArClusterSplitting
                       LArTwoDReco/LArCosmicRay
                       LArTwoDReco/LArSeedFinding
                       LArUtility
                       LArVertex )
install_source( SUBDIRS LArCalculators
                        LArCheating
                        LArHelpers
                        LArMonitoring
                        LArObjects
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
                        LArUtility
                        LArVertex  )

