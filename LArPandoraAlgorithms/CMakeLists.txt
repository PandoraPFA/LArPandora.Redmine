
cet_add_compiler_flags( CXX -DMONITORING=1 )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{PANDORA_FQ_DIR}/PandoraSDK/include )
include_directories( $ENV{PANDORA_FQ_DIR}/PandoraMonitoring/include )

art_make( BASENAME_ONLY
          LIBRARY_NAME  LArPandoraAlgorithms
          SUBDIRS LArCalculators
                  LArClusterAssociation
                  LArClusterSeedAssociation
                  LArClusterSplitting
                  LArClustering
                  LArHelpers
                  LArMonitoring
                  LArObjects
                  LArReclustering
                  LArThreeDSeed
                  LArTwoDSeed
                  LArUtility
                  LArVertex
          LIB_LIBRARIES ${PANDORASDK}
                        ${PANDORAMONITORING}
          )

install_headers( SUBDIRS LArCalculators
                  LArClusterAssociation
                  LArClusterSeedAssociation
                  LArClusterSplitting
                  LArClustering
                  LArHelpers
                  LArMonitoring
                  LArObjects
                  LArReclustering
                  LArThreeDSeed
                  LArTwoDSeed
                  LArUtility
                  LArVertex )
install_fhicl( SUBDIRS LArCalculators
                  LArClusterAssociation
                  LArClusterSeedAssociation
                  LArClusterSplitting
                  LArClustering
                  LArHelpers
                  LArMonitoring
                  LArObjects
                  LArReclustering
                  LArThreeDSeed
                  LArTwoDSeed
                  LArUtility
                  LArVertex )
install_source( SUBDIRS LArCalculators
                  LArClusterAssociation
                  LArClusterSeedAssociation
                  LArClusterSplitting
                  LArClustering
                  LArHelpers
                  LArMonitoring
                  LArObjects
                  LArReclustering
                  LArThreeDSeed
                  LArTwoDSeed
                  LArUtility
                  LArVertex )


