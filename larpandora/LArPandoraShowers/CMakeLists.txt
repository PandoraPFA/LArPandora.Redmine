include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )
include_directories( $ENV{EIGEN_INC} )

art_make( 
          LIB_LIBRARIES ${PANDORASDK}
                        ${PANDORAMONITORING}
                        LArPandoraContent
                        MODULE_LIBRARIES larpandora_LArPandoraShowers
          )

install_headers()
install_fhicl()
install_source()
