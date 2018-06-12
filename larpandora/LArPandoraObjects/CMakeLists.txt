include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{PANDORA_INC} )

art_make( BASENAME_ONLY
          LIBRARY_NAME larpandoraobj
          LIB_LIBRARIES
            ${PANDORASDK}
            cetlib cetlib_except
          DICT_LIBRARIES
            lardataobj_RecoBase
          )

install_headers()
install_fhicl()
install_source()
