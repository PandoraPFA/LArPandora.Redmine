
# included for Mac OSX using XQuartz for X11
if(APPLE)
  include_directories( /opt/X11/include )
endif(APPLE)

# and this is because of an error that is triggered by a X11 header
# (that seems to be impossible to disable with #pragma's)
cet_set_compiler_flags(EXTRA_CXX_FLAGS -Wno-variadic-macros)


art_make_library( LIBRARY_NAME larpandora_RootGraphicsEnabler
                  SOURCE RootGraphicsEnabler.cxx 
                  LIBRARIES
                    ${ROOT_CORE}
                    ${ROOT_GX11}
                    ${ROOT_GUI}
                )
art_make(
          EXCLUDE RootGraphicsEnabler.cxx 
          SERVICE_LIBRARIES
            larpandora_RootGraphicsEnabler
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_SERVICES_REGISTRY} 
        )

install_headers()
install_fhicl()
install_source()
