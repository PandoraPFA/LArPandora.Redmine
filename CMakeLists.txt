# ======================================================================
#  larpandora main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larpandora/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

project(larpandora)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
		      $ENV{CETBUILDTOOLS_DIR}/Modules
		      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
  )

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larcoreobj )
find_ups_product( larcorealg )
find_ups_product( larcore v1_00_00 )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( lardata v1_00_00 )
find_ups_product( larevt v1_00_00 )
find_ups_product( larsim v1_00_00 )
find_ups_product( larreco v1_00_00 )
find_ups_product( nutools v0_00_00 )
find_ups_product( art v1_09_00 )
find_ups_product( pandora v02_00_00 )
find_ups_product( larpandoracontent v02_04_00 )
find_ups_product( cetbuildtools v4_09_00 )
find_ups_product( postgresql v9_1_5 )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

cet_find_library( PANDORASDK NAMES PandoraSDK PATHS ENV PANDORA_LIB )
cet_find_library( PANDORAMONITORING NAMES PandoraMonitoring PATHS ENV PANDORA_LIB )

# find larpandoracontent headers if building at the same time
#message(STATUS "larpandora: checking for MRB_SOURCE")
set( mrb_source $ENV{MRB_SOURCE} )
if( LARPANDORACONTENT_INC )
  #message(STATUS "larpandora: found LARPANDORACONTENT_INC: $ENV{LARPANDORACONTENT_INC}")
  # note that mrb v1_07_00 or later will define LARPANDORACONTENT_INC
  if( mrb_source )
    include_directories( ${LARPANDORACONTENT_INC} )
  endif()
else()
  if( mrb_source )
    include_directories( ${mrb_source}/larpandoracontent )
  else()
    message(FATAL_ERROR "LARPANDORACONTENT_INC is undefined, please setup larpandoracontent")
  endif()
endif()

# source
add_subdirectory(larpandora)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)

