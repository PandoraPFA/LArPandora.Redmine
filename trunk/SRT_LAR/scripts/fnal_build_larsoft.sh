#!/bin/sh

# Wrapper Script to build the LArSoft software.
# This sources the development release of the software
# and updates all of the packages.
# It then calls a perl script to make a clean build 
# in standard and debug modes.
#
# February 13 2011 - brebel@fnal.gov
#

source /grid/fermiapp/lbne/lar/code/larsoft/setup/setup_larsoft_fnal.sh
cd ${SRT_DIST}/setup
svn update

UPDATE_LOG=${LARHOME}/lar_update.log

echo "Update FNAL LArSoft at " >& ${UPDATE_LOG}
date >> ${UPDATE_LOG}

${SRT_PUBLIC_CONTEXT}/SRT_LAR/scripts/lar_update_rel -rel development >> ${UPDATE_LOG}

conflict=`grep -i 'C ' ${UPDATE_LOG}`

if [ "${conflict}" -ne "1" ]; then
   echo "Possible merge conflicts detected.  Will not build" >> ${UPDATE_LOG}
fi

${SRT_PUBLIC_CONTEXT}/SRT_LAR/scripts/lar_build -rel development -log $LARHOME/lar_make.log
. srt_environment_lar -X -s SRT_QUAL=debug
${SRT_PUBLIC_CONTEXT}/SRT_LAR/scripts/lar_build -rel development -debug -log $LARHOME/lar_debug_make.log
