#!/bin/bash
#
# Author: echurch@fnal.gov from dbox@fnal.gov
# .fcl-specific ART codework by joshua.spitz@yale.edu
# A script to run argoneut framework (ART) jobs on the local cluster/grid.
#
# It takes 4 arguments, the latter two of which are specified. First two are inherited.:
#   cluster    - the condor job cluster number
#   process    - the condor job process number, within the cluster
#   user       - the username of the person who submitted the job
#   submitdir  - the directory from which the job was submitted ( not used at this time).
#
# Outputs:
#  - All output files are created in the grid scratch space.  At the end of the job
#    all files in this directory will be copied to:
#      /grid/data/argoneut/outstage/$user/${cluster}_${process}
#    This includes a copy of the input files.
#
# Notes:
#
# 1) For documentation on using the grid, see
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtm
#    For details on cpn and outstage see:
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml#cpn
#      http://mu2e.fnal.gov/atwork/computing/fermigrid.shtml#outstage
#
sleep `echo $((RANDOM%100+0))`
# sleep 10000
verbose=T
source /grid/fermiapp/minos/minossoft/setup/setup_minossoft_FNALU.sh R2.2
# Copy arguments into meaningful names.
cluster=${CLUSTER}
process=${PROCESS}
user=$1
submitdir=$2
echo "Input arguments:"
echo "Cluster:    " $cluster
echo "Process:    " $process
echo "User:       " $user
echo "SubmitDir:  " $submitdir
echo " "

# Do not change this section.
# It creates a temporary working directory that automatically cleans up all
# leftover files at the end.
# Next 3 commented-out lines are outdated cruft for grid, now replaced at 
# Dennis's suggestion. EC, 23-Nov-2010. 
ORIGDIR=`pwd`
#TMP=`mktemp -d ${OSG_WN_TMP:-/var/tmp}/working_dir.XXXXXXXXXX`
#TMP=${TMP:-${OSG_WN_TMP:-/var/tmp}/working_dir.$$}
TMP=`mktemp -d ${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${_CONDOR_SCRATCH_DIR:-/var/tmp}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { cd ; rm -rf \"$TMP\"; }" 0


outstage=/argoneut/data/outstage/$user

j=`echo $RANDOM`


cd $TMP
# cp -r $ORIGDIR/* .


for i in `ls -d /argoneut/data/outstage/saima/878844_*|sort`; 
do 

if [ -f $i/loon -o -f $i/MINOSin.sntp.root ]
then
echo $i exists
else
echo Working on $i;
touch "$i"/loon
/grid/fermiapp/minos/sim/minerva/process_minerva_mc.sh -i "$i"/MINOSin.txt -s /argoneut/data/outstage/saima/scratch -o "$i" -x 0,0,0,0 photonnoisewindow=2.0e-6 

fi
done


for i in `ls -d /argoneut/data/outstage/saima/878844_*|sort`; 
do 

if [ -f $i/loon2 -o -f $i/MINOSout_fhc.root ]
then
echo $i/MINOSout_fhc.root exists
else
problem=`grep -r -i sighup "$i"/MINOSin.detsim.log`
if [ -z "$problem" ]
then
if [  -f $i/MINOSin.sntp.root ]
then
touch "$i"/loon2
loon -bq /argoneut/app/users/saima/myloadlibs.C '/argoneut/app/users/saima/microDst.C+("'$i'/MINOSin.sntp.root","'$i'/MINOSout_fhc.root")'
fi
fi

fi

done



unset LD_LIBRARY_PATH

# Establish environment and run the job.
#source /grid/fermiapp/products/mu2e/setupmu2e.sh
#source /grid/fermiapp/mu2e/Offline/v0_2_6/setup.sh
export GROUP=argoneut
export EXPERIMENT=argoneut
export HOME=/argoneut/app/users/saima/svn_commit/sam_batch
## This sets all the needed FW and SRT and LD_LIBRARY_PATH envt variables. 
## Then we cd back to our TMP area. EC, 23-Nov-2010.
source /grid/fermiapp/lbne/lar/code/larsoft/releases/development/setup/setup_larsoft_fnal.sh
cd /argoneut/app/users/saima/svn_commit; srt_setup -a



for i in `ls -d /argoneut/data/outstage/saima/878844_*|sort `; 
do 

if [ -f $i/mergesim.root -o -f $i/loon3 ]
then
continue
fi

if [ -f $i/MINOSout_fhc.root ]
then
sleep `echo $((RANDOM%10+0))`
export DIR="$i"
echo $DIR
cd $TMP
touch "$i"/loon3
lar -c job/mergemc.fcl -s "$i"/t962gana.root >&thisjob3.log


chmod -R g+rw *
stat=`echo $?`

if [ "$stat" == "0" ]
then
/grid/fermiapp/minos/scripts/cpn * "$i"
fi


fi

done




# fw thisjob.py >& thisjob.log


# echo $i
# echo $i>>/argoneut/app/users/spitz7/scanminoslist
# Make sure the user's output staging area exists.
# test -e $outstage || mkdir $outstage
# if [ ! -d $outstage ];then
#    echo "File exists but is not a directory."
#    outstage=/grid/data/argoneut/outstage/nobody
#    echo "Changing outstage directory to: " $outstage 
#    exit
# fi
# 
# # Make a directory in the outstage area to hold all files from this job.
# mkdir ${outstage}/${cluster}_${process}
# rm /argoneut/app/users/spitz7/larsoft_ART2/batch/minos_scan_merge_"$j".fcl
# Copy all files from the working directory to the output staging area.
# /grid/fermiapp/minos/scripts/cpn * ${outstage}/${cluster}_${process}

exit 0