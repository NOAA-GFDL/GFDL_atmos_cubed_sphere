#!/bin/sh -xe

##############################################################################
## User set up variables
## Root directory for CI
dirRoot=/contrib/fv3
## Intel version to be used
intelVersion=2023.2.0
##############################################################################
## HPC-ME container
container=/contrib/containers/noaa-intel-prototype_2023.09.25.sif
container_env_script=/contrib/containers/load_spack_noaa-intel.sh
##############################################################################

#Parse Arguments
branch=main
commit=none
while [[ $# -gt 0 ]]; do
  case $1 in
    -b|--branch)
      branch="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--hash)
      commit="$2"
      shift # past argument
      shift # past value
      ;;
    -c|--config)
      config="$2"
      shift # past argument
      shift # past value
      ;;
    --hydro)
      hydro="$2"
      shift # past argument
      shift # past value
      ;;
    --bit)
      bit="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--mode)
      mode="$2"
      shift # past argument
      shift # past value
      ;;
    *)
      echo "unknown argument"
      exit 1
      ;;
  esac
done

if [ -z $mode ] || [ -z $bit ] || [ -z $hydro ] || [ -z $config ]
  then
    echo "must specify config, hydro, bit, and mode options for compile"
    exit 1
fi

echo "branch is $branch"
echo "commit is $commit"
echo "mode is $mode"
echo "bit is $bit"
echo "hydro is $hydro"
echo "config is $config"

if [ $hydro = "sw" ] && [ $config = "shield" ]
  then
    echo "this combination should not be tested"
  else
    ## Set up the directories
    testDir=${dirRoot}/${intelVersion}/GFDL_atmos_cubed_sphere/${branch}/${commit}
    logDir=${testDir}/log
    # Set up build
    cd ${testDir}/SHiELD_build/Build
    #Define External Libs path
    export EXTERNAL_LIBS=${dirRoot}/externallibs
    # Build SHiELD
    set -o pipefail
    singularity exec -B /contrib ${container} ${container_env_script} "./COMPILE ${config} ${hydro} ${bit} ${mode} intel clean" |& tee ${logDir}/compile_${config}_${hydro}_${bit}_${mode}_intel.out
fi
