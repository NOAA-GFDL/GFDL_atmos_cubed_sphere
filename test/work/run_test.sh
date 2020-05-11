#!/bin/bash

set -x
set -e

cd /work

if [ -d ./rundir ] ; then
    echo "ERROR: rundir already exists"
    exit 1
fi

#./create_rundir.py
curl -X GET 'https://www.googleapis.com/storage/v1/b/vcm-ml-public/o/data_for_demos%2FTravisCI%2Frundir.tar.gz?alt=media' | tar xvz

if [ ! -d ./rundir -o ! -f ./rundir/INPUT/gfs_data.tile1.nc ] ; then
    echo "ERROR: problem creating rundir"
    exit 1
fi

cd /work/rundir
ln -s /FV3/fv3.exe .
mpirun --allow-run-as-root --mca btl_vader_single_copy_mechanism none \
       --oversubscribe --merge-stderr-to-stdout --tag-output --timestamp-output \
       ./fv3.exe 2>&1 | tee fv3.out

exit 0
