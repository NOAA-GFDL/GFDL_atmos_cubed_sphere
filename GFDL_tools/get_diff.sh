#!/bin/bash

rdir='/ncrc/home1/Mikyung.Lee/awg/xanadu/cm4p12_xanadu00/src/atmos_cubed_sphere/GFDL_tools'

ofile='ALLDIFF_JAN03'
if [ -f $ofile ] ; then 
    rm $ofile 
fi
touch -a $ofile


for myfile in *.*90 ; do
    echo "-----------------------$myfile-----------------------" >> $ofile
    diff -w $myfile $rdir/$myfile >> $ofile
    echo $myfile
done



    
