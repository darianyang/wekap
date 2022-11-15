#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(v00 v01 v02 threshv00)
SYSTEMS=(v00 v01)
# $3 in kinetics pipeline
#WEST_DIR=2kod_c2_80
WEST_DIR=2kod_oa

for SYS in ${SYSTEMS[@]} ; do 
    for ANGLE in {40..52..1} ; do
    #for ANGLE in {36..20..1} ; do
        bash kinetics_pipeline_2d_oa.sh $ANGLE $SYS $WEST_DIR
        wait
    done
done
