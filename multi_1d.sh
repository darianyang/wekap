#!/bin/bash
# run multiple kinetics analysis pipelines

SYSTEMS=(v00 v01 v02 threshv00)
#SYSTEMS=(v00)
# $3 in kinetics pipeline
#WEST_DIR=2kod_c2_80
WEST_DIR=1a43_c2_we

for SYS in ${SYSTEMS[@]} ; do 
    for ANGLE in {54..68..1} ; do
    #for ANGLE in {36..20..1} ; do
        bash kinetics_pipeline_1d_1a43.sh $ANGLE $SYS $WEST_DIR
        #bash kinetics_pipeline_1d_2kod.sh $ANGLE $SYS $WEST_DIR
        wait
    done
done
