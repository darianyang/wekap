#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(v00 v01 v02 threshv00)
SYSTEMS=(v00 v01 v02)
# $3 in kinetics pipeline
#WEST_DIR=2kod_c2_80
#WEST_DIR=2kod_oa_65c2
#WEST_DIR=2kod_oa_60c2

for SYS in ${SYSTEMS[@]} ; do 
    for C2 in {65..75} ; do
        WEST_DIR=2kod_oa_${C2}c2
    # o_angle
    for ANGLE in {45..55..1} ; do
    #for ANGLE in {36..20..1} ; do
        bash kinetics_pipeline_2d_oa.sh $ANGLE $SYS $WEST_DIR $C2
        wait
    done
    done
done
