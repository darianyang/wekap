#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(v00 v01 v02 threshv00)
#SYSTEMS=(2kod_v00 lo_pH_v00 150end_v00)
SYSTEMS=(2kod_v00 lo_pH_v00)
#SYSTEMS=(150end_v00)
# $3 in kinetics pipeline
WEST_DIR=gt2500oapdt_gt65c2

for SYS in ${SYSTEMS[@]} ; do 
    for C2 in {65..65} ; do
        #WEST_DIR=2kod_oa_${C2}c2
    # o_angle
    for ANGLE in {2500..2500} ; do
    #for ANGLE in {36..20..1} ; do
        bash kinetics_pipeline_2d_oapdt.sh $ANGLE $SYS $WEST_DIR $C2
        wait
    done
    done
done
