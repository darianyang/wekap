#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(2kod_v00 lo_pH_v00 150end_v00)
#SYSTEMS=(2kod_v00)
#SYSTEMS=(150end_v00)
SYSTEMS=(WT_v00)

# $3 in kinetics pipeline
WEST_DIR=oamax_c2_2dgrid_2

for SYS in ${SYSTEMS[@]} ; do 
    #for C2 in {66..72} ; do
    for C2 in {72..80} ; do
        #WEST_DIR=2kod_oa_${C2}c2
    # o_angle
    for ANGLE in {48..54} ; do
        bash kinetics_pipeline_2d_oamax.sh $ANGLE $SYS $WEST_DIR $C2
        wait
    done
    done
done
