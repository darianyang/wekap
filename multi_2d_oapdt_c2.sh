#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(2kod_v00 lo_pH_v00 150end_v00)
#SYSTEMS=(2kod_v00)
#SYSTEMS=(150end_v00)
SYSTEMS=(4F_v00 7F_v00)

# $3 in kinetics pipeline
WEST_DIR=oapdt_c2_2dgrid

for SYS in ${SYSTEMS[@]} ; do 
    #for C2 in {66..72} ; do
    for C2 in {70..70} ; do
        #WEST_DIR=2kod_oa_${C2}c2
    # o_angle
    #for ANGLE in {2100..2600..100} ; do
    for ANGLE in {2500..2500} ; do
    #for ANGLE in {36..20..1} ; do
        bash kinetics_pipeline_2d_oapdt.sh $ANGLE $SYS $WEST_DIR $C2
        wait
    done
    done
done
