#!/bin/bash
# run multiple kinetics analysis pipelines

SYSTEMS=(WT_v00)

# $3 in kinetics pipeline
WEST_DIR=oa1_oa2_c2

for SYS in ${SYSTEMS[@]} ; do 
    for C2 in {72..72} ; do
    # o_angle
    for ANGLE in {60..60} ; do
        bash kinetics_pipeline_3d.sh $ANGLE $SYS $WEST_DIR $C2
        wait
    done
    done
done
