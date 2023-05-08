#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(WT_v00)
#SYSTEMS=(WT_v01 WT_v02 WT_v03 WT_v04 4F_v00 7F_v00)
#SYSTEMS=(WT_v01 WT_v02 WT_v03 WT_v04)
SYSTEMS=(4F_v00 7F_v00)

# $3 in kinetics pipeline
WEST_DIR=oa1_oa2_c2

for SYS in ${SYSTEMS[@]} ; do 
    for C2 in {72..72} ; do
    # o_angle : 52, 55, 60
    for ANGLE in {55..55} ; do
        bash kinetics_pipeline_3d.sh $ANGLE $SYS $WEST_DIR $C2
        wait
    done
    done
done
