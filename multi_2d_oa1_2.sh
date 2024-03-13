#!/bin/bash
# run multiple kinetics analysis pipelines

#SYSTEMS=(WT_v00)
#SYSTEMS=(WT_v00 WT_v01 WT_v02 WT_v03 WT_v04 4F_v00 7F_v00)
#SYSTEMS=(WT_v00 WT_v01 WT_v02 WT_v03 WT_v04)
#SYSTEMS=(4F_v00 7F_v00)
#SYSTEMS=(4F_v01 4F_v02 4F_v03 4F_v04 7F_v01 7F_v02 7F_v03 7F_v04)
#SYSTEMS=(WT_v00 WT_v01 WT_v02 WT_v03 WT_v04 4F_v00 4F_v01 4F_v02 4F_v03 4F_v04 7F_v00 7F_v01 7F_v02 7F_v03 7F_v04)
#SYSTEMS=(3d_v00 3d_v01)
SYSTEMS=(7F_v03)

for SYS in ${SYSTEMS[@]} ; do 
    # o_angle 
    #for ANGLE in {12..16..2} ; do
    for ANGLE in {14..14} ; do
        bash kinetics_pipeline_2d_oa1_2.sh $ANGLE $SYS 
        wait
    done
done
