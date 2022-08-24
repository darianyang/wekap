#!/bin/bash
# run multiple kinetics analysis pipelines

SYSTEMS=(v00 v01 v02 threshv00)

for SYS in ${SYSTEMS[@]} ; do 
    for ANGLE in {56..68..1} ; do
        bash kinetics_pipeline_1d.sh $ANGLE $SYS
        wait
    done
done
