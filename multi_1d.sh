#!/bin/bash
# run multiple kinetics analysis pipelines

for ANGLE in {46..56..2} ; do
    bash kinetics_pipeline_1d.sh $ANGLE
    wait
done
