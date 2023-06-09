#!/bin/bash
# kinetics_pipeline.sh
# Reproduces the results of w_ipa -ao, but manually with self defined aux values
# w_ipa is essentially: w_assign and w_direct
# w_direct all is: w_kinetics, w_kinavg, and w_stateprobs
# Pipeline of commands to get kinetics information from west.h5 file.
# Run in a directory with target west.h5 file
# currently configured to run kinetics for a 2D scheme with 2 aux values

# 2kod v02 crawled file parent dir
#WDIR=/Users/darian/Google/MBSB/Research/Projects/hiv1_capsid/ctd_1d_we_sim/kinetics
#WEST=west_i200_crawled.h5
#WDIR=/Users/darian/bigfiles/multi_west
#WEST=multi_2kod.h5
WDIR=/home/dty7/we_data
WEST=v00.h5
AUX_A="pcoord"
#AUX_B="c2_angle"
#SCHEME="o_and_c2_angle_2kod_v00" 

TS_DEF=$1
SYSTEM=$2
WEST_DIR=$3

WEST=$WEST_DIR/${SYSTEM}.h5

SCHEME="$WEST_DIR/1d_${TS_DEF}_${SYSTEM}" 

# make scheme dir and fill with current pipeline run file
mkdir -pv $SCHEME
cp -v kinetics_pipeline_1d_oa.sh $SCHEME
cd $SCHEME

# define bins and states with yaml files
cat << EOF > BINS
---
bins:
    type: RectilinearBinMapper
    # add buffer region: strict starting state definition
    boundaries: [[0.0, 20.0, 44.0, ${TS_DEF}, 65.0, 'inf']] 
EOF
cat << EOF > STATES
---
states:
  - label: a
    coords:
      - [25]

  - label: b
    coords:
      - [64]
EOF

# create module.py file to process 1D or 2D scheme
cat << EOF > module.py
#!/usr/bin/env python

import numpy

def pull_data_1d(n_iter, iter_group):
    '''
    This function reshapes auxiliary data for each iteration and returns it.
    '''
    #auxdata = iter_group['auxdata']['${AUX_A}'][...]
    auxdata = iter_group['${AUX_A}'][...]
    #data = auxdata[:,:,numpy.newaxis]
    data = numpy.atleast_3d(auxdata)
    return data

def pull_data_2d(n_iter, iter_group):
    '''
    This function reshapes 2 auxiliary datasets for each iteration and returns it.
    '''
    #auxdata1 = iter_group['auxdata']['${AUX_A}'][...]
    # special for pcoord
    auxdata1 = iter_group['${AUX_A}'][...]
    auxdata2 = iter_group['auxdata']['${AUX_B}'][...]
    data = numpy.dstack((auxdata1, auxdata2))
    return data
EOF

# run w_assign to assign macrostates based off of defined BINS and STATES
w_assign -W $WDIR/$WEST --bins-from-file BINS --states-from-file STATES -o assign.h5 --construct-dataset module.pull_data_1d --serial &&
echo "    w_assign finished: assign.h5 file created in $SCHEME directory"

# run w_direct to then calculate multiple flux values and rate constants
# 'all' is equivalent to w_kinetics, w_kinavg, and w_stateprobs
# step size of 1 is used, and can later be window averaged based off of ACF decay lag time, or can use default
# w_direct by default will set 'correl' = True, which calculates the ACF lag time and window averages
#w_direct kinetics -a assign.h5 -W $WDIR/$WEST -o direct.h5 --evolution-mode cumulative --step-iter 1 &&
w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 --evolution-mode cumulative --step-iter 1 &&
#w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 --step-iter 1 &&
#w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 &&
echo "    w_direct finished: direct.h5 file created in $SCHEME directory"

