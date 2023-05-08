#!/bin/bash
# kinetics_pipeline.sh
# Reproduces the results of w_ipa -ao, but manually with self defined aux values
# w_ipa is essentially: w_assign and w_direct
# w_direct all is: w_kinetics, w_kinavg, and w_stateprobs
# Pipeline of commands to get kinetics information from west.h5 file.
# Run in a directory with target west.h5 file
# currently configured to run kinetics for a 2D scheme with 2 aux values

AUX_A="o_angle_m1"
AUX_B="o_angle_m2"
AUX_C="c2_angle"

TS_DEF=$1
SYSTEM=$2
WEST_DIR=$3

C2=$4

WDIR=~/we_data/1d_oamax
WEST=${SYSTEM}.h5

SCHEME="$WEST_DIR/${SYSTEM}/${TS_DEF}oamax_${C2}c2" 
#SCHEME="$WEST_DIR" 

# make scheme dir and fill with current pipeline run file
mkdir -pv $SCHEME
cp -v kinetics_pipeline_3d.sh $SCHEME
cd $SCHEME

# define bins and states with yaml files
cat << EOF > BINS
---
bins:
    type: RectilinearBinMapper
    # add buffer region: strict starting state definition
    boundaries: [[0.0, 50.0, ${TS_DEF}, 'inf'], 
                 [0.0, 50.0, ${TS_DEF}, 'inf'],
                 [0.0, 65.0, ${C2}, 'inf'],
                ]
EOF
cat << EOF > STATES
---
states:
  - label: a
    coords:
      - [25, 25, 50]

  - label: b
    coords:
      - [65, 25, 100]

  - label: c
    coords:
      - [25, 65, 100]

  - label: bc
    coords:
      - [65, 25, 100]
      - [25, 65, 100]
EOF

# create module.py file to process 1D or 2D scheme
cat << EOF > module.py
#!/usr/bin/env python

import numpy

def pull_data_1d(n_iter, iter_group):
    '''
    This function reshapes auxiliary data for each iteration and returns it.
    '''
    auxdata = iter_group['auxdata']['${AUX_A}'][...]
    data = auxdata[:,:,numpy.newaxis]
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

def pull_data_3d(n_iter, iter_group):
    '''
    This function reshapes 3 auxiliary datasets for each iteration and returns it.
    '''
    # special for pcoord
    #auxdata1 = iter_group['${AUX_A}'][...]

    auxdata1 = iter_group['auxdata']['${AUX_A}'][...]
    auxdata2 = iter_group['auxdata']['${AUX_B}'][...]
    auxdata3 = iter_group['auxdata']['${AUX_C}'][...]
    data = numpy.dstack((auxdata1, auxdata2, auxdata3))
    return data
EOF

# run w_assign to assign macrostates based off of defined BINS and STATES
w_assign -W $WDIR/$WEST --bins-from-file BINS --states-from-file STATES -o assign.h5 --construct-dataset module.pull_data_3d --serial &&
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

