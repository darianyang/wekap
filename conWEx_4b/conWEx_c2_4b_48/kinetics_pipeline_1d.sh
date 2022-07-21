#!/bin/bash
# kinetics_pipeline.sh
# Reproduces the results of w_ipa -ao, but manually with self defined aux values
# w_ipa is essentially: w_assign and w_direct
# w_direct all is: w_kinetics, w_kinavg, and w_stateprobs
# Pipeline of commands to get kinetics information from west.h5 file.
# Run in a directory with target west.h5 file

# optional arg1 to set state bounds

# 2kod v02 crawled file parent dir
#WDIR=/Users/darian/Google/MBSB/Research/Projects/hiv1_capsid/ctd_1d_we_sim/kinetics
#WEST=west_i200_crawled.h5
#WDIR=/Users/darian/bigfiles/multi_west
#WEST=multi_2kod.h5
#WDIR=/Users/darian/github/wedap/wedap/data
WDIR=/home/dty7/we_data
#WEST=west_c2.h5
#WEST=west_c2x.h5
WEST=west_c2x_4b.h5
AUX_A="1_75_39_c2"
#AUX_A="rms_bb_xtal" # > 3.4, 3.5, 3.6
SCHEME="conWEx_4b/conWEx_c2_4b_$1" 
#SCHEME="v02_2kod_lt32C2" 

# make scheme dir and fill with current pipeline run file
mkdir -pv $SCHEME
cp -v kinetics_pipeline_1d.sh $SCHEME
cd $SCHEME

# define bins and states with yaml files
cat << EOF > BINS
---
bins:
    type: RectilinearBinMapper
    # add buffer region: strict starting state definition
    boundaries: [[0.0, 30.0, 40.0, $1.0, 66.0, 'inf']]
    #boundaries: [[0.0, 2.8, 3.6, 5.0, 'inf']]
EOF
cat << EOF > STATES
---
states:
  - label: a
    coords:
      - [35]

  - label: b
    coords:
      - [63]
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
    auxdata1 = iter_group['auxdata']['${AUX_A}'][...]
    auxdata2 = iter_group['auxdata']['${AUX_B}'][...]
    data = numpy.dstack((auxdata1, auxdata2))
    return data
EOF

# run w_assign to assign macrostates based off of defined BINS and STATES
w_assign -W $WDIR/$WEST --bins-from-file BINS --states-from-file STATES -o assign.h5 --construct-dataset module.pull_data_1d --serial &&
echo "    w_assign finished: assign.h5 file created in $SCHEME directory"

# run w_direct to then calculate multiple flux values and rate constants
# 'all' is equivalent to w_kinetics, w_kinavg, and w_stateprobs
    # or w_direct init/kinetics/probs
# step size of 1 is used, and can later be window averaged based off of ACF decay lag time, or can use default
# w_direct by default will set 'correl' = True, which calculates the ACF lag time and window averages
#w_direct init -a assign.h5 -W $WDIR/$WEST -o direct.h5 --step-iter 1 &&
#w_direct kinetics -a assign.h5 -W $WDIR/$WEST -o direct.h5 --evolution-mode cumulative --step-iter 1 &&

w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 --evolution-mode cumulative --step-iter 1 &&
#w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 --step-iter 1 &&
#w_direct all -a assign.h5 -W $WDIR/$WEST -o direct.h5 &&
echo "    w_direct finished: direct.h5 file created in $SCHEME directory"

