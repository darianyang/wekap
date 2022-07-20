#!/usr/bin/env python

import numpy

def pull_data_1d(n_iter, iter_group):
    '''
    This function reshapes the progress coordinate and 
    auxiliary data for each iteration and returns it to
    the tool.
    '''
    auxdata = iter_group['auxdata']['1_75_39_c2'][...]
    data = auxdata[:,:,numpy.newaxis]
    return data

def pull_data_2d(n_iter, iter_group):
    '''
    This function reshapes the progress coordinate and 
    auxiliary data for each iteration and returns it to
    the tool.
    '''
    auxdata1 = iter_group['auxdata']['1_75_39_c2'][...]
    auxdata2 = iter_group['auxdata']['fit_m1_rms_heavy_m2'][...]
    data = numpy.dstack((auxdata1, auxdata2))
    return data
