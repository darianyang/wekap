#!/usr/bin/env python

import numpy

def pull_data_1d(n_iter, iter_group):
    '''
    This function reshapes auxiliary data for each iteration and returns it.
    '''
    auxdata = iter_group['auxdata']['o_angle_m1'][...]
    data = auxdata[:,:,numpy.newaxis]
    return data

def pull_data_2d(n_iter, iter_group):
    '''
    This function reshapes 2 auxiliary datasets for each iteration and returns it.
    '''
    #auxdata1 = iter_group['auxdata']['o_angle_m1'][...]
    # special for pcoord
    auxdata1 = iter_group['o_angle_m1'][...]
    auxdata2 = iter_group['auxdata']['o_angle_m2'][...]
    data = numpy.dstack((auxdata1, auxdata2))
    return data

def pull_data_3d(n_iter, iter_group):
    '''
    This function reshapes 3 auxiliary datasets for each iteration and returns it.
    '''
    # special for pcoord
    #auxdata1 = iter_group['o_angle_m1'][...]

    auxdata1 = iter_group['auxdata']['o_angle_m1'][...]
    auxdata2 = iter_group['auxdata']['o_angle_m2'][...]
    auxdata3 = iter_group['auxdata']['c2_angle'][...]
    data = numpy.dstack((auxdata1, auxdata2, auxdata3))
    return data
