#!/usr/bin/env python

import numpy

def pull_data_1d(n_iter, iter_group):
    '''
    This function reshapes auxiliary data for each iteration and returns it.
    '''
    #auxdata = iter_group['auxdata']['pcoord'][...]
    auxdata = iter_group['pcoord'][...]
    data = auxdata[:,:,numpy.newaxis]
    return data

def pull_data_2d(n_iter, iter_group):
    '''
    This function reshapes 2 auxiliary datasets for each iteration and returns it.
    '''
    #auxdata1 = iter_group['auxdata']['pcoord'][...]
    # special for pcoord
    auxdata1 = iter_group['pcoord'][...]
    auxdata2 = iter_group['auxdata'][''][...]
    data = numpy.dstack((auxdata1, auxdata2))
    return data
