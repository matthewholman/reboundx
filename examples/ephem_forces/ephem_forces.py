"""This is a python wrapper to an ephemeris-quality integrator function.
This wrapper uses ctypes to access a reboundx c library that contains an
extension for carrying out highly accurate integrations of test
particles moving in the field of the sun, planets, moon, and massive
asteroids, with positions and velocities supplied by JPL ephemeris 
files.  
"""

import ctypes
from ctypes import *
import numpy as np
from os import path as osp

rebx_location = osp.join(osp.dirname(osp.realpath(__file__)), 'libreboundx.so')
rebx_lib = CDLL(rebx_location)

class TimeState(Structure):
    """
    A ctypes mapping to the structure populated by integration_function.
    """
    _fields_ = [
        ('t', POINTER(c_double)),
        ('state', POINTER(c_double)),
        ('n_alloc', c_int),        
        ('n_particles', c_int)
    ]

def integration_function(tstart, tstep, trange,
                         geocentric,
                         n_particles,
                         instate_arr):

    n_alloc = 1000

    tsize = (n_alloc*8+1)    
    array_of_tsize_doubles = c_double*tsize

    ssize = (n_alloc*8+1)*6*n_particles*7
    array_of_ssize_doubles = c_double*ssize    

    outtime  = array_of_tsize_doubles()
    outstate = array_of_ssize_doubles()

    n_out = c_int()    

    # Set up call to integration_function
    _integration_function = rebx_lib.integration_function

    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_int,
                                      POINTER(c_double),
                                      c_int,
                                      POINTER(c_int),
                                      POINTER(array_of_tsize_doubles),
                                      POINTER(array_of_ssize_doubles))
    
    _integration_function.restype = c_int

    return_value = _integration_function(tstart, tstep, trange,
                                         geocentric,
                                         n_particles,
                                         instate_arr.ctypes.data_as(POINTER(c_double)),
                                         n_alloc,
                                         byref(n_out),
                                         byref(outtime),
                                         byref(outstate))


    # Parse and restructure the results

    #times  = np.ctypeslib.as_array(timestate.t, shape=(n_out*8+1,))
    #states = np.ctypeslib.as_array(timestate.state, shape=(n_out*8+1, 7*n_particles, 6))
    #n_particles = timestate.n_particles

    return outtime, outstate, n_out.value, n_particles, return_value

