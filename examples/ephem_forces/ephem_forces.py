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
                         instate_arr,
                         epsilon = 1e-8):

    # Set up call to integration_function
    _integration_function = rebx_lib.integration_function

    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_int,
                                      c_double,
                                      POINTER(c_double),
                                      c_int,
                                      POINTER(c_int),
                                      POINTER(c_double),
                                      POINTER(c_double))
    
    _integration_function.restype = c_int

    n_alloc = int(abs(trange/tstep)*1.2)
    return_value = 5
    fac = 1.0

    while(return_value == 5):

        n_alloc = int(fac*n_alloc)
        tsize = (n_alloc*8+1)    
        ssize = (n_alloc*8+1)*6*n_particles*7

        outtime = np.zeros((tsize), dtype=np.double)
        outstate = np.zeros((ssize), dtype=np.double)

        n_out = c_int()    

        return_value = _integration_function(tstart, tstep, trange,
                                             geocentric,
                                             n_particles,
                                             epsilon,
                                             instate_arr.ctypes.data_as(POINTER(c_double)),
                                             n_alloc,
                                             byref(n_out),
                                             outtime.ctypes.data_as(POINTER(c_double)),
                                             outstate.ctypes.data_as(POINTER(c_double)))

        outstate = np.reshape(outstate, (-1, 7*n_particles, 6))
        outstate = outstate[:8*n_out.value+1]
        outtime = outtime[:8*n_out.value+1]

        fac = int(1.5*abs(trange/(outtime[-1]-outtime[0])))        

    return outtime, outstate, n_out.value, n_particles, n_alloc, return_value    
    

