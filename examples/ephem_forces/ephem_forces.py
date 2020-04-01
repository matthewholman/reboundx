from ctypes import *
import numpy as np

class TimeState(Structure):
    _fields_ = [
        ('t', POINTER(c_double)),
        ('state', POINTER(c_double)),
        ('n_out', c_int),        
        ('n_particles', c_int)
    ]
    
rebx_lib = CDLL('/Users/mholman/reboundx/examples/ephem_forces/libreboundx.so')

def integration_function(tstart, tstep, trange,
                         geocentric,
                         n_particles,
                         instate_arr):

    InStateArray = 6*n_particles*c_double

    timestate = TimeState()

    _integration_function = rebx_lib.integration_function
    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_int,
                                      POINTER(c_double),
                                      POINTER(TimeState))

    _integration_function.restype = None


    return_value = _integration_function(tstart, tstep, trange, geocentric,
                                         n_particles,
                                         instate_arr.ctypes.data_as(POINTER(c_double)),
                                         byref(timestate))

    n_out = timestate.n_out
    times  = np.ctypeslib.as_array(timestate.t, shape=(n_out,))
    states = np.ctypeslib.as_array(timestate.state, shape=(n_particles, 6, n_out))
    n_particles = timestate.n_particles

    return times, states, n_out, n_particles

