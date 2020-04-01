from ctypes import *
import numpy as np

class Tstate(Structure):
    _fields_ = [
        ('t', c_double),        
        ('x', c_double),
        ('y', c_double),
        ('z', c_double),
        ('vx', c_double),
	('vy', c_double),
	('vz', c_double),
        ('ax', c_double),
	('ay', c_double),
	('az', c_double)
    ]

class TimeState(Structure):
    _fields_ = [
        ('t', c_double),
        ('n_particles', c_int),
        ('n_out', c_int),        
        ('state', c_double)
    ]
    

rebx_lib = CDLL('/Users/mholman/reboundx/examples/ephem_forces/libreboundx.so')

def integration_function(tstart, tstep, trange,
                         geocentric,
                         n_particles,
                         instate_arr):

    fac = 5
    n_out = int(fac*8*np.abs(trange/tstep))
    nout = c_int(n_out)

    InStateArray = 6*n_particles*c_double

    States = TimeState()
    instate = InStateArray()
    #outtime_arr = np.zeros((n_out), dtype=np.double)        
    #outstate_arr = np.zeros((6*n_particles*n_out), dtype=np.double)    

    _integration_function = rebx_lib.integration_function
    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_int,
                                      POINTER(c_double),
                                      byref(instate))

    _integration_function.restype = None

    return_value = _integration_function(tstart, tstep, trange, geocentric,
                                         n_particles,
                                         instate_arr.ctypes.data_as(POINTER(c_double)),
                                         byref(timestate))

    return timestate

