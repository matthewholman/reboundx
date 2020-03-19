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

rebx_lib = CDLL('/Users/mholman/reboundx/examples/ephem_forces/libreboundx.so')

def integration_function(tstart, tstep, trange, geocentric,
                         xi, yi, zi,
                         vxi, vyi, vzi):

    n_out = int(8*np.abs(trange/tstep));
    nout = c_int(n_out)

    OutArray = Tstate * n_out
    outstate_arr = OutArray()


    _integration_function = rebx_lib.integration_function
    _integration_function.argtypes = (c_double, c_double, c_double, c_int,
                                      c_double, c_double, c_double,
                                      c_double, c_double, c_double,                                      
                                      POINTER(OutArray), POINTER(c_int))
    _integration_function.restype = None

    return_value = _integration_function(tstart, tstep, trange, geocentric, 
                                         xi, yi, zi, vxi, vyi, vzi,
                                         byref(outstate_arr),
                                         byref(nout))

    return outstate_arr, nout.value

