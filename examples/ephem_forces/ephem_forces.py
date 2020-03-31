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

def keplerians(num, GM, state_arr):
    """
    Computes arrays of Keplerian orbital elements a, e, incl, longnode,
    argperi, and meananom, given a GM constant and an array of input states.
    
    *Returns*
    numpy arrays of a, e, incl, longnode, argperi, meananom
    
    """

    StateArray = State * num

    a_arr = np.zeros((num), dtype=np.double)
    e_arr = np.zeros((num), dtype=np.double)
    incl_arr = np.zeros((num), dtype=np.double)
    longnode_arr = np.zeros((num), dtype=np.double)
    argperi_arr = np.zeros((num), dtype=np.double)
    meananom_arr =np.zeros((num), dtype=np.double)

    _keplerians = lib.keplerians
    _keplerians.argtypes = (c_int, c_double, POINTER(StateArray), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double))
    _keplerians.restype = None


    return_value = _keplerians(num, GM, byref(state_arr),
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)))

    return a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr


def integration_function(tstart, tstep, trange,
                         geocentric,
                         n_particles,
                         instate_arr):

    fac = 5
    n_out = int(fac*8*np.abs(trange/tstep));
    nout = c_int(n_out)

    InStateArray = 6*n_particles*c_double

    outtime_arr = np.zeros((n_out), dtype=np.double)        
    outstate_arr = np.zeros((6*n_particles*n_out), dtype=np.double)    


    _integration_function = rebx_lib.integration_function
    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_int,
                                      POINTER(c_double),
                                      POINTER(c_int),
                                      POINTER(c_double),
                                      POINTER(c_double))

    _integration_function.restype = None

    return_value = _integration_function(tstart, tstep, trange, geocentric,
                                         n_particles,
                                         instate_arr.ctypes.data_as(POINTER(c_double)),
                                         byref(nout),
                                         outtime_arr.ctypes.data_as(POINTER(c_double)),
                                         outstate_arr.ctypes.data_as(POINTER(c_double)))

    return outtime_arr, outstate_arr, nout.value

