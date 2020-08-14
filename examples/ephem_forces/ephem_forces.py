"""This is a python wrapper to an ephemeris-quality integrator function.
This wrapper uses ctypes to access a reboundx c library that contains an
extension for carrying out highly accurate integrations of test
particles moving in the field of the sun, planets, moon, and massive
asteroids, with positions and velocities supplied by JPL ephemeris 
files.  
"""

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
        ('n_out', c_int),        
        ('n_particles', c_int)
    ]

def integration_function(tstart, tstep, trange,
                         geocentric,
                         n_particles,
                         instate_arr):
    """Wrapper for the c integration_function.
    Parameters
    ----------
    tstart : float
         The starting time of the integration in JD (TDB).
    tstep :  float
         Suggested time step for the integration, also in JD (TDB).
    trange : float
         The total amount of time to be integrated.
    geocentric : int
         Indicates if the initial conditions are geocentric.
    n_particles : int
         Number of test particles to be integrated.
    instate_arr : numpy array of float (compatible with c doubles).
         This a 1-d array of floats, giving initial (x, y, z, vx, vy, vz) 
         for each test particle.
    Returns
    -------
    np.array 
         an array of the times (JD TDB) in the output.
    np.array 
         an array of the dynamical states (positions and velocities) of
         the test particles at each of the output times.  The array has
         shape (number of times, n_particles, 6).
    int
         number of output times         
    int
         number of particles
    """

    # Instantiate a TimeState structure to hold the integration results.
    timestate = TimeState()

    # Set up call to integration_function
    _integration_function = rebx_lib.integration_function

    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_int,
                                      POINTER(c_double),
                                      POINTER(TimeState))

    #_integration_function.restype = None

    return_value = _integration_function(tstart, tstep, trange, geocentric,
                                         n_particles,
                                         instate_arr.ctypes.data_as(POINTER(c_double)),
                                         byref(timestate))

    # Parse and restructure the results
    n_out = timestate.n_out
    times  = np.ctypeslib.as_array(timestate.t, shape=(n_out,))
    states = np.ctypeslib.as_array(timestate.state, shape=(n_out, 7*n_particles, 6))
    n_particles = timestate.n_particles

    return times, states, n_out, n_particles
