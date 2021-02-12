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
                         n_var,
                         invar_part,                         
                         invar,
                         epsilon = 1e-8):

    # Set up call to integration_function
    _integration_function = rebx_lib.integration_function

    _integration_function.argtypes = (c_double, c_double, c_double,
                                      c_int,
                                      c_double,
                                      c_int,
                                      POINTER(c_double),
                                      c_int,
                                      POINTER(c_int),
                                      POINTER(c_double),
                                      c_int,
                                      POINTER(c_int),
                                      POINTER(c_double),
                                      POINTER(c_double))
    
    _integration_function.restype = c_int

    n_alloc = int(abs(trange/tstep)*1.2)
    return_value = 5
    fac = 1.0
    iters = 0
    
    while(return_value == 5 and iters<5):

        n_alloc = int(fac*n_alloc)
        tsize = (n_alloc*8+1)    
        ssize = (n_alloc*8+1)*6*(n_particles+n_var)

        outtime = np.zeros((tsize), dtype=np.double)
        outstate = np.zeros((ssize), dtype=np.double)

        n_out = c_int()    

        return_value = _integration_function(tstart, tstep, trange,
                                             geocentric,
                                             epsilon,
                                             n_particles,
                                             instate_arr.ctypes.data_as(POINTER(c_double)),
                                             n_var,
                                             invar_part.ctypes.data_as(POINTER(c_int)),                                             
                                             invar.ctypes.data_as(POINTER(c_double)),
                                             n_alloc,
                                             byref(n_out),
                                             outtime.ctypes.data_as(POINTER(c_double)),
                                             outstate.ctypes.data_as(POINTER(c_double)))


        fac = int(1.5*abs(trange/(outtime[-1]-outtime[0])))
        iters += 1

    outstate = np.reshape(outstate, (-1, (n_particles+n_var), 6))
    outstate = outstate[:8*n_out.value+1]
    outtime = outtime[:8*n_out.value+1]

    states = outstate[:,0:n_particles,:]
    var_state = outstate[:,n_particles:,:]
    var_ng = None

    return outtime, states, var_state, var_ng, return_value    
    

def production_integration_function_wrapper(
        tstart,
        tstep,
        trange,
        epoch,
        geocentric,
        n_particles,
        instate_arr,
        non_grav_dict_list = None,
        epsilon=1e-8):
    """
    Standardized wrapper for calling integration_function
    Intended for "production" usage by MPC
     - Sets a bunch of standard defaults (e.g. tangent-eqn evolution vector specificaton)
    
    In the future it will handle non-gravs, but for now it does *** NOT ***
    Non-grav dicts will look like ...
        "nongrav_data": {
            "non_gravs": false,
            "booleans": {
                "yarkovski": false,
                "srp": false,
                "marsden": false,
                "yc": false,
                "yabushita": false,
                "A1": false,
                "A2": false,
                "A3": false,
                "DT": false
            },
            "coefficients": {
                "yarkovski": null,
                "srp": null,
                "A1": null,
                "A2": null,
                "A3": null,
                "DT": null
            }
        }
    """
    
    # We may have non-gravs ...
    if isinstance(non_grav_dict_list, list):

        pass
    
        # *** FOR NOW LET'S NOT WORRY ABOUT IMPLEMENTING THIS !!! ***
        #assert len(non_grav_dict_list) == n_particles
        #assert np.all( [isinstance(_, dict) for _ in non_grav_dict_list ] )
        ###do something to interpret each dict
        ###do they all need to be the same model?
        
    # variational eqn stuff
    n_var       = 6*n_particles
    invar_part  = np.repeat(np.arange(n_particles),6)
    invar       = np.concatenate([np.identity(6) for i in np.arange(n_particles)])

    # call the integrator
    # For non-gravs:
    # pass in an identifier for the model
    # pass in the parameters for that model
    # pass in variational particle states for the
    # non-grav parameters, in addition to the usual
    # variational particles
    return integration_function(tstart,
                                tstep,
                                trange,
                                geocentric,
                                n_particles,
                                instate_arr,
                                n_var,
                                invar_part,
                                invar,
                                epsilon = epsilon)
