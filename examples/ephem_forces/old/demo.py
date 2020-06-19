import numpy as np
import pandas as pd
import ephem_forces

df = pd.read_csv('out_states.txt', sep='\s+', names=['jd','x','y','z','vx','vy','vz'])

tstart, tstep, trange = 2458849.5, 20.0, 500
geocentric = 0
xi, yi, zi = 3.338876057509365E+00, -9.176517956664152E-01, -5.038590450387491E-01
vxi, vyi, vzi = 2.805663678557796E-03, 7.550408259144305E-03, 2.980028369986096E-03
states, nout = ephem_forces.integration_function(tstart, tstep, trange, geocentric,
                                                 xi, yi, zi,
                                                 vxi, vyi, vzi)
for i, state in enumerate(states):
    print(i, state.t, state.x, state.y, state.z)




