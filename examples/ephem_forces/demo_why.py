import numpy as np
import ephem_forces

tstart, tstep, trange = 2458849.5, 20.0, 500
geocentric = 0
n_particles = 1
instates = np.array([3.338876057509365E+00, -9.176517956664152E-01, -5.038590450387491E-01,
                     2.805663678557796E-03, 7.550408259144305E-03, 2.980028369986096E-03])

times, states, n_out, n_particles = ephem_forces.integration_function(tstart, tstep, trange, geocentric, n_particles, instates)

outstate0 = states[0, 0, :]
print(instates - outstate0)



