/**
 * Ephemeris-quality integrations
 *
 * This example uses the IAS15 integrator to integrate
 * the orbits of test particles in the field of the sun, moon, planets
 * and massive asteroids.  The positions and velocities of the
 * massive bodies are taken from JPL ephemeris files.  Solar GR
 * is included.
 *
 * To do:
 * 
 * 0. Write a wrapper function that takes an initial time, the initial positions and
 *    velocities of a set of one or more test particles, an integration time span or 
 *    final time, and time step.  The function should call the appropriate rebound
 *    routines and load the results into arrays that are accessible upon return (through
 *    pointers.  The return value of the function should represent success, failure, etc.
 * 
 * 1. Modify the code so that the initial conditions of the particle, the time
 *    span of the integration, and the time step come from a file.  We probably want to 
 *    allow the user to specific barycentric or geocentric.
 * 
 * 2. Rearrange the ephem() function so that it returns all the positions in one shot.
 * 
 * 3. Check position of the moon.  DONE.
 * 
 * 4. Separate ephem() and ast_ephem() from the rest of ephemeris_forces code.  DONE.
 * 
 * 5. Streamline ephem() function.  DONE.
 * 
 * 6. Put in earth J2 and J4.  This is a low priority.
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

double tmax;
void heartbeat(struct reb_simulation* r);

void ephem(const double G, const int i, const double t, double* const m,
	   double* const x, double* const y, double* const z,
	   double* const vx, double* const vy, double* const vz,
	   double* const ax, double* const ay, double* const az);


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->G = 0.295912208285591100E-03; // Gravitational constant (AU, solar masses, days)
    r->dt = 0.1;                    // time step in days
    r->integrator = REB_INTEGRATOR_IAS15;
    r->heartbeat = heartbeat;
    r->display_data = NULL;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity = REB_GRAVITY_NONE;
    r->usleep = 1000.;

    r->t = 2458849.5; // set simulation internal time to the time of test particle initial conditions.

    struct reb_particle tp = {0};

    // Initial conditions for Ceres from JPL Horizons
    // Barycentric equatorial frame
    
    /* Ceres
    tp.x = 1.003810191255037E+00;
    tp.y =-2.383247476572548E+00;
    tp.z =-1.329143955066118E+00;
    tp.vx = 9.193372298255753E-03;
    tp.vy = 3.368462462504294E-03;
    tp.vz = -2.856144947515055E-04;
    */

    // Holman
    /*
    tp.x =  3.338876057509365E+00;
    tp.y =  -9.176517956664152E-01;
    tp.z = -5.038590450387491E-01;
    tp.vx = 2.805663678557796E-03;
    tp.vy = 7.550408259144305E-03;
    tp.vz = 2.980028369986096E-03;
    */

    // 2020 CD3
    /* heliocentric
    tp.x  = -1.743654025843151E-01;
    tp.y  = 8.883058030225525E-01;
    tp.z  = 3.953036087753029E-01;
    tp.vx = -1.722157767734381E-02;
    tp.vy = -2.728117293034965E-03;
    tp.vz = -1.109706882773078E-03;

    // Shift to geocenter
    double xe, ye, ze, vxe, vye, vze, axe, aye, aze, m;
    ephem(r->G, 3, r->t, &m, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze); // Get position and mass of the earth wrt barycenter
    printf("%le %le %le %le %le %le %le %le %le\n", xe, ye, ze, vxe, vye, vze, axe, aye, aze);
    tp.x -= xe;
    tp.y -= ye;
    tp.z -= ze;
    tp.vx -= vxe;
    tp.vy -= vye;
    tp.vz -= vze;

    */

    // 2020 CD3
    /* geocentric */
    tp.x  =  -4.222728794111765E-03;
    tp.y  =  -7.672121501277942E-03;
    tp.z  =   6.861481718947984E-03;
    tp.vx =   2.589170139913478E-05;
    tp.vy =   9.226590401423815E-06;
    tp.vz =   7.662476921885610E-05;
    
    reb_add(r, tp);

    struct rebx_extras* rebx = rebx_attach(r);

    // Also add "ephemeris_forces" 
    struct rebx_force* ephem_forces = rebx_load_force(rebx, "ephemeris_forces");
    rebx_add_force(rebx, ephem_forces);

    // Have to set speed of light in right units (set by G & initial conditions).  Here we use default units of AU/(yr/2pi)
    rebx_set_param_int(rebx, &ephem_forces->ap, "N_ephem", 11);
    rebx_set_param_int(rebx, &ephem_forces->ap, "N_ast", 16);
    rebx_set_param_double(rebx, &ephem_forces->ap, "c", 173.144632674);

    tmax            = r->t + 200.;

    // clearing out the file
    FILE* g = fopen("states.txt","w");
    fclose(g);
    
    reb_integrate(r, tmax);

}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 0.5)){
        reb_output_timing(r, tmax);
        reb_integrator_synchronize(r);

	double xe, ye, ze, vxe, vye, vze, axe, aye, aze, m;	
	ephem(r->G, 3, r->t, &m, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze); // Get position and mass of the earth wrt barycenter	
        FILE* g = fopen("states.txt","a");
        fprintf(g,"%lf\t",r->t);
	const int N = r->N;	
	for (int i=0;i<N;i++){
	  struct reb_particle p = r->particles[i];
	  fprintf(g,"%e\t%e\t%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
	  //fprintf(g,"%e\t%e\t%e\t%e\t%e\t%e\n",
	  //p.x-xe, p.y-ye, p.z-ze, p.vx-vxe, p.vy-vye, p.vz-vze);	  
	}
        fclose(g);
	//reb_output_ascii(r, "states.txt");
    }
}


