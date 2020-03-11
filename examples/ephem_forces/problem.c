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
 *    allow the user to specific barycentric or geocentric. DONE.
 * 
 * 2. Rearrange the ephem() function so that it returns all the positions in one shot.
 * 
 * 3. Check position of the moon.  DONE.
 * 
 * 4. Separate ephem() and ast_ephem() from the rest of ephemeris_forces code.  DONE.
 * 
 * 5. Streamline ephem() function.  DONE.
 * 
 * 6. Put in earth J2 and J4.  DONE.  Could put in the orientation of the spin 
 *    axis.  Can't include both J2/J4 of earth and sun without this.
 * 
 * 7. Put in J2 of the sun.  This will require thinking about the orientation of the spin 
 *    axis.
 *
 *
 */
#include <stdio.h>
#include <string.h>
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

void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* xi, double* yi, double* zi,
		 double* vxi, double* vyi, double* vzi);


void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* xi, double* yi, double* zi,
		 double* vxi, double* vyi, double* vzi){


     char label[100]; /* hardwired for length */  
     FILE* fp;

     if((fp = fopen(filename, "r")) != NULL){

      while(fscanf(fp, "%s", label) != EOF){
        if(!strcmp(label, "tstart")){
  	 fscanf(fp, "%lf", tstart);     
        } else if(!strcmp(label, "tstep")){
 	 fscanf(fp, "%lf", tstep);
        } else if(!strcmp(label, "trange")){
 	 fscanf(fp, "%lf", trange);
        } else if(!strcmp(label, "geocentric")){
 	 fscanf(fp, "%d", geocentric);
        } else if(!strcmp(label, "state")){
 	 fscanf(fp, "%lf%lf%lf", xi, yi, zi);	 
 	 fscanf(fp, "%lf%lf%lf", vxi, vyi, vzi);
        } else {
 	 printf("No label: %s\n", label);
	 exit(EXIT_FAILURE);
        }
      }

      fclose(fp);

     }else{
       exit(EXIT_FAILURE);       
     }

}


int main(int argc, char* argv[]){

    struct reb_particle* outstate = malloc(10000*sizeof(struct reb_particle));
    double* outtime = malloc(10000*sizeof(double));

    outtime[100] = 1.0;
    outstate[100].x = 1.0;    

    int n_out;
  
    void integration_function(double tstart, double tstep, double trange, int geocentric,
			      double xi, double yi, double zi, double vxi, double vyi, double vzi,
			      struct reb_particle *outstate,
			      double* outtime,
			      int* n_out);			      
  
    // Read ICs & integration params from file
    double tstart, tstep, trange;
    double xi, yi, zi, vxi, vyi, vzi;
    int geocentric;

    if(argc >=2){
      read_inputs(argv[1], &tstart, &tstep, &trange, &geocentric, &xi, &yi, &zi, &vxi, &vyi, &vzi);
    }else{
      read_inputs("initial_conditions.txt", &tstart, &tstep, &trange, &geocentric, &xi, &yi, &zi, &vxi, &vyi, &vzi);
    }

    integration_function(tstart, tstep, trange, geocentric,
			 xi, yi, zi, vxi, vyi, vzi,
			 outstate, outtime, &n_out);

    // clearing out the file
    FILE* g = fopen("out_states.txt","w");

    for(int i=0; i<n_out; i++){
      fprintf(g,"%lf\t",outtime[i]);
      struct reb_particle p = outstate[i];
      fprintf(g,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz);
    }
    fclose(g);

}


void integration_function(double tstart, double tstep, double trange, int geocentric,
			   double xi, double yi, double zi, double vxi, double vyi, double vzi,
			  struct reb_particle *outstate,
			  double* outtime,
			  int* n_out){
			  
    struct reb_simulation* r = reb_create_simulation();

    // Set up simulation constants
    r->G = 0.295912208285591100E-03; // Gravitational constant (AU, solar masses, days)
    r->integrator = REB_INTEGRATOR_IAS15;
    r->heartbeat = heartbeat;
    r->display_data = NULL;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity = REB_GRAVITY_NONE;
    r->usleep = 20000.;

    struct rebx_extras* rebx = rebx_attach(r);

    // Also add "ephemeris_forces" 
    struct rebx_force* ephem_forces = rebx_load_force(rebx, "ephemeris_forces");
    rebx_add_force(rebx, ephem_forces);

    rebx_set_param_int(rebx, &ephem_forces->ap, "geocentric", geocentric);

    // Set number of ephemeris bodies
    rebx_set_param_int(rebx, &ephem_forces->ap, "N_ephem", 11);

    // Set number of massive asteroids
    rebx_set_param_int(rebx, &ephem_forces->ap, "N_ast", 16);

    // Set speed of light in right units (set by G & initial conditions).
    // Here we use default units of AU/(yr/2pi)
    rebx_set_param_double(rebx, &ephem_forces->ap, "c", 173.144632674);

    rebx_set_param_int(rebx, &ephem_forces->ap, "n_out", 0);
    rebx_set_param_pointer(rebx, &ephem_forces->ap, "outstate", outstate);
    rebx_set_param_pointer(rebx, &ephem_forces->ap, "outtime", outtime);    

    struct reb_particle tp = {0};

    tp.x  =  xi;
    tp.y  =  yi;
    tp.z  =  zi;
    tp.vx =  vxi;
    tp.vy =  vyi;
    tp.vz =  vzi;
    
    reb_add(r, tp);

    r->t = tstart;    // set simulation internal time to the time of test particle initial conditions.
    tmax  = r->t + trange;
    r->dt = tstep;                   // time step in days

    reb_integrate(r, tmax);

    int *no = rebx_get_param(rebx, ephem_forces->ap, "n_out");    
    *n_out = *no;

}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 0.5)){
        reb_output_timing(r, tmax);
        //reb_integrator_synchronize(r);

	struct rebx_extras* rebx = r->extras;
	struct rebx_force* ephem_forces = rebx_get_force(rebx, "ephemeris_forces");

	int *n_out = rebx_get_param(rebx, ephem_forces->ap, "n_out");
	struct reb_particle* outstate = rebx_get_param(rebx, ephem_forces->ap, "outstate");
	double* outtime  = rebx_get_param(rebx, ephem_forces->ap, "outtime");
	
	outstate[*n_out].x = r->particles[0].x;
	outstate[*n_out].y = r->particles[0].y;
	outstate[*n_out].z = r->particles[0].z;
	outstate[*n_out].vx = r->particles[0].vx;
	outstate[*n_out].vy = r->particles[0].vy;
	outstate[*n_out].vz = r->particles[0].vz;

	outtime[*n_out] = r->t;

	*n_out += 1;

	rebx_set_param_int(rebx, &ephem_forces->ap, "n_out", *n_out);


    }
}


