#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

typedef struct {
  double t, x, y, z, vx, vy, vz, ax, ay, az;
} tstate;

void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* xi, double* yi, double* zi,
		 double* vxi, double* vyi, double* vzi,
		 int* n_particles);

int main(int argc, char* argv[]){

    // Need to clean this up to make the array size
    // dynamic.
    tstate* outstate = malloc(100000*sizeof(tstate));
    
    int n_out;

    int integration_function(double tstart, double tstep, double trange, int geocentric,
			     double* xi, double* yi, double* zi, double* vxi, double* vyi, double* vzi,
			     int n_particles,
			     tstate* outstate, int* n_out);			      

    // Read ICs & integration params from file
    double tstart, tstep, trange;
    double *xi, *yi, *zi;
    double *vxi, *vyi, *vzi;
    int geocentric;
    int n_particles;

    //harcoded allocation for max particles
    xi = (double*)calloc(10,sizeof(double));
    yi = (double*)calloc(10,sizeof(double));
    zi = (double*)calloc(10,sizeof(double));
    vxi = (double*)calloc(10,sizeof(double));
    vyi = (double*)calloc(10,sizeof(double));
    vzi = (double*)calloc(10,sizeof(double));

    if(argc >=2){
      read_inputs(argv[1], &tstart, &tstep, &trange, &geocentric, xi, yi, zi, vxi, vyi, vzi, &n_particles);
    }else{
      read_inputs("initial_conditions.txt", &tstart, &tstep, &trange, &geocentric, xi, yi, zi, vxi, vyi, vzi, &n_particles);
    }

    //deallocate unused space
    xi = realloc(xi, ((n_particles)*sizeof(double)));
    yi = realloc(yi, ((n_particles)*sizeof(double)));
    zi = realloc(zi, ((n_particles)*sizeof(double)));
    vxi = realloc(vxi, ((n_particles)*sizeof(double)));
    vyi = realloc(vyi, ((n_particles)*sizeof(double)));
    vzi = realloc(vzi, ((n_particles)*sizeof(double)));

    integration_function(tstart, tstep, trange, geocentric,
			 xi, yi, zi, vxi, vyi, vzi,
			 n_particles,
			 outstate, &n_out);

    // clearing out the file
    FILE* g = fopen("out_states.txt","w");

    printf("%d %d\n", n_particles, n_out);
    //for(int i=0; i<n_out; i++){
    for(int i=0; i<n_out*n_particles; i++){    
      tstate p = outstate[i];
      fprintf(g,"%lf %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",p.t, p.x,p.y,p.z,p.vx,p.vy,p.vz);
    }
    fclose(g);    

}

void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* xi, double* yi, double* zi,
		 double* vxi, double* vyi, double* vzi, int* n_particles){

     char label[100]; /* hardwired for length */
     FILE* fp;

     *n_particles = 0;
     
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
         fscanf(fp, "%lf%lf%lf", &xi[*n_particles], &yi[*n_particles], &zi[*n_particles]);
         fscanf(fp, "%lf%lf%lf", &vxi[*n_particles], &vyi[*n_particles], &vzi[*n_particles]);
         (*n_particles)++;
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

