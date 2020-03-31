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
		 double* instate,
		 int* n_particles);

int main(int argc, char* argv[]){

    // Need to clean this up to make the array size
    // dynamic.
    double* instate = malloc(100000*6*sizeof(double));
    double* outstate = malloc(100000*6*sizeof(double));
    double* outtime  = malloc(100000*1*sizeof(double));    
    
    int n_out;

    // Read ICs & integration params from file
    double tstart, tstep, trange;
    //double *xi, *yi, *zi;
    //double *vxi, *vyi, *vzi;
    int geocentric;
    int n_particles;

    //harcoded allocation for max particles
    /*
    xi = (double*)calloc(1000, sizeof(double));
    yi = (double*)calloc(1000, sizeof(double));
    zi = (double*)calloc(1000, sizeof(double));
    vxi = (double*)calloc(1000, sizeof(double));
    vyi = (double*)calloc(1000, sizeof(double));
    vzi = (double*)calloc(1000, sizeof(double));
    */

    if(argc >=2){
	read_inputs(argv[1], &tstart, &tstep, &trange, &geocentric, instate, &n_particles);
    }else{
	read_inputs("initial_conditions.txt", &tstart, &tstep, &trange, &geocentric, instate, &n_particles);
    }

    //deallocate unused space
    instate = realloc(instate, n_particles*6*sizeof(double));

    int integration_function(double tstart, double tstep, double trange,
			     int geocentric,
			     int n_particles,
			     double* instate,
			     int* n_out, double* outtime, double* outstate);
    
    integration_function(tstart, tstep, trange,
			 geocentric,
			 n_particles,
			 instate,
			 &n_out, outtime, outstate);

    // clearing out the file
    FILE* g = fopen("out_states.txt","w");

    for(int i=0; i<n_out; i++){
	for(int j=0; j<n_particles; j++){
	    fprintf(g,"%lf ", outtime[i]);
	    fprintf(g,"%d ", j);
	    int offset = (i*n_particles+j)*6;
	    for(int k=0; k<6; k++){
		fprintf(g,"%16.8e ", outstate[offset+k]);
	    }
	    fprintf(g,"\n");
	}
    }
    fclose(g);    

}

void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* instate,
		 int* n_particles){

     char label[100]; /* hardwired for length */
     FILE* fp;

     int np = 0;
     
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
         fscanf(fp, "%lf%lf%lf", &instate[6*np+0], &instate[6*np+1], &instate[6*np+2]);
         fscanf(fp, "%lf%lf%lf", &instate[6*np+3], &instate[6*np+4], &instate[6*np+5]);	 
         np++;
        } else {
         printf("No label: %s\n", label);
         exit(EXIT_FAILURE);
        }
      }

     *n_particles = np;
      
      fclose(fp);

     }else{
       exit(EXIT_FAILURE);
     }

}

