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
		 double* vxi, double* vyi, double* vzi);


int main(int argc, char* argv[]){

    // Need to clean this up to make the array size
    // dynamic.
    tstate* outstate = malloc(1000000*sizeof(tstate));
    
    int n_out;

    int integration_function(double tstart, double tstep, double trange, int geocentric,
			      double xi, double yi, double zi, double vxi, double vyi, double vzi,
			      tstate *outstate, int* n_out);			      
  
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
			 outstate, &n_out);

    // clearing out the file
    FILE* g = fopen("out_states.txt","w");

    for(int i=0; i<n_out; i++){
      tstate p = outstate[i];
      fprintf(g,"%lf %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",p.t, p.x,p.y,p.z,p.vx,p.vy,p.vz);
    }
    fclose(g);

}



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
