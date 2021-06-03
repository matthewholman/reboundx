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

typedef struct {
    double* t;
    double* state;
    int n_out;
    int n_particles;
} timestate;

void read_inputs(char *filename, double* tepoch, double* tstart, double* tend, double* tstep, 
		 int *geocentric,
		 double *epsilon,
		 int *n_particles,		 
		 double **instate,
		 int *n_var,
		 int **invar_part,		 
		 double **invar,
		 double **cov_mat);

int main(int argc, char* argv[]){

    //timestate ts;
    //timestate* tsp = (timestate*)malloc(sizeof(timestate));

    double *instate;
    int *invar_part;            
    double *invar;        
    double* outstate = NULL;
    double* outtime = NULL;
    double* cov_mat;

    int n_alloc = 10000;
    int N = 10;
    outstate = (double *) malloc((8*n_alloc+1)*N*6*sizeof(double));
    outtime  = (double *) malloc((8*n_alloc+1)*sizeof(double));
    
    //tsp->t = outtime;
    //tsp->state = outstate;
    //tsp->n_out = n_alloc;    

    int n_out;

    // Read ICs & integration params from file
    double tepoch, tstart, tend, tstep;
    int geocentric;
    int n_particles;
    int n_var;
    double epsilon;

    if(argc >=2){
	read_inputs(argv[1], &tepoch, &tstart, &tend, &tstep,
		    &geocentric, &epsilon,
		    &n_particles,		    
		    &instate,
		    &n_var,
		    &invar_part,
		    &invar,
		    &cov_mat);
    }else{
//	read_inputs("initial_conditions.txt", &tepoch, &tstart, &tstep, &trange, &geocentric, &instate, &n_particles);
        printf("No Input File\n");
        exit(EXIT_FAILURE);
    }

    double scale;

    sscanf(argv[2], "%lf", &scale);


    int integration_function(double tstart, double tend, double tstep, 
			     int geocentric,
			     double epsilon,
			     int n_particles,
			     double* instate,
			     int n_var,
			     int* invar_part,			 
			     double* invar,
			     int n_alloc,			 
			     int *n_out,
			     double* outtime,
			     double* outstate);			 

    // clearing out the file
    FILE* g = fopen("out_states.txt","w");


    if(tstart >= tepoch){
	int status;

	status = integration_function(tepoch, tend, tstep, 
				      geocentric,
				      epsilon,
				      n_particles,
				      instate,
				      n_var,
				      invar_part,				      
				      invar,
				      n_alloc,
				      &n_out,
				      outtime,
				      outstate);

	for(int i=0; i<(8*n_out+1); i++){
	    int offset = i*(n_particles+n_var)*6; //XYZ
	    for(int j=0; j<(n_particles+n_var); j++){
		fprintf(g,"%lf ", outtime[i]);
		fprintf(g,"%3d ", j);
		for(int k=0; k<6; k++){
		    fprintf(g,"%23.16e ", outstate[offset+6*j+k]);
		}
		fprintf(g,"\n");
		/*
		fprintf(g,"%lf ", outtime[i]);
		fprintf(g,"%3d ", j);
		for(int k=3; k<6; k++){
		    fprintf(g,"%23.16e ", outstate[offset+6*j+k]);
		}
		fprintf(g,"\n");
		*/
	    }
	}
     
    }else{

	int status = integration_function(tepoch, -tstep, tstart,
					  geocentric,
					  epsilon,
					  n_particles,
					  instate,
					  n_var,
					  invar_part,
					  invar,
					  n_alloc,					  
					  &n_out,
					  outtime,
					  outstate);

	for(int i=(8*n_out); i>0; i--){
	    int offset = i*(n_particles+n_var)*6; //XYZ
	    for(int j=0; j<(n_particles+n_var); j++){
		fprintf(g,"%lf ", outtime[i]);
		fprintf(g,"%3d ", j);
		for(int k=0; k<6; k++){
		    fprintf(g,"%23.16e ", outstate[offset+6*j+k]);
		}
		fprintf(g,"\n");
		/*
		fprintf(g,"%lf ", outtime[i]);
		fprintf(g,"%3d ", j);
		for(int k=3; k<6; k++){
		    fprintf(g,"%23.16e ", outstate[offset+6*j+k]);
		}
		fprintf(g,"\n");
		*/
	    }
	}

	status = integration_function(tepoch, tstep, tend,
				      geocentric,
				      epsilon,
				      n_particles,
				      instate,
				      n_var,
				      invar_part,
				      invar,
				      n_alloc,					  				      
				      &n_out,
				      outtime,
				      outstate);

	for(int i=0; i<(8*n_out+1); i++){
	    int offset = i*(n_particles+n_var)*6; //XYZ
	    for(int j=0; j<(n_particles+n_var); j++){
		fprintf(g,"%lf ", outtime[i]);
		fprintf(g,"%3d ", j);
		for(int k=0; k<6; k++){
		    fprintf(g,"%23.16e ", outstate[offset+6*j+k]);
		}
		fprintf(g,"\n");
		/*
		fprintf(g,"%lf ", outtime[i]);
		fprintf(g,"%3d ", j);
		for(int k=3; k<6; k++){
		    fprintf(g,"%23.16e ", outstate[offset+6*j+k]);
		}
		fprintf(g,"\n");
		*/
		
	    }
	
	}
    }

    fclose(g);    

}

void read_inputs(char *filename, double* tepoch, double* tstart, double* tend, double* tstep, 
		 int *geocentric,
		 double *epsilon,
		 int *n_particles,		 
		 double **instate,
		 int *n_var,
		 int **invar_part,		 
		 double **invar,
		 double **cov_mat){

     char label[100]; /* hardwired for length */
     FILE* fp;

     int np = 0;
     int nvar = 0;

     int n_alloc = 1;
     double* state = malloc(n_alloc*6*sizeof(double));

     int n_var_alloc = 1;
     double* var = malloc(n_var_alloc*6*sizeof(double));
     int* var_part = malloc(n_var_alloc*sizeof(int));

     double* cov = malloc(36*sizeof(double));
     
     if((fp = fopen(filename, "r")) != NULL){

      while(fscanf(fp, "%s", label) != EOF){
        if(!strcmp(label, "tepoch")){
	  fscanf(fp, "%lf", tepoch);     
        } else if(!strcmp(label, "tstart")){
          fscanf(fp, "%lf", tstart);
        } else if(!strcmp(label, "tstep")){
	  fscanf(fp, "%lf", tstep);
        } else if(!strcmp(label, "tend")){
	  fscanf(fp, "%lf", tend);
        } else if(!strcmp(label, "epsilon")){
         fscanf(fp, "%lf", epsilon);
        } else if(!strcmp(label, "geocentric")){
         fscanf(fp, "%d", geocentric);
        } else if(!strcmp(label, "state")){
         fscanf(fp, "%lf%lf%lf", &state[6*np+0], &state[6*np+1], &state[6*np+2]);
         fscanf(fp, "%lf%lf%lf", &state[6*np+3], &state[6*np+4], &state[6*np+5]);	 
         np++;

	 // Resize the array, if needed.
	 if(np==n_alloc){
	     n_alloc *= 2;
	     state = realloc(state, n_alloc*6*sizeof(double));
	 }
	 
        } else if(!strcmp(label, "var")){
	    fscanf(fp, "%d", &var_part[nvar]);
	    fscanf(fp, "%lf%lf%lf", &var[6*nvar+0], &var[6*nvar+1], &var[6*nvar+2]);
	    fscanf(fp, "%lf%lf%lf", &var[6*nvar+3], &var[6*nvar+4], &var[6*nvar+5]);	 
	    nvar++;

	 // Resize the array, if needed.
	 if(nvar==n_var_alloc){
	     n_var_alloc *= 2;
	     var = realloc(var, n_var_alloc*6*sizeof(double));
	     var_part = realloc(var_part, n_var_alloc*sizeof(int));	     
	 }
	 
        } else if(!strcmp(label, "covariance")){
         fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[0], &cov[1], &cov[2],&cov[3], &cov[4], &cov[5]);
         fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[6], &cov[7], &cov[8],&cov[9], &cov[10], &cov[11]);
         fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[12], &cov[13], &cov[14],&cov[15], &cov[16], &cov[17]);
         fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[18], &cov[19], &cov[20],&cov[21], &cov[22], &cov[23]);
         fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[24], &cov[25], &cov[26],&cov[27], &cov[28], &cov[29]);
         fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[30], &cov[31], &cov[32],&cov[33], &cov[34], &cov[35]);
 
        } else {
         printf("No label: %s\n", label);
         exit(EXIT_FAILURE);
        }
      }

      //deallocate unused space
      state = realloc(state, np*6*sizeof(double));
      var   = realloc(var, nvar*6*sizeof(double));
      var_part   = realloc(var_part, nvar*6*sizeof(int));            

      *n_particles = np;
      *n_var = nvar;      
      *instate = state;      
      *invar = var;
      *invar_part = var_part;            
      *cov_mat = cov;

      fclose(fp);

     }else{
       exit(EXIT_FAILURE);
     }

     return;
}

