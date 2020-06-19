#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* xi, double* yi, double* zi,
		 double* vxi, double* vyi, double* vzi, int* n_testicles);


int main(int argc, char* argv[]){

    // Read ICs & integration params from file
    double tstart, tstep, trange;
    double *xi, *yi, *zi;
    double *vxi, *vyi, *vzi;
    int geocentric;
    int n_testicles;
    
    xi = (double*)calloc(1,sizeof(double));
    yi = (double*)calloc(1,sizeof(double));
    zi = (double*)calloc(1,sizeof(double));
    vxi = (double*)calloc(1,sizeof(double));
    vyi = (double*)calloc(1,sizeof(double));
    vzi = (double*)calloc(1,sizeof(double));

    printf("Pointer: %p\n", xi);    

    if(argc >=2){
      read_inputs(argv[1], &tstart, &tstep, &trange, &geocentric, xi, yi, zi, vxi, vyi, vzi, &n_testicles);
    }else{
      read_inputs("initial_conditions.txt", &tstart, &tstep, &trange, &geocentric, xi, yi, zi, vxi, vyi, vzi, &n_testicles);
    }

    printf("Pointer: %p\n", xi);        

      printf("in calling routine:\n");
      printf("n_testicles: %d\n", n_testicles);
      printf("\n");
      for(int i = 0; i<n_testicles; i++){
       printf("%d %16.8e %16.8e %16.8e\n", i, xi[i], yi[i], zi[i]);
       printf("%d %16.8e %16.8e %16.8e\n", i, vxi[i], vyi[i], vzi[i]);
       printf("\n");
      }
}



void read_inputs(char *filename, double* tstart, double* tstep, double* trange,
		 int* geocentric, 
		 double* xi, double* yi, double* zi,
		 double* vxi, double* vyi, double* vzi, int* n_testicles){


     char label[100]; /* hardwired for length */  
     FILE* fp;
     int array_size = 0;

     printf("in read_inputs routine:\n");
     *n_testicles = 0;
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
	    array_size++;	    
	    if(array_size > 1){
		printf("realloc\n");
		printf("inside: %p\n", xi);        		
		xi = (double *)realloc(xi, ((array_size)*sizeof(double)));
		yi = (double *)realloc(yi, ((array_size)*sizeof(double)));
		zi = (double *)realloc(zi, ((array_size)*sizeof(double)));
		vxi = (double *)realloc(vxi, ((array_size)*sizeof(double)));
		vyi = (double *)realloc(vyi, ((array_size)*sizeof(double)));
		vzi = (double *)realloc(vzi, ((array_size)*sizeof(double)));

		printf("inside: %p\n", xi);        				

	    }
	    
	    fscanf(fp, "%lf%lf%lf", &xi[(array_size)-1], &yi[(array_size)-1], &zi[(array_size)-1]);	 
	    fscanf(fp, "%lf%lf%lf", &vxi[(array_size)-1], &vyi[(array_size)-1], &vzi[(array_size)-1]);
	    printf("array size: %d\n", array_size);
	    for(int i=0; i<array_size; i++){
		printf("%le %le %le %le %le %le\n", xi[i], yi[i], zi[i], vxi[i], vyi[i], vzi[i]);
	    }
	    printf("\n");
	    
        } else {
	    printf("No label: %s\n", label);
	    exit(EXIT_FAILURE);
        }
      }

      fclose(fp);

     }else{
       exit(EXIT_FAILURE);       
     }
     
     for(int i = 0; i<array_size; i++){
	 printf("%d %16.8e %16.8e %16.8e\n", i, xi[i], yi[i], zi[i]);
	 printf("%d %16.8e %16.8e %16.8e\n", i, vxi[i], vyi[i], vzi[i]);
	 printf("\n");
     }

     *n_testicles = array_size;

}
