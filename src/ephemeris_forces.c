/** * @file ephemeris_forces.c
 * @brief   A routine for ephemeris-quality force calculations.
 * @author  Matthew Holman <mholman@cfa.harvard.edu>
 * @author  Arya Akmal <akmala@gmail.com>
 *
 *
 * Ephemeris-quality integrations
 *
 * This example uses the IAS15 integrator to integrate
 * the orbits of test particles in the field of the sun, moon, planets
 * and massive asteroids.  The positions and velocities of the
 * massive bodies are taken from JPL ephemeris files.  Solar GR, solar J2,
 * and earth J2/J4 are included.
 *
 * This is being developed to be incorporated into the reboundx package.
 *
 * Contributors:
 *
 * Robert Weryk <weryk@hawaii.edu>
 * Daniel Tamayo <dtamayo@astro.princeton.edu>
 * Matthew Payne <mpayne@cfa.harvard.edu>
 * David Hernandez <dmhernandez@cfa.harvard.edu>
 * Hanno Rein <hanno.rein@utoronto.ca>
 * Davide Farnocchia <davide.farnocchia@jpl.nasa.gov>
 * Jon Giorgini <jon.giorgini@jpl.nasa.gov>
 * 
 * To do:
 * 
 * 0. Write a wrapper function that takes an initial time, the initial positions and
 *    velocities of a set of one or more test particles, an integration time span or 
 *    final time, and time step.  The function should call the appropriate rebound
 *    routines and load the results into arrays that are accessible upon return (through
 *    pointers.  The return value of the function should represent success, failure, etc.
 *    DONE--mostly.
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
 *    axis.  Can't include both J2/J4 of earth and sun without this.  DONE.
 * 
 * 7. Put in J2 of the sun.  This will require thinking about the orientation of the spin 
 *    axis.  DONE.
 *
 * 8. Fix loop over particles.
 * 
 * 9. Develop sensible code that transitions to and from geocentric system.
 *
 * 10. Document a bunch of tests of the variational equation sections of the code. 
 *     The straight nbody routines for the planets and massive asteroids are likely
 *     ok, because they were copied from other parts of the rebound code.  The Earth J2/J4,
 *     solar J2, and solar GR sections need to be checked carefully.  I am not sure
 *     what kinds of tests to conduct.
 *
 *     Here's an idea: isolate each accelaration calculation and its associated 
 *     variational equation section and test them separately.
 *     
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "rebound.h"
#include "reboundx.h"

#include "spk.h"
#include "planets.h"

int ebody[11] = {
        PLAN_SOL,                       // Sun (in barycentric)
        PLAN_MER,                       // Mercury center
        PLAN_VEN,                       // Venus center
        PLAN_EAR,                       // Earth center
        PLAN_LUN,                       // Moon center
        PLAN_MAR,                       // Mars center
        PLAN_JUP,                       // ...
        PLAN_SAT,
        PLAN_URA,
        PLAN_NEP,
        PLAN_PLU
};

// Added gravitational constant G (2020 Feb 26)
// Added vx, vy, vz for GR stuff (2020 Feb 27)
// Consolidated the routine, removing the if block.
//
void ephem(const int i, const double jde, double* const GM,
	   double* const x, double* const y, double* const z,
	   double* const vx, double* const vy, double* const vz,
	   double* const ax, double* const ay, double* const az){

    static int initialized = 0;

    static struct _jpl_s *pl;
    struct mpos_s now;

    // The values below are G*mass.
    // Units are solar masses, au, days.
    const static double JPL_GM[11] =
	{
	    0.295912208285591100E-03, // 0  sun  
	    0.491248045036476000E-10, // 1  mercury
	    0.724345233264412000E-09, // 2  venus
	    0.888769244512563400E-09, // 3  earth
	    0.109318945074237400E-10, // 4  moon
	    0.954954869555077000E-10, // 5  mars
	    0.282534584083387000E-06, // 6  jupiter
	    0.845970607324503000E-07, // 7  saturn
	    0.129202482578296000E-07, // 8  uranus
	    0.152435734788511000E-07, // 9  neptune
	    0.217844105197418000E-11, // 10 pluto
	};

    if(i<0 || i>10){
      fprintf(stderr, "body out of range\n");
      exit(EXIT_FAILURE);
    }

    if (initialized == 0){
      
      if ((pl = jpl_init()) == NULL) {
	fprintf(stderr, "could not load DE430 file, fool!\n");
	exit(EXIT_FAILURE);
      }


      initialized = 1;

    }

    // Get position, velocity, and mass of body i in barycentric coords. 
    
    *GM = JPL_GM[i];

    jpl_calc(pl, &now, jde, ebody[i], PLAN_BAR); 

    // Convert to au/day and au/day^2
    vecpos_div(now.u, pl->cau);
    vecpos_div(now.v, pl->cau/86400.);
    vecpos_div(now.w, pl->cau/(86400.*86400.));

    *x = now.u[0];
    *y = now.u[1];
    *z = now.u[2];
    *vx = now.v[0];
    *vy = now.v[1];
    *vz = now.v[2];
    *ax = now.w[0];
    *ay = now.w[1];
    *az = now.w[2];
    
}

static void ast_ephem(const int i, const double jde, double* const GM, double* const x, double* const y, double* const z){

    static int initialized = 0;

    static struct spk_s *spl;
    struct mpos_s pos;

    // 1 Ceres, 4 Vesta, 2 Pallas, 10 Hygiea, 31 Euphrosyne, 704 Interamnia,
    // 511 Davida, 15 Eunomia, 3 Juno, 16 Psyche, 65 Cybele, 88 Thisbe, 
    // 48 Doris, 52 Europa, 451 Patientia, 87 Sylvia

    // The values below are G*mass.
    // Units are solar masses, au, days.      
    const static double JPL_GM[16] =
	{
	    1.400476556172344e-13, // ceres
	    3.854750187808810e-14, // vesta
	    3.104448198938713e-14, // pallas
	    1.235800787294125e-14, // hygiea
	    6.343280473648602e-15, // euphrosyne
	    5.256168678493662e-15, // interamnia
	    5.198126979457498e-15, // davida
	    4.678307418350905e-15, // eunomia
	    3.617538317147937e-15, // juno
	    3.411586826193812e-15, // psyche
	    3.180659282652541e-15, // cybele
	    2.577114127311047e-15, // thisbe
	    2.531091726015068e-15, // doris
	    2.476788101255867e-15, // europa
	    2.295559390637462e-15, // patientia
	    2.199295173574073e-15, // sylvia
	};
    

    if(i<0 || i>15){
      fprintf(stderr, "asteroid out of range\n");
      exit(EXIT_FAILURE);
    }

    if (initialized == 0){
      
      if ((spl = spk_init("sb431-n16s.bsp")) == NULL) {
	fprintf(stderr, "could not load sb431-n16 file, fool!\n");
	exit(EXIT_FAILURE);
      }

      
      initialized = 1;

    }

    *GM = JPL_GM[i];
    spk_calc(spl, i, jde, &pos);          
    *x = pos.u[0];
    *y = pos.u[1];
    *z = pos.u[2];
    
}

void rebx_ephemeris_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){

    const double G = sim->G;
    const double t = sim->t;
    const int N_var = sim->N_var;
    const int N_real= N;
//  const int N_var = sim->N - N;   //The N passed in is the number of real particles (sim->N - sim->N_var).

    const int* const N_ephem = rebx_get_param(sim->extras, force->ap, "N_ephem");
    if (N_ephem == NULL){
        fprintf(stderr, "REBOUNDx Error: Need to set N_ephem for ephemeris_forces\n");
        return;
    }
    
    const int* const N_ast = rebx_get_param(sim->extras, force->ap, "N_ast");
    if (N_ephem == NULL){
        fprintf(stderr, "REBOUNDx Error: Need to set N_ast for ephemeris_forces\n");
        return;
    }

    double* c = rebx_get_param(sim->extras, force->ap, "c");
    if (c == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n");
        return;
    }

    int* geo = rebx_get_param(sim->extras, force->ap, "geocentric"); // Make sure there is a default set.
    if (geo == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set geo flag.  See examples in documentation.\n");
        return;
    }

    const double C2 = (*c)*(*c);  // This could be stored as C2.
    
    double GM;
    double x, y, z, vx, vy, vz, ax, ay, az;
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
    double xe, ye, ze, vxe, vye, vze, axe, aye, aze;
    double xo, yo, zo, vxo, vyo, vzo;
    double xr, yr, zr, vxr, vyr, vzr;

    // Get mass, position, velocity, and acceleration of the Earth and Sun
    // for later use
    ephem(3, t, &GM, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze);
    ephem(0, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);     

    // The offset position is used to adjust the particle positions.
    if(*geo == 1){
      xo = xe;    yo = ye;    zo = ze;
      vxo = vxe;  vyo = vye;  vzo = vze;
    }else{
      xo = 0.0;  yo = 0.0;  zo = 0.0;
      vxo = 0.0; vyo = 0.0; vzo = 0.0;      
    }

    // Calculate acceleration due to sun and planets
    // Rearrange to get all the planet positions at once
    for (int i=0; i<*N_ephem; i++){

        // Get position and mass of massive body i.
        ephem(i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az); 

        for (int j=0; j<N_real; j++){
  	  // Compute position vector of test particle j relative to massive body i.
	  const double dx =  particles[j].x + (xo - x); 
	  const double dy =  particles[j].y + (yo - y);
	  const double dz =  particles[j].z + (zo - z);
	  const double _r = sqrt(dx*dx + dy*dy + dz*dz);
	  const double prefac = GM/(_r*_r*_r);

  	  particles[j].ax -= prefac*dx;
  	  particles[j].ay -= prefac*dy;
  	  particles[j].az -= prefac*dz;

        }
    }

    // Calculate acceleration of variational particles due to sun and planets
    for (int i=0; i<*N_ephem; i++){

        ephem(i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az); 

        for (int j=0; j<N_real; j++){ //loop over test particles

	    const double dx = particles[j].x + (xo - x);
	    const double dy = particles[j].y + (yo - y);
	    const double dz = particles[j].z + (zo - z);
	    const double r2 = dx*dx + dy*dy + dz*dz;
	    const double _r  = sqrt(r2);
	    const double r3inv = 1./(r2*_r);
	    const double r5inv = 3.*r3inv/r2;

	    // Coefficients for variational equations
	    const double dxdx = dx*dx*r5inv - r3inv;
	    const double dydy = dy*dy*r5inv - r3inv;
	    const double dzdz = dz*dz*r5inv - r3inv;
	    const double dxdy = dx*dy*r5inv;
	    const double dxdz = dx*dz*r5inv;
	    const double dydz = dy*dz*r5inv;

	    // Loop over variational particles
	    for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

		// Variational particle coords
		const double ddx = particles[v].x;
		const double ddy = particles[v].y;
		const double ddz = particles[v].z;
		const double Gmi = GM;

		// Matrix multiplication
		const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
		const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
		const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

		// No variational mass contributions for test particles!

		// Accumulate acceleration terms
		particles[v].ax += Gmi * dax; 
		particles[v].ay += Gmi * day; 
		particles[v].az += Gmi * daz; 

	    }
        }
    }


    // Calculate acceleration due to massive asteroids
    // Rearrange to get all the asteroid positions at once
    for (int i=0; i<*N_ast; i++){

        ast_ephem(i, t, &GM, &x, &y, &z); // Get position and mass of asteroid i.

	// Translate massive asteroids from heliocentric to barycentric.
	x += xs; y += ys; z += zs;
	
        for (int j=0; j<N_real; j++){
  	  // Compute position vector of test particle j relative to massive body i.
	    const double dx = particles[j].x + (xo - x);
	    const double dy = particles[j].y + (yo - y);
	    const double dz = particles[j].z + (zo - z); 	    	    
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
	    
            const double prefac = GM/(_r*_r*_r);
            particles[j].ax -= prefac*dx;
            particles[j].ay -= prefac*dy;
            particles[j].az -= prefac*dz;

        }
    }

    // Calculate acceleration of variational particles due to massive asteroids
    for (int i=0; i<*N_ast; i++){

        ast_ephem(i, t, &GM, &x, &y, &z); 

        for (int j=0; j<N_real; j++){ //loop over test particles

           const double dx = particles[j].x + (xo - x);
           const double dy = particles[j].y + (yo - y);
           const double dz = particles[j].z + (zo - z);
           const double r2 = dx*dx + dy*dy + dz*dz;
           const double _r  = sqrt(r2);
           const double r3inv = 1./(r2*_r);
           const double r5inv = 3.*r3inv/r2;

	   // Coefficients for variational equations
	   const double dxdx = dx*dx*r5inv - r3inv;
	   const double dydy = dy*dy*r5inv - r3inv;
	   const double dzdz = dz*dz*r5inv - r3inv;
	   const double dxdy = dx*dy*r5inv;
	   const double dxdz = dx*dz*r5inv;
	   const double dydz = dy*dz*r5inv;

	   // Loop over variational particles
           for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

	       // Variational particle coords
	       const double ddx = particles[v].x;
	       const double ddy = particles[v].y;
	       const double ddz = particles[v].z;
	       const double Gmi = GM;

	       // Matrix multiplication
	       const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
	       const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
	       const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

	       // No variational mass contributions for test particles!

	       particles[v].ax += Gmi * dax; 
	       particles[v].ay += Gmi * day; 
	       particles[v].az += Gmi * daz; 

           }
        }
    }

    // Here is the treatment of the Earth's J2 and J4.
    // Borrowed code from gravitational_harmonics example.
    // Assumes the coordinates are geocentric.
    // Also assuming that Earth's pole is along the z
    // axis.  This is only precisely true at the J2000
    // epoch.
    //

    // The geocenter is the reference for the J2/J4 calculations.
    xr = xe;  yr = ye;  zr = ze;

    // Hard-coded constants.  BEWARE!
    // Clean up on aisle 3!
    const double GMearth = 0.888769244512563400E-09;
    const double J2e = 0.00108262545;
    const double J4e = -0.000001616*0.0;
    const double au = 149597870.700;
    const double Re_eq = 6378.1263/au;
    // Unit vector to equatorial pole at the epoch
    // Clean this up!

    double RAs =  359.87123273*M_PI/180.;
    double Decs =  89.88809752*M_PI/180.;

    //double xp = cos(Decs)*cos(RAs);
    //double yp = cos(Decs)*sin(RAs);
    //double zp = sin(Decs);

    double xp =  0.0019111736356920146;
    double yp = -1.2513100974355823e-05;
    double zp =   0.9999981736277104;

    double incl = acos(zp);
    double longnode;
    if(xp != 0.0 || yp !=0.0) {    
      longnode = atan2(xp, -yp);
    } else {
      longnode = 0.0;
    }

    // Rearrange this loop for efficiency
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
	
	// Rotate to Earth equatorial frame
	// This could be a single rotation

	// Rotate around z by RA
	double cosr = cos(-longnode);
	double sinr = sin(-longnode);

	double dxp =  dx * cosr - dy * sinr;
	double dyp =  dx * sinr + dy * cosr;
	double dzp =  dz;

	// Rotate around x by Dec
	double cosd = cos(-incl);
	double sind = sin(-incl);
	
	dx =  dxp;
	dy =  dyp * cosd - dzp * sind;
	dz =  dyp * sind + dzp * cosd;

	// Calculate acceleration in
	// Earth equatorial frame	

	// J2 terms
        const double costheta2 = dz*dz/r2;
        const double J2e_prefac = 3.*J2e*Re_eq*Re_eq/r2/r2/r/2.;
        const double J2e_fac = 5.*costheta2-1.;
        const double J2e_fac2 = 7.*costheta2-1.;
        const double J2e_fac3 = 35.*costheta2*costheta2-30.*costheta2+3.;

	double resx = GMearth*J2e_prefac*J2e_fac*dx;
	double resy = GMearth*J2e_prefac*J2e_fac*dy;
	double resz = GMearth*J2e_prefac*(J2e_fac-2.)*dz;	

	// J4 terms
        const double J4e_prefac = 5.*J4e*Re_eq*Re_eq*Re_eq*Re_eq/r2/r2/r2/r/8.;
        const double J4e_fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;
        const double J4e_fac2= 33.*costheta2*costheta2-18.*costheta2 + 1.;
        const double J4e_fac3= 33.*costheta2*costheta2-30.*costheta2 + 5.;
        const double J4e_fac4= 231.*costheta2*costheta2*costheta2-315.*costheta2*costheta2+105.*costheta2 - 5.;

        resx += GMearth*J4e_prefac*J4e_fac*dx;
        resy += GMearth*J4e_prefac*J4e_fac*dy;
        resz += GMearth*J4e_prefac*(J4e_fac+12.-28.*costheta2)*dz;

	// Rotate back to original frame
	// Rotate around x by -Dec
	double resxp =  resx;
	double resyp =  resy * cosd + resz * sind;
	double reszp = -resy * sind + resz * cosd;
	
	// Rotate around z by -RA
	resx =  resxp * cosr + resyp * sinr;
	resy = -resxp * sinr + resyp * cosr;
	resz =  reszp;

	// Accumulate final acceleration terms
  	particles[j].ax += resx;
        particles[j].ay += resy; 
        particles[j].az += resz;

	// Constants for variational equations
	// J2 terms
	const double dxdx = GMearth*J2e_prefac*(J2e_fac-5.*J2e_fac2*dx*dx/r2);
	const double dydy = GMearth*J2e_prefac*(J2e_fac-5.*J2e_fac2*dy*dy/r2);
	const double dzdz = GMearth*J2e_prefac*(-1.)*J2e_fac3;
	const double dxdy = GMearth*J2e_prefac*(-5.)*J2e_fac2*dx*dy/r2;
	const double dydz = GMearth*J2e_prefac*(-5.)*(J2e_fac2-2.)*dy*dz/r2;
	const double dxdz = GMearth*J2e_prefac*(-5.)*(J2e_fac2-2.)*dx*dz/r2;
	// J4 terms
	const double dxdxJ4 = GMearth*J4e_prefac*(J4e_fac-21.*J4e_fac2*dx*dx/r2);
	const double dydyJ4 = GMearth*J4e_prefac*(J4e_fac-21.*J4e_fac2*dy*dy/r2);
	const double dzdzJ4 = GMearth*J4e_prefac*(-3.)*J4e_fac4;
	const double dxdyJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac2*dx*dy/r2;
	const double dydzJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac3*dy*dz/r2;
	const double dxdzJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac3*dx*dz/r2;

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){
          double ddx = particles[v].x;
          double ddy = particles[v].y;
          double ddz = particles[v].z;

	  // Rotate to Earth equatorial frame
  	  double ddxp =  ddx * cosr - ddy * sinr;
  	  double ddyp =  ddx * sinr + ddy * cosr;
  	  double ddzp =  ddz;
  	  ddx =  ddxp;
  	  ddy =  ddyp * cosd - ddzp * sind;
  	  ddz =  ddyp * sind + ddzp * cosd;

          // Matrix multiplication
          double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
          double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
          double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

          dax +=   ddx * dxdxJ4 + ddy * dxdyJ4 + ddz * dxdzJ4;
          day +=   ddx * dxdyJ4 + ddy * dydyJ4 + ddz * dydzJ4;
          daz +=   ddx * dxdzJ4 + ddy * dydzJ4 + ddz * dzdzJ4;

	  // Rotate back
  	  double daxp =  dax;
  	  double dayp =  day * cosd + daz * sind;
  	  double dazp = -day * sind + daz * cosd;
  	  dax =  daxp * cosr + dayp * sinr;
  	  day = -daxp * sinr + dayp * cosr;
  	  daz =  dazp;

	  // Accumulate acceleration terms
          particles[v].ax += dax;
          particles[v].ay += day;
          particles[v].az += daz;

        }
    }

    // Here is the treatment of the Sun's J2.
    // Borrowed code from gravitational_harmonics.
    // Also assuming that Sun's pole is along the z
    // axis.  This is not true.  Will correct soon.

    // The Sun center is reference for these calculations.
    xr = xs;  yr = ys;  zr = zs;

    // Hard-coded constants.  BEWARE!
    // Clean up on aisle 3!
    // Mass of sun in solar masses.    
    const double Msun = 1.0;  // hard-code parameter.
    const double Rs_eq = 696000.0/au;
    const double J2s = 2.1106088532726840e-07;

    RAs = 268.13*M_PI/180.;
    Decs = 63.87*M_PI/180.;

    xp = cos(Decs)*cos(RAs);
    yp = cos(Decs)*sin(RAs);
    zp = sin(Decs);

    incl = acos(zp);
    if(xp != 0.0 || yp !=0.0) {    
      longnode = atan2(xp, -yp);
    } else {
      longnode = 0.0;
    }
    
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);

	// Rotate to solar equatorial frame

	// Rotate around z by RA
	double cosr = cos(-longnode);
	double sinr = sin(-longnode);
	
	// Rotate around z by RA
	double dxp =  dx * cosr - dy * sinr;
	double dyp =  dx * sinr + dy * cosr;
	double dzp =  dz;

	// Rotate around x by Dec
	double cosd = cos(-incl);
	double sind = sin(-incl);

	dx =  dxp;
	dy =  dyp * cosd - dzp * sind;
	dz =  dyp * sind + dzp * cosd;

	const double costheta2 = dz*dz/r2;
        const double J2s_prefac = 3.*J2s*Rs_eq*Rs_eq/r2/r2/r/2.;
        const double J2s_fac = 5.*costheta2-1.;
        const double J2s_fac2 = 7.*costheta2-1.;
        const double J2s_fac3 = 35.*costheta2*costheta2-30.*costheta2+3.;

	// Calculate acceleration
	double resx = G*Msun*J2s_prefac*J2s_fac*dx;
	double resy = G*Msun*J2s_prefac*J2s_fac*dy;
	double resz = G*Msun*J2s_prefac*(J2s_fac-2.)*dz;

        // Variational equations

	// Rotate back to original frame
	// Rotate around x by -Dec
	double resxp =  resx;
	double resyp =  resy * cosd + resz * sind;
	double reszp = -resy * sind + resz * cosd;
	
	// Rotate around z by -RA
	resx =  resxp * cosr + resyp * sinr;
	resy = -resxp * sinr + resyp * cosr;
	resz =  reszp;

        particles[j].ax += resx;
        particles[j].ay += resy;
        particles[j].az += resz;

	// Constants for variational equations
	const double dxdx = G*Msun*J2s_prefac*(J2s_fac-5.*J2s_fac2*dx*dx/r2);
	const double dydy = G*Msun*J2s_prefac*(J2s_fac-5.*J2s_fac2*dy*dy/r2);
	const double dzdz = G*Msun*J2s_prefac*(-1.)*J2s_fac3;
	const double dxdy = G*Msun*J2s_prefac*(-5.)*J2s_fac2*dx*dy/r2;
	const double dydz = G*Msun*J2s_prefac*(-5.)*(J2s_fac2-2.)*dy*dz/r2;
	const double dxdz = G*Msun*J2s_prefac*(-5.)*(J2s_fac2-2.)*dx*dz/r2;

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

          double ddx = particles[v].x;
          double ddy = particles[v].y;
          double ddz = particles[v].z;

	  // Rotate to solar equatorial frame
	  double ddxp =  ddx * cosr - ddy * sinr;
	  double ddyp =  ddx * sinr + ddy * cosr;
	  double ddzp =  ddz;

          ddx =  ddxp;
	  ddy =  ddyp * cosd - ddzp * sind;
	  ddz =  ddyp * sind + ddzp * cosd;

          double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
          double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
          double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

          // Rotate back to original frame
	  double daxp =  dax;
	  double dayp =  day * cosd + daz * sind;
	  double dazp = -day * sind + daz * cosd;
	  dax =  daxp * cosr + dayp * sinr;
	  day = -daxp * sinr + dayp * cosr;
	  daz =  dazp;

	  // Accumulate acceleration terms
          particles[v].ax += dax;
          particles[v].ay += day;
          particles[v].az += daz;

        } 
       
    }

    // Here is the Solar GR treatment
    // The Sun is the reference for these calculations.    
    xr  = xs;  yr  = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;

    const double mu = G*Msun; 
    const int max_iterations = 10; // hard-coded parameter.
    for (int j=0; j<N_real; j++){

        struct reb_particle p = particles[j];
        struct reb_vec3d vi;

	p.x += (xo - xr);
	p.y += (yo - yr);
	p.z += (zo - zr);
	p.vx += (vxo - vxr);
	p.vy += (vyo - vyr);
	p.vz += (vzo - vzr);
	
        vi.x = p.vx;
        vi.y = p.vy;
        vi.z = p.vz;
        double vi2=vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	
        int q = 0;
        double A = (0.5*vi2 + 3.*mu/ri)/C2;
        struct reb_vec3d old_v;
        for(q=0; q<max_iterations; q++){
            old_v.x = vi.x;
            old_v.y = vi.y;
            old_v.z = vi.z;
            vi.x = p.vx/(1.-A);
            vi.y = p.vy/(1.-A);
            vi.z = p.vz/(1.-A);
            vi2 =vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
            A = (0.5*vi2 + 3.*mu/ri)/C2;
            const double dvx = vi.x - old_v.x;
            const double dvy = vi.y - old_v.y;
            const double dvz = vi.z - old_v.z;
            if ((dvx*dvx + dvy*dvy + dvz*dvz)/vi2 < DBL_EPSILON*DBL_EPSILON){
                break;
            }
        }
        const int default_max_iterations = 10;
        if(q==default_max_iterations){
            reb_warning(sim, "REBOUNDx Warning: 10 iterations in ephemeris forces failed to converge. This is typically because the perturbation is too strong for the current implementation.");
        }
  
        const double B = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/C2;
        const double rdotrdot = p.x*p.vx + p.y*p.vy + p.z*p.vz;
        
        struct reb_vec3d vidot;
        vidot.x = p.ax + B*p.x;
        vidot.y = p.ay + B*p.y;
        vidot.z = p.az + B*p.z;
        
        const double vdotvdot = vi.x*vidot.x + vi.y*vidot.y + vi.z*vidot.z;
        const double D = (vdotvdot - 3.*mu/(ri*ri*ri)*rdotrdot)/C2;

        particles[j].ax += B*(1.-A)*p.x - A*p.ax - D*vi.x;
        particles[j].ay += B*(1.-A)*p.y - A*p.ay - D*vi.y;
        particles[j].az += B*(1.-A)*p.z - A*p.az - D*vi.z;

	const double prefac = mu/(ri*ri*ri)/C2;
	const double rdotv = p.x*p.vx+p.y*p.vy+p.z*p.vz;
	const double fac1 = mu/ri-vi2;
	const double fac2 = 3.*vi2/ri/ri-4.*mu/ri/ri/ri;;
	const double fac3 = 12.*rdotv/ri/ri;

	const double dxdx = prefac*(fac1+fac2*p.x*p.x+4.*p.vx*p.vx-fac3*p.vx*p.x);
	const double dydy = prefac*(fac1+fac2*p.y*p.y+4.*p.vy*p.vy-fac3*p.vy*p.y);
	const double dzdz = prefac*(fac1+fac2*p.z*p.z+4.*p.vz*p.vz-fac3*p.vz*p.z);

	const double dxdy = prefac*(fac2*p.x*p.y+4.*p.vx*p.vy-fac3*p.vx*p.y);
	const double dydx = prefac*(fac2*p.y*p.x+4.*p.vy*p.vx-fac3*p.vy*p.x);
	const double dxdz = prefac*(fac2*p.x*p.z+4.*p.vx*p.vz-fac3*p.vx*p.z);

	const double dzdx = prefac*(fac2*p.z*p.x+4.*p.vz*p.vx-fac3*p.vz*p.x);
	const double dydz = prefac*(fac2*p.y*p.z+4.*p.vy*p.vz-fac3*p.vy*p.z);
	const double dzdy = prefac*(fac2*p.z*p.y+4.*p.vz*p.vy-fac3*p.vz*p.y);

	const double dxdvx = prefac*(4.*rdotv-2.*p.x*p.vx+4.*p.x*p.vx);
	const double dydvy = prefac*(4.*rdotv-2.*p.y*p.vy+4.*p.y*p.vy);
	const double dzdvz = prefac*(4.*rdotv-2.*p.z*p.vz+4.*p.z*p.vz);

	const double dxdvy = prefac*(-2.*p.x*p.vy+4.*p.y*p.vx);
	const double dydvx = prefac*(-2.*p.y*p.vx+4.*p.x*p.vy);
	const double dxdvz = prefac*(-2.*p.x*p.vz+4.*p.z*p.vx);

	const double dzdvx = prefac*(-2.*p.z*p.vx+4.*p.x*p.vz);
	const double dydvz = prefac*(-2.*p.y*p.vz+4.*p.z*p.vy);
	const double dzdvy = prefac*(-2.*p.z*p.vy+4.*p.y*p.vz);

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

	    // variational particle coords
	    double ddx = particles[v].x;
	    double ddy = particles[v].y;
	    double ddz = particles[v].z;
	    double ddvx = particles[v].vx;
	    double ddvy = particles[v].vy;
	    double ddvz = particles[v].vz;

	    // Matrix multiplication
	    const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		+   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
	    const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		+   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
	    const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		+   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;

	    // Accumulate acceleration terms
	    particles[v].ax += dax;
	    particles[v].ay += day;
	    particles[v].az += daz;

        }
    }

    if(*geo == 1){

      for (int j=0; j<N_real; j++){    

	particles[j].ax -= axe;
	particles[j].ay -= aye;
	particles[j].az -= aze;

      }

    }

}

/**
 * @brief Struct containing pointers to intermediate values
 */
struct reb_dpconst7 {
    double* const restrict p0;  ///< Temporary values at intermediate step 0 
    double* const restrict p1;  ///< Temporary values at intermediate step 1 
    double* const restrict p2;  ///< Temporary values at intermediate step 2 
    double* const restrict p3;  ///< Temporary values at intermediate step 3 
    double* const restrict p4;  ///< Temporary values at intermediate step 4 
    double* const restrict p5;  ///< Temporary values at intermediate step 5 
    double* const restrict p6;  ///< Temporary values at intermediate step 6 
};

static struct reb_dpconst7 dpcast(struct reb_dp7 dp){
    struct reb_dpconst7 dpc = {
        .p0 = dp.p0, 
        .p1 = dp.p1, 
        .p2 = dp.p2, 
        .p3 = dp.p3, 
        .p4 = dp.p4, 
        .p5 = dp.p5, 
        .p6 = dp.p6, 
    };
    return dpc;
}

typedef struct {
  double t, x, y, z, vx, vy, vz, ax, ay, az;
} tstate;

typedef struct {
    double* t;
    double* state;
    int n_out;
    int n_particles;
} timestate;

// Gauss Radau spacings
static const double h[9]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626, 1.0};

int integration_function(double tstart, double tstep, double trange,
			 int geocentric,
			 int n_particles,
			 double* instate,
			 timestate *ts){

    void store_function(struct reb_simulation* r, int n_out, int n_particles, tstate* last, double* outtime, double* outstate);
    struct reb_simulation* r = reb_create_simulation();

    // Set up simulation constants
    r->G = 0.295912208285591100E-03; // Gravitational constant (AU, solar masses, days)
    r->integrator = REB_INTEGRATOR_IAS15;
    r->heartbeat = NULL;
    r->display_data = NULL;
    r->collision = REB_COLLISION_NONE;  // This is important and needs to be considered carefully.
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity = REB_GRAVITY_NONE;
    
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

    int n_out = (int)8*fabs(trange/tstep);

    // Allocate arrays
    static double* outstate = NULL;
    static double* outtime = NULL;

    if(outstate == NULL){
	outstate = (double *) malloc(7*n_out*n_particles*6*sizeof(double));
    }else{
	outstate = (double *) realloc(outstate, 7*n_out*n_particles*6*sizeof(double));
    }

    if(outtime == NULL){    
	outtime  = (double *) malloc(n_out*sizeof(double));
    }else{
	outtime  = (double *) realloc(outtime, n_out*sizeof(double));	
    }

    rebx_set_param_int(rebx, &ephem_forces->ap, "n_out", 0);
    rebx_set_param_pointer(rebx, &ephem_forces->ap, "outstate", outstate);

    // Add and initialize particles    
    for(int i=0; i<n_particles; i++){

	struct reb_particle tp = {0};

	tp.x  =  instate[6*i+0];
	tp.y  =  instate[6*i+1];
	tp.z  =  instate[6*i+2];
	tp.vx =  instate[6*i+3];
	tp.vy =  instate[6*i+4];
	tp.vz =  instate[6*i+5];

	reb_add(r, tp);
    }

    // Add and initialize variational particles

    int var_i;
    for(int i=0; i<n_particles; i++){
 
        var_i = reb_add_var_1st_order(r, i); // The i corresponds to the index of the testparticle that we vary.
        r->particles[var_i].x = 1.; // perturbation of x-component test particles.
        r->particles[var_i].y = 0.; // lines superfluous; all var. particle ICs are zero by default.
        r->particles[var_i].z = 0.;
        r->particles[var_i].vx = 0.;
        r->particles[var_i].vy = 0.;
        r->particles[var_i].vz = 0.;
 
        var_i = reb_add_var_1st_order(r, i); 
        r->particles[var_i].x = 0.; 
        r->particles[var_i].y = 1.; 
        r->particles[var_i].z = 0.;
        r->particles[var_i].vx = 0.;
        r->particles[var_i].vy = 0.;
        r->particles[var_i].vz = 0.;
 
        var_i = reb_add_var_1st_order(r, i); 
        r->particles[var_i].x = 0.; 
        r->particles[var_i].y = 0.; 
        r->particles[var_i].z = 1.;
        r->particles[var_i].vx = 0.;
        r->particles[var_i].vy = 0.;
        r->particles[var_i].vz = 0.;
 
        var_i = reb_add_var_1st_order(r, i); 
        r->particles[var_i].x = 0.; 
        r->particles[var_i].y = 0.; 
        r->particles[var_i].z = 0.;
        r->particles[var_i].vx = 1.;
        r->particles[var_i].vy = 0.;
        r->particles[var_i].vz = 0.;
 
        var_i = reb_add_var_1st_order(r, i); 
        r->particles[var_i].x = 0.; 
        r->particles[var_i].y = 0.; 
        r->particles[var_i].z = 0.;
        r->particles[var_i].vx = 0.;
        r->particles[var_i].vy = 1.;
        r->particles[var_i].vz = 0.;
 
        var_i = reb_add_var_1st_order(r, i); 
        r->particles[var_i].x = 0.; 
        r->particles[var_i].y = 0.; 
        r->particles[var_i].z = 0.;
        r->particles[var_i].vx = 0.;
        r->particles[var_i].vy = 0.;
        r->particles[var_i].vz = 1.;

    }

    int N = r->N; // N includes real+variational particles

    r->t = tstart;    // set simulation internal time to the time of test particle initial conditions.
    r->dt = tstep;    // time step in days, this is just an initial value.

    outtime[0] = r->t;	
    
    for(int j=0; j<N; j++){
	int offset = j*6;
	outstate[offset+0] = r->particles[j].x;
	outstate[offset+1] = r->particles[j].y;
	outstate[offset+2] = r->particles[j].z;
	outstate[offset+3] = r->particles[j].vx;
	outstate[offset+4] = r->particles[j].vy;
	outstate[offset+5] = r->particles[j].vz;
    }

    tstate last[N]; 

    //reb_integrate(r, times[0]); // Not sure this is needed.
    reb_update_acceleration(r); // This is needed to save the acceleration.
 
    double tmax = tstart+trange;
    int i = 1;
    const double dtsign = copysign(1.,r->dt);   // Used to determine integration direction

    while((r->t)*dtsign<tmax*dtsign){ 

	fflush(stdout);
	// This could be a helper function
	// Save the previous completed state.
	for(int j=0; j<N; j++){ 
	    last[j].t = r->t;	
	    last[j].x = r->particles[j].x;
	    last[j].y = r->particles[j].y;
	    last[j].z = r->particles[j].z;
	    last[j].vx = r->particles[j].vx;
	    last[j].vy = r->particles[j].vy;
	    last[j].vz = r->particles[j].vz;
	    last[j].ax = r->particles[j].ax;
	    last[j].ay = r->particles[j].ay;
	    last[j].az = r->particles[j].az;
	}

	// Integrate a step.  
	reb_step(r);

	// Store the results, including the substeps
	store_function(r, 8*(i-1), N, last, outtime, outstate); //XYZ

	if((n_out - i*8) < 8){
	    n_out *= 2;
	    outstate = (double *) realloc(outstate, 7*n_out*n_particles*6*sizeof(double));
	    outtime  = (double *) realloc(outtime, n_out*sizeof(double));	    
	}

	reb_update_acceleration(r); // This is needed to save the acceleration.

	i++;

    }

    // Reallocate the arrays to trim off unneede space.
    //outstate = (double *) realloc(outstate, 7*8*(i-1)*n_particles*6*sizeof(double));
    //outtime  = (double *) realloc(outtime, 8*(i-1)*sizeof(double));	    

    // save pointers to the results in the timestate struct
    ts->t = outtime;
    ts->state = outstate;
    ts->n_particles = n_particles;
    ts->n_out = (i-1)*8;

    // explicitly free all the memory allocated by REBOUNDx

    //rebx_free(rebx);    
    //reb_free_simulation(r);

    fflush(stdout);
    return(1);
}

// This function is doing two related things:
// 1. Calculating the positions and velocities at the substeps
// 2. Storing the times and positions/velocities in the arrays
//    that are provided.
// For this to work, we need:
// * the last valid state for all particles,
// * the b coefficients for all the particles,
// * the last time step

void store_function(struct reb_simulation* r, int n_out, int n_particles, tstate* last, double* outtime, double* outstate){
    
    int N = r->N;
    int N3 = 3*N;
    double s[9]; // Summation coefficients

    // Store the time and state for each particle for last
    // completed step.
    outtime[n_out] = last[0].t;

    for(int j=0; j<n_particles; j++){
	// Rather than computing an offset, perhaps there should be
	// an index for the next open slot.
	int offset = (n_out*n_particles+j)*6;
	outstate[offset+0] = last[j].x;
	outstate[offset+1] = last[j].y;	
	outstate[offset+2] = last[j].z;
	outstate[offset+3] = last[j].vx;	
	outstate[offset+4] = last[j].vy;
	outstate[offset+5] = last[j].vz;	
    }
    
    // Convenience variable.  The 'br' field contains the 
    // set of coefficients from the last completed step.
    const struct reb_dpconst7 b  = dpcast(r->ri_ias15.br);

    double* x0 = malloc(sizeof(double)*N3);
    double* v0 = malloc(sizeof(double)*N3);
    double* a0 = malloc(sizeof(double)*N3);

    for(int j=0;j<N;j++) {

	const int k0 = 3*j+0;
	const int k1 = 3*j+1;
	const int k2 = 3*j+2;

	x0[k0] = last[j].x;
	x0[k1] = last[j].y;
	x0[k2] = last[j].z;	

	v0[k0] = last[j].vx;
	v0[k1] = last[j].vy;
	v0[k2] = last[j].vz;	

	a0[k0] = last[j].ax;
	a0[k1] = last[j].ay;
	a0[k2] = last[j].az;

    }

    // Loop over intervals using Gauss-Radau spacings      
    for(int n=1;n<8;n++) {                          

	s[0] = r->dt_last_done * h[n];

	s[1] = s[0] * s[0] / 2.;
	s[2] = s[1] * h[n] / 3.;
	s[3] = s[2] * h[n] / 2.;
	s[4] = 3. * s[3] * h[n] / 5.;
	s[5] = 2. * s[4] * h[n] / 3.;
	s[6] = 5. * s[5] * h[n] / 7.;
	s[7] = 3. * s[6] * h[n] / 4.;
	s[8] = 7. * s[7] * h[n] / 9.;

	double t = r->t + r->dt_last_done * (-1.0 + h[n]);
	outtime[n_out+n] = t;	

	// Predict positions at interval n using b values
	// for all the particles
	for(int j=0;j<N;j++) {  
	  //int mj = j;
	  const int k0 = 3*j+0;
	  const int k1 = 3*j+1;
	  const int k2 = 3*j+2;

	  double xx0 = x0[k0] + (s[8]*b.p6[k0] + s[7]*b.p5[k0] + s[6]*b.p4[k0] + s[5]*b.p3[k0] + s[4]*b.p2[k0] + s[3]*b.p1[k0] + s[2]*b.p0[k0] + s[1]*a0[k0] + s[0]*v0[k0] );
	  double xy0 = x0[k1] + (s[8]*b.p6[k1] + s[7]*b.p5[k1] + s[6]*b.p4[k1] + s[5]*b.p3[k1] + s[4]*b.p2[k1] + s[3]*b.p1[k1] + s[2]*b.p0[k1] + s[1]*a0[k1] + s[0]*v0[k1] );
	  double xz0 = x0[k2] + (s[8]*b.p6[k2] + s[7]*b.p5[k2] + s[6]*b.p4[k2] + s[5]*b.p3[k2] + s[4]*b.p2[k2] + s[3]*b.p1[k2] + s[2]*b.p0[k2] + s[1]*a0[k2] + s[0]*v0[k2] );

	  // Store the results
	  int offset = ((n_out+n)*n_particles+j)*6;
	  outstate[offset+0] = xx0;
	  outstate[offset+1] = xy0;	  	  
	  outstate[offset+2] = xz0;
	  
	}

	s[0] = r->dt_last_done * h[n];
	s[1] =      s[0] * h[n] / 2.;
	s[2] = 2. * s[1] * h[n] / 3.;
	s[3] = 3. * s[2] * h[n] / 4.;
	s[4] = 4. * s[3] * h[n] / 5.;
	s[5] = 5. * s[4] * h[n] / 6.;
	s[6] = 6. * s[5] * h[n] / 7.;
	s[7] = 7. * s[6] * h[n] / 8.;

	// Predict velocities at interval n using b values
	// for all the particles
	for(int j=0;j<N;j++) {

	  const int k0 = 3*j+0;
	  const int k1 = 3*j+1;
	  const int k2 = 3*j+2;

	  double vx0 = v0[k0] + s[7]*b.p6[k0] + s[6]*b.p5[k0] + s[5]*b.p4[k0] + s[4]*b.p3[k0] + s[3]*b.p2[k0] + s[2]*b.p1[k0] + s[1]*b.p0[k0] + s[0]*a0[k0];
	  double vy0 = v0[k1] + s[7]*b.p6[k1] + s[6]*b.p5[k1] + s[5]*b.p4[k1] + s[4]*b.p3[k1] + s[3]*b.p2[k1] + s[2]*b.p1[k1] + s[1]*b.p0[k1] + s[0]*a0[k1];
	  double vz0 = v0[k2] + s[7]*b.p6[k2] + s[6]*b.p5[k2] + s[5]*b.p4[k2] + s[4]*b.p3[k2] + s[3]*b.p2[k2] + s[2]*b.p1[k2] + s[1]*b.p0[k2] + s[0]*a0[k2];

	  // Store the results
	  int offset = ((n_out+n)*n_particles+j)*6;	  
	  outstate[offset+3] = vx0;
	  outstate[offset+4] = vy0;	  	  
	  outstate[offset+5] = vz0;

	}
    }

    free(x0);
    free(v0);
    free(a0);

}
