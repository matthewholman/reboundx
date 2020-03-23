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
void ephem(const double G, const int i, const double jde, double* const m,
	   double* const x, double* const y, double* const z,
	   double* const vx, double* const vy, double* const vz,
	   double* const ax, double* const ay, double* const az){

    static int initialized = 0;

    static struct _jpl_s *pl;
    struct mpos_s now;

    static double M[11];

    if(i<0 || i>10){
      fprintf(stderr, "body out of range\n");
      exit(EXIT_FAILURE);
    }

    if (initialized == 0){
      
      if ((pl = jpl_init()) == NULL) {
	fprintf(stderr, "could not load DE430 file, fool!\n");
	exit(EXIT_FAILURE);
      }

      // The values below are G*mass, so we need to divide by G.
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

      for(int k=0; k<11; k++){
	M[k] = JPL_GM[k]/G;
      }
      
      initialized = 1;

    }

    // Get position, velocity, and mass of body i in barycentric coords. 
    
    *m = M[i];

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

static void ast_ephem(const double G, const int i, const double jde, double* const m, double* const x, double* const y, double* const z){

    static int initialized = 0;

    static struct spk_s *spl;
    struct mpos_s pos;

    static double M[11];    

    if(i<0 || i>15){
      fprintf(stderr, "asteroid out of range\n");
      exit(EXIT_FAILURE);
    }

    if (initialized == 0){
      
      if ((spl = spk_init("sb431-n16s.bsp")) == NULL) {
	fprintf(stderr, "could not load sb431-n16 file, fool!\n");
	exit(EXIT_FAILURE);
      }

      // 1 Ceres, 4 Vesta, 2 Pallas, 10 Hygiea, 31 Euphrosyne, 704 Interamnia,
      // 511 Davida, 15 Eunomia, 3 Juno, 16 Psyche, 65 Cybele, 88 Thisbe, 
      // 48 Doris, 52 Europa, 451 Patientia, 87 Sylvia

      // The values below are G*mass, so we need to divide by G.
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
      
      for(int k=0; k<16; k++){
	M[k] = JPL_GM[i]/ G;
      }
      
      initialized = 1;

    }

    *m = M[i];
    spk_calc(spl, i, jde, &pos);          
    *x = pos.u[0];
    *y = pos.u[1];
    *z = pos.u[2];
    
}

void rebx_ephemeris_forces(struct reb_simulation* const sim, struct rebx_force* const force, struct reb_particle* const particles, const int N){
    const int* const N_ephem = rebx_get_param(sim->extras, force->ap, "N_ephem");
    const int* const N_ast = rebx_get_param(sim->extras, force->ap, "N_ast");
    if (N_ephem == NULL){
        fprintf(stderr, "REBOUNDx Error: Need to set N_ephem for ephemeris_forces\n");
        return;
    }

    const double G = sim->G;
    const double t = sim->t;

    double* c = rebx_get_param(sim->extras, force->ap, "c");
    if (c == NULL){
        reb_error(sim, "REBOUNDx Error: Need to set speed of light in gr effect.  See examples in documentation.\n");
        return;
    }
    const double C2 = (*c)*(*c);

    double m, x, y, z, vx, vy, vz, ax, ay, az;
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
    double xe, ye, ze, vxe, vye, vze, axe, aye, aze;
    double xo, yo, zo, vxo, vyo, vzo;
    double xr, yr, zr, vxr, vyr, vzr;

    int* geo = rebx_get_param(sim->extras, force->ap, "geocentric");

    // Get mass, position, velocity, and acceleration of the Earth and Sun
    // for later use
    ephem(G, 3, t, &m, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze);
    ephem(G, 0, t, &m, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);     

    // The offset position is used to adjust the particle positions.
    if(*geo == 1){
      xo = xe;    yo = ye;    zo = ze;
      vxo = vxe;  vyo = vye;  vzo = vze;
    }else{
      xo = 0.0;  yo = 0.0;  zo = 0.0;
      vxo = 0.0; vyo = 0.0; vzo = 0.0;      
    }
     
    // Calculate acceleration due to sun and planets
    for (int i=0; i<*N_ephem; i++){
        // Get position and mass of massive body i.
        ephem(G, i, t, &m, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az); 
        for (int j=0; j<N; j++){
  	  // Compute position vector of test particle j relative to massive body i.
	  const double dx =  particles[j].x + (xo - x); 
	  const double dy =  particles[j].y + (yo - y);
	  const double dz =  particles[j].z + (zo - z);
	  const double _r = sqrt(dx*dx + dy*dy + dz*dz);
	  const double prefac = G*m/(_r*_r*_r);

	  //printf("%le %le %le\n", dx, dy, dz);
	  
	  particles[j].ax -= prefac*dx;
	  particles[j].ay -= prefac*dy;
	  particles[j].az -= prefac*dz;
        }
    }


    // Calculate acceleration due to massive asteroids
    for (int i=0; i<*N_ast; i++){

        ast_ephem(G, i, t, &m, &x, &y, &z); // Get position and mass of asteroid i.

	// Translate massive asteroids from heliocentric to barycenter.
	x += xs; y += ys; z += zs;
	
        for (int j=0; j<N; j++){
  	  // Compute position vector of test particle j relative to massive body i.
	    const double dx = particles[j].x + (xo - x);
	    const double dy = particles[j].y + (yo - y);
	    const double dz = particles[j].z + (zo - z); 	    	    
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
	    
            const double prefac = G*m/(_r*_r*_r);
            particles[j].ax -= prefac*dx;
            particles[j].ay -= prefac*dy;
            particles[j].az -= prefac*dz;
        }
    }

    // Here is the treatment of the Earth's J2 and J4.
    // Borrowed code from gravitational_harmonics.
    // Assuming the coordinates are geocentric.
    // Also assuming that Earth's pole is along the z
    // axis.  This is only precisely true at the J2000
    // epoch.
    //

    // The geocenter is the reference for the J2/J4 calculations.
    xr = xe;  yr = ye;  zr = ze;

    // Hard-coded constants.  BEWARE!
    // Clean up on aisle 3!
    const double Mearth = 0.888769244512563400E-09/G;
    const double J2e = 0.00108262545*1.001;
    const double J4e = -0.000001616;
    const double au = 149597870.700;
    const double Re_eq = 6378.1263/au;

    // Unit vector to equatorial pole at the epoch

    double RAs =  359.87123273*M_PI/180.;
    double Decs =  89.88809752*M_PI/180.;

    //double xp = cos(Decs)*cos(RAs);
    //double yp = cos(Decs)*sin(RAs);
    //double zp = sin(Decs);

    //printf("%lf %lf %lf\n", xp, yp, zp);

    double xp =  0.0019111736356920146;
    double yp = -1.2513100974355823e-05;
    double zp =   0.9999981736277104;

    //xp = 0.0;
    //yp = 0.0;
    //zp = 1.0;
    
    double incl = acos(zp);
    double longnode;
    if(xp != 0.0 || yp !=0.0) {    
      longnode = atan2(xp, -yp);
    } else {
      longnode = 0.0;
    }

    for (int i=0; i<N; i++){
        const struct reb_particle p = particles[i];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
	
	// Rotate to Earth equatorial frame

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

	// Calculate acceleration in body frame
	
        const double costheta2 = dz*dz/r2;
        const double J2e_prefac = 3.*J2e*Re_eq*Re_eq/r2/r2/r/2.;
        const double J2e_fac = 5.*costheta2-1.;

	double resx = G*Mearth*J2e_prefac*J2e_fac*dx;
	double resy = G*Mearth*J2e_prefac*J2e_fac*dy;
	double resz = G*Mearth*J2e_prefac*(J2e_fac-2.)*dz;	

        const double J4e_prefac = 5.*J4e*Re_eq*Re_eq*Re_eq*Re_eq/r2/r2/r2/r/8.;
        const double J4e_fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;

        resx += G*Mearth*J4e_prefac*J4e_fac*dx;
        resy += G*Mearth*J4e_prefac*J4e_fac*dy;
        resz += G*Mearth*J4e_prefac*(J4e_fac+12.-28.*costheta2)*dz;

	// Rotate back to original frame
	// Rotate around x by -Dec
	double resxp =  resx;
	double resyp =  resy * cosd + resz * sind;
	double reszp = -resy * sind + resz * cosd;
	
	// Rotate around z by -RA
	resx =  resxp * cosr + resyp * sinr;
	resy = -resxp * sinr + resyp * cosr;
	resz =  reszp;

	particles[i].ax += resx;
        particles[i].ay += resy; 
        particles[i].az += resz;
	
	
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
    
    for (int i=0; i<N; i++){
        const struct reb_particle p = particles[i];
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

	// Calculate acceleration
	double resx = G*Msun*J2s_prefac*J2s_fac*dx;
	double resy = G*Msun*J2s_prefac*J2s_fac*dy;
	double resz = G*Msun*J2s_prefac*(J2s_fac-2.)*dz;

	// Rotate back to original frame
	// Rotate around x by -Dec
	double resxp =  resx;
	double resyp =  resy * cosd + resz * sind;
	double reszp = -resy * sind + resz * cosd;
	
	// Rotate around z by -RA
	resx =  resxp * cosr + resyp * sinr;
	resy = -resxp * sinr + resyp * cosr;
	resz =  reszp;

        particles[i].ax += resx;
        particles[i].ay += resy;
        particles[i].az += resz;
	
    }

    // Here is the Solar GR treatment
    // The Sun is the reference for these calculations.    
    xr  = xs;  yr  = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;

    //xr = 0.0;  yr = 0.0;  zr = 0.0;
    //vxr = 0.0; vyr = 0.0; vzr = 0.0;      

    const double mu = G*Msun; 
    const int max_iterations = 10; // hard-coded parameter.
    for (int i=0; i<N; i++){
        struct reb_particle p = particles[i];
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
	
        particles[i].ax += B*(1.-A)*p.x - A*p.ax - D*vi.x;
        particles[i].ay += B*(1.-A)*p.y - A*p.ay - D*vi.y;
        particles[i].az += B*(1.-A)*p.z - A*p.az - D*vi.z;

    }

    // The expressions below are in here for another purpose.
    /*
    double ae = sqrt(axe*axe + aye*aye + aze*aze);
    double rho = sqrt(G*Msun/ae);
    */

    if(*geo == 1){

      for (int i=0; i<N; i++){    

	particles[i].ax -= axe;
	particles[i].ay -= aye;
	particles[i].az -= aze;

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

// Gauss Radau spacings
static const double h[9]    = { 0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626, 1.0};

int integration_function(double tstart, double tstep, double trange, int geocentric,
			  double xi, double yi, double zi, double vxi, double vyi, double vzi,
			  tstate *outstate, int* n_out){

    void store_function(struct reb_simulation* r, tstate* outstate, tstate last, int n_out);
    struct reb_simulation* r = reb_create_simulation();

    // Set up simulation constants
    r->G = 0.295912208285591100E-03; // Gravitational constant (AU, solar masses, days)
    r->integrator = REB_INTEGRATOR_IAS15;
    r->heartbeat = NULL;
    r->display_data = NULL;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->gravity = REB_GRAVITY_NONE;
    //r->usleep = 20000.;
    
    int nsteps = fabs(trange/tstep);
    int nouts = 8*fabs(trange/tstep);

    double* times = malloc(nouts*sizeof(double));

    for(int i=0; i<nsteps; i++){
      times[i] = tstart + i*tstep;
    }

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

    struct reb_particle tp = {0};

    tp.x  =  xi;
    tp.y  =  yi;
    tp.z  =  zi;
    tp.vx =  vxi;
    tp.vy =  vyi;
    tp.vz =  vzi;
    
    reb_add(r, tp);

    r->t = tstart;    // set simulation internal time to the time of test particle initial conditions.
    //tmax  = r->t + trange;
    r->dt = tstep;    // time step in days

    outstate[0].t = r->t;	
    outstate[0].x = r->particles[0].x;
    outstate[0].y = r->particles[0].y;
    outstate[0].z = r->particles[0].z;
    outstate[0].vx = r->particles[0].vx;
    outstate[0].vy = r->particles[0].vy;
    outstate[0].vz = r->particles[0].vz;

    tstate last;

    reb_integrate(r, times[0]);    
    reb_update_acceleration(r); // This should not be needed but is.
    for(int j=1; j<nsteps; j++){

      last.t = r->t;	
      last.x = r->particles[0].x;
      last.y = r->particles[0].y;
      last.z = r->particles[0].z;
      last.vx = r->particles[0].vx;
      last.vy = r->particles[0].vy;
      last.vz = r->particles[0].vz;
      last.ax = r->particles[0].ax;
      last.ay = r->particles[0].ay;
      last.az = r->particles[0].az;

      reb_integrate(r, times[j]);
      //reb_step(r);
      store_function(r, outstate, last, 8*(j-1));
      reb_update_acceleration(r);

    }
    
    *n_out = (nsteps-1)*8;

    free(times);

    return(1);
}

void store_function(struct reb_simulation* r, tstate* outstate, tstate last, int n_out){

      int N = r->N;
      int N3 = 3*N;
      double s[9]; // Summation coefficients

      outstate[n_out].t = last.t;
      outstate[n_out].x = last.x;
      outstate[n_out].y = last.y;
      outstate[n_out].z = last.z;
      outstate[n_out].vx = last.vx;
      outstate[n_out].vy = last.vy;
      outstate[n_out].vz = last.vz;

      // Convenience variable.  The 'br' field contains the 
      // set of coefficients from the last completed step.
      const struct reb_dpconst7 b  = dpcast(r->ri_ias15.br);

      double* x0 = malloc(sizeof(double)*N3);
      double* v0 = malloc(sizeof(double)*N3);
      double* a0 = malloc(sizeof(double)*N3);

      for(int i=0;i<N;i++) {

	const int k0 = 3*i+0;
	const int k1 = 3*i+1;
	const int k2 = 3*i+2;

	x0[k0] = last.x;
	x0[k1] = last.y;
	x0[k2] = last.z;	

	v0[k0] = last.vx;
	v0[k1] = last.vy;
	v0[k2] = last.vz;	

	a0[k0] = last.ax;
	a0[k1] = last.ay;
	a0[k2] = last.az;

      }

      // Loop over interval using Gauss-Radau spacings      
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

	// Predict positions at interval n using b values	
	for(int i=0;i<N;i++) {  
	  //int mi = i;
	  const int k0 = 3*i+0;
	  const int k1 = 3*i+1;
	  const int k2 = 3*i+2;

	  double xx0 = x0[k0] + (s[8]*b.p6[k0] + s[7]*b.p5[k0] + s[6]*b.p4[k0] + s[5]*b.p3[k0] + s[4]*b.p2[k0] + s[3]*b.p1[k0] + s[2]*b.p0[k0] + s[1]*a0[k0] + s[0]*v0[k0] );
	  double xy0 = x0[k1] + (s[8]*b.p6[k1] + s[7]*b.p5[k1] + s[6]*b.p4[k1] + s[5]*b.p3[k1] + s[4]*b.p2[k1] + s[3]*b.p1[k1] + s[2]*b.p0[k1] + s[1]*a0[k1] + s[0]*v0[k1] );
	  double xz0 = x0[k2] + (s[8]*b.p6[k2] + s[7]*b.p5[k2] + s[6]*b.p4[k2] + s[5]*b.p3[k2] + s[4]*b.p2[k2] + s[3]*b.p1[k2] + s[2]*b.p0[k2] + s[1]*a0[k2] + s[0]*v0[k2] );

	  double t = r->t + r->dt_last_done * (h[n] - 1.0);
	  
	  outstate[n_out+n].t = t;
	  outstate[n_out+n].x = xx0;
	  outstate[n_out+n].y = xy0;
	  outstate[n_out+n].z = xz0;
	  
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
	for(int i=0;i<N;i++) {

	  const int k0 = 3*i+0;
	  const int k1 = 3*i+1;
	  const int k2 = 3*i+2;

	  double vx0 = v0[k0] + s[7]*b.p6[k0] + s[6]*b.p5[k0] + s[5]*b.p4[k0] + s[4]*b.p3[k0] + s[3]*b.p2[k0] + s[2]*b.p1[k0] + s[1]*b.p0[k0] + s[0]*a0[k0];
	  double vy0 = v0[k1] + s[7]*b.p6[k1] + s[6]*b.p5[k1] + s[5]*b.p4[k1] + s[4]*b.p3[k1] + s[3]*b.p2[k1] + s[2]*b.p1[k1] + s[1]*b.p0[k1] + s[0]*a0[k1];
	  double vz0 = v0[k2] + s[7]*b.p6[k2] + s[6]*b.p5[k2] + s[5]*b.p4[k2] + s[4]*b.p3[k2] + s[3]*b.p2[k2] + s[2]*b.p1[k2] + s[1]*b.p0[k2] + s[0]*a0[k2];

	  outstate[n_out+n].vx = vx0;
	  outstate[n_out+n].vy = vy0;
	  outstate[n_out+n].vz = vz0;
	  
	}
	
      }

      free(x0);
      free(v0);
      free(a0);

}
