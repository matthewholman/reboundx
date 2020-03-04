/** * @file central_force.c
 * @brief   A general central force.
 * @author  Dan Tamayo <tamayo.daniel@gmail.com>
 * 
 * @section     LICENSE
 * Copyright (c) 2015 Dan Tamayo, Hanno Rein
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.  All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Central Force$       // Effect category (must be the first non-blank line after dollar signs and between dollar signs to be detected by script).
 *
 * ======================= ===============================================
 * Authors                 D. Tamayo
 * Implementation Paper    *In progress*
 * Based on                None
 * C Example               :ref:`c_example_central_force`
 * Python Example          `CentralForce.ipynb <https://github.com/dtamayo/reboundx/blob/master/ipython_examples/CentralForce.ipynb>`_.
 * ======================= ===============================================
 * 
 * Adds a general central acceleration of the form a=Acentral*r^gammacentral, outward along the direction from a central particle to the body.
 * Effect is turned on by adding Acentral and gammacentral parameters to a particle, which will act as the central body for the effect,
 * and will act on all other particles.
 *
 * **Effect Parameters**
 * 
 * None
 *
 * **Particle Parameters**
 *
 * ============================ =========== ==================================================================
 * Field (C type)               Required    Description
 * ============================ =========== ==================================================================
 * Acentral (double)             Yes         Normalization for central acceleration.
 * gammacentral (double)         Yes         Power index for central acceleration.
 * ============================ =========== ==================================================================
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
void ephem(const double G, const int i, const double jde, double* const m,
	   double* const x, double* const y, double* const z,
	   double* const vx, double* const vy, double* const vz,
	   double* const ax, double* const ay, double* const az){

    static int initialized = 0;

    static struct _jpl_s *pl;
    struct mpos_s now;

    double M[11] =
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
      M[k] /= G;
    }

    if (initialized == 0){
      
      if ((pl = jpl_init()) == NULL) {
	fprintf(stderr, "could not load DE430 file, fool!\n");
	exit(EXIT_FAILURE);
      }

      printf("initialization complete\n");

      initialized = 1;

    }

    // Get position, velocity, and mass of body i in barycentric coords. 
    
    *m = M[i];

    jpl_calc(pl, &now, jde, ebody[i], PLAN_BAR); 

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

    if(0){ // make geocentric
      jpl_calc(pl, &now, jde, PLAN_EAR, PLAN_BAR); //earth in barycentric coords. 

      vecpos_div(now.u, pl->cau);
      vecpos_div(now.v, (pl->cau/86400.));
      vecpos_div(now.w, (pl->cau/86400.)*(pl->cau/86400.));
      
      *x -= now.u[0];
      *y -= now.u[1];
      *z -= now.u[2];
      *vx -= now.v[0];
      *vy -= now.v[1];
      *vz -= now.v[2];
      *ax -= now.w[0];
      *ay -= now.w[1];
      *az -= now.w[2];

    }
    
}

static void ast_ephem(const double G, const int i, const double jde, double* const m, double* const x, double* const y, double* const z){

    static int initialized = 0;

    static struct spk_s *spl;
    struct mpos_s pos;    

    if(i<0 || i>15){
      fprintf(stderr, "asteroid out of range\n");
      exit(EXIT_FAILURE);
    }

    // 1 Ceres, 4 Vesta, 2 Pallas, 10 Hygiea, 31 Euphrosyne, 704 Interamnia,
    // 511 Davida, 15 Eunomia, 3 Juno, 16 Psyche, 65 Cybele, 88 Thisbe, 
    // 48 Doris, 52 Europa, 451 Patientia, 87 Sylvia
    
    double M[16] =
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
      M[k] /= G;
    }

    if (initialized == 0){
      
      if ((spl = spk_init("sb431-n16s.bsp")) == NULL) {
	fprintf(stderr, "could not load sb431-n16 file, fool!\n");
	exit(EXIT_FAILURE);
      }
      printf("asteroid initialization complete\n");

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

    ephem(G, 3, t, &m, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze); // Get position and mass of earth

    // Calculate acceleration due to sun and planets
    for (int i=0; i<*N_ephem; i++){
        // Get position and mass of massive body i.
        ephem(G, i, t, &m, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az); 
        for (int j=0; j<N; j++){
  	  // Compute position vector of test particle j relative to massive body i.
            const double dx = xe + particles[j].x - x; 
            const double dy = ye + particles[j].y - y;
            const double dz = ze + particles[j].z - z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double prefac = G*m/(_r*_r*_r);
            particles[j].ax -= prefac*dx;
            particles[j].ay -= prefac*dy;
            particles[j].az -= prefac*dz;
        }
    }

    // Get position, velocity, and mass of the sun wrt barycenter
    // in order to translate heliocentric asteroids to barycenter.
    ephem(G, 0, t, &m, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs); 

    // Now calculate acceleration due to massive asteroids
    for (int i=0; i<*N_ast; i++){
        ast_ephem(G, i, t, &m, &x, &y, &z); // Get position and mass of asteroid i.
	x += xs;
	y += ys;
	z += zs;
        for (int j=0; j<N; j++){
  	  // Compute position vector of test particle j relative to massive body i.
            const double dx = xe + particles[j].x - x; 
            const double dy = ye + particles[j].y - y;
            const double dz = ze + particles[j].z - z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double prefac = G*m/(_r*_r*_r);
            particles[j].ax -= prefac*dx;
            particles[j].ay -= prefac*dy;
            particles[j].az -= prefac*dz;
        }
    }

    // Here is the GR treatment
    const double Msun = 1.0; // mass of sun in solar masses.
    const double mu = G*Msun; // careful here.  We are assuming that the central body is at the barycenter.
    const int max_iterations = 10; // careful of hard-coded parameter.
    for (int i=0; i<N; i++){
        struct reb_particle p = particles[i];
        struct reb_vec3d vi;

	p.x += xe;
	p.y += ye;
	p.z += ze;	
	p.vx += vxe;
	p.vy += vye;
	p.vz += vze;	
	
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

    // These are in here for another purpose.
    /*
    double ae = sqrt(axe*axe + aye*aye + aze*aze);
    double rho = sqrt(G*Msun/ae);
    */

    for (int i=0; i<N; i++){    

      particles[i].ax -= axe;
      particles[i].ay -= aye;
      particles[i].az -= aze;

    }
    //printf("\n");

}
