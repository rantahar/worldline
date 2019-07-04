#ifndef WORLDLINE_H
#define WORLDLINE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <time.h>
#include <sys/time.h>

#include "lattice.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

/* Lattice size and dimensions */
EXTERN int NX;
EXTERN int NT;


#define VOLUME (NT*NX)

#define ANTISYMMETRIC //Antisymmetric boundaries
//#define SYMMETRIC   //implemented in Thirring_hop
//#define OPENX       //open in space, (anti)symmetric in time

#define SECTOR_STEP 0.01

/* Enumerate possible values for a field */
#define MONOMER 1
#define LINK_TUP (2+TUP)   // Links enumerated as 2+direction 
#define LINK_XUP (2+XUP)
#define LINK_TDN (2+TDN)
#define LINK_XDN (2+XDN)
#define SOURCE_MONOMER (2+NDIRS)
#define EMPTY -100   //A meta value for sites that don't exist

#define CG_ACCURACY 1e-30
#define CG_MAX_ITER 10000

/* Propability of exiting in the monomer moving worm update */
#define flip_exit_propability 0.2

#define WITH_MASS_MONOMERS
#define MAX_SECTOR 301





/* storage */
EXTERN int    ***eta;   //Staggered eta matrix
EXTERN double m;
EXTERN double U;
EXTERN double mu;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
EXTERN int **field;
EXTERN int **diraclink;
EXTERN int max_changes;
EXTERN char configuration_filename[100];






/* Worldline functions */
void setup_lattice(long seed);
int configuration_sign();
int update_config( int nsteps );
void thermalise( int nsteps );
void write_configuration(char * filename);
void read_configuration(char * filename);
void save_config();
void restore_config();
int negative_loops();

#endif
