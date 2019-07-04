/******************************************************************
 * Simulate the 1+1d Thirring model using the fermion bag algorithm 
 * (arXiv:0910.5736) using a worldline representation. The mass
 * term is represented as a field of monomers (occupied sites) 
 * and the four fermion term is represented as dimers 
 * (occupied links). 
 *
 ******************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "worldline.h"

/* For saving a field before update with accept/reject */
int **old_field;
int **old_diraclink;

/* Saves the current configuration in a file */
void write_configuration(char * filename){
  FILE * config_file;

  int * buffer = malloc(2*NX*NT*sizeof(int));
  for (int t=0; t<NT; t++){
    for (int x=0; x<NX; x++){
      buffer[t*NX+x] = field[t][x];
      buffer[NX*NT+t*NX+x] = diraclink[t][x];
    }
  }
  
  config_file = fopen(filename,"wb");
  if (config_file){
    fwrite(buffer, 2*VOLUME, sizeof(int), config_file);
    fclose(config_file);
  } else {
    printf("Filename %s\n", filename);
    printf("Could not write configuration\n");
    exit(1);
  }
  free(buffer);
}

void read_configuration(char * filename){
  FILE * config_file;
  
  config_file = fopen(filename,"rb");
  int * buffer = malloc(2*NX*NT*sizeof(int));
  if (config_file){
    printf("Reading configuration\n");
    fread(buffer, 2*VOLUME, sizeof(int), config_file);
    fclose(config_file);

    for (int t=0; t<NT; t++){
      for (int x=0; x<NX; x++) {
        field[t][x] = buffer[t*NX+x];
        diraclink[t][x] = buffer[NT*NX+t*NX+x];
      }
    }
  } else {
    /* Create a cold configuration */
    printf("No configuration file, starting from a COLD state\n");
    for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
      field[t][x] = 0;
      diraclink[t][x] = 0;
    }

    // Thermalise
    printf("Thermalising\n");
    thermalise(NX);
  }
  free(buffer);
}



/* Make a copy of the fields */
void save_config(){
  for( int x=0; x<NX; x++ ) for( int t=0; t<NT; t++ ){
    old_field[t][x] = field[t][x];
    old_diraclink[t][x] = diraclink[t][x];
  }
}

/* Restore the saved field */
void restore_config(){
  for( int x=0; x<NX; x++ ) for( int t=0; t<NT; t++ ){
    field[t][x] = old_field[t][x];
    diraclink[t][x] = old_diraclink[t][x];
  }
}


/* Turn a link on */
static inline void link_on(int t, int x, int dir){
  int t2 = tdir(t, dir), x2 = xdir(x, dir);
  if ( field[t][x] == 0 && field[t2][x2] == 0 ){
    field[t][x] = dir+2;
    field[t2][x2] = opp_dir(dir)+2;
#ifdef DEBUG
    printf("Turned on link (%d,%d,%d) (%d,%d,%d)\n",t,x,dir,t2,x2,opp_dir(dir));
#endif
  } else {
    errormessage("Link already occupied\n");
  }
}

/* Add two monomers */
static inline void monomers_on(int t, int x, int dir){
  int t2 = tdir(t, dir), x2 = xdir(x, dir);
  if ( field[t][x] == 0 && field[t2][x2] == 0 ){
    field[t][x] = MONOMER;
    field[t2][x2] = MONOMER;
#ifdef DEBUG
    printf("Turned on monomers at (%d,%d) (%d,%d)\n",t,x,t2,x2);
#endif
  } else {
    errormessage("Sites already occupied\n");
  }
}

/* Turn a link off */
static inline void link_off(int t, int x, int dir){
  int t2 = tdir(t, dir);
  int x2 = xdir(x, dir);
  if ( field[t][x] != 0 && field[t2][x2] != 0 ){
    field[t][x] = 0;
    field[t2][x2] = 0;
#ifdef DEBUG
    printf("Turned off link (%d,%d) (%d,%d)\n",t,x,t2,x2);
#endif
  } else {
    errormessage("Link already off\n");
  }
}

/* Remove two monomers */
static inline void monomers_off(int t, int x, int dir){
  int t2 = tdir(t, dir);
  int x2 = xdir(x, dir);
  if ( field[t][x] == MONOMER && field[t2][x2] == MONOMER ){
    field[t][x] = 0;
    field[t2][x2] = 0;
#ifdef DEBUG
    printf("Turned off monomers at (%d,%d) (%d,%d)\n",t,x,t2,x2);
#endif
  } else {
    errormessage("Monomer already off\n");
  }
}




/* Try to add or remove a link at a given site */
int update_link_at(int s, int dir2)
{
  int success = 0;
  int t = s%NT, x = s/NT;

  if( field[t][x] >= LINK_TUP ) {
    int dir = field[t][x] - LINK_TUP;
    if( dir == dir2 ){
      int t2 = tdir(t,dir), x2 = xdir(x,dir);
      /* Remove link at t,x */
      //Note the factor of 4 from the dirac operator,
      //each dirac link has 0.5
      if( mersenne() < 1/(4*U) ) {
        link_off(t,x,dir);
        /* Replace with 2 opposing arrows */
        diraclink[t][x] = dir;
        diraclink[t2][x2] = opp_dir(dir);
        success = 1;
      }
    }
  } else if( field[t][x] == 0 ) {
    /* No link, add if possible */
    int dir = diraclink[t][x];
    if( dir == dir2 ){
      int t2 = tdir(t,dir), x2 = xdir(x,dir);
      int nu = diraclink[t2][x2];
      
      if( dir<NDIRS && dir == opp_dir(nu) ){
        //Two opposing arrows, easy to add
        //Note the factor of 4 from the dirac operator,
        //each dirac link has 0.5
        if( mersenne() < 4*U ) {
          link_on(t,x,dir);
          diraclink[t][x] = NDIRS;
          diraclink[t2][x2] = NDIRS;
          success = 1;
        }
      }
    }
  }
  return success;
}


int update_link()
{
  return update_link_at( (int) (mersenne()*VOLUME), (int) (mersenne()*ND) );
}


/* Try to add or remove monomers at a given site */
int update_monomers_at(int s, int dir)
{
  int success = 0;
  int t = s%NT, x = s/NT;

  if( field[t][x] == MONOMER ) {
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    if( field[t2][x2] == MONOMER ){
      /* Remove link at t,x */
      //Note the factor of 4 from the dirac operator,
      //each dirac link has 0.5
      if( mersenne() < 1/(4*m*m) ) {
        monomers_off(t,x,dir);
        /* Replace with 2 opposing arrows */
        diraclink[t][x] = dir;
        diraclink[t2][x2] = opp_dir(dir);
      }
    }
  } else if( field[t][x] == 0 ) {
    /* No link, add if possible */
    if( dir == diraclink[t][x] ){
      int t2 = tdir(t,dir), x2 = xdir(x,dir);
      int nu = diraclink[t2][x2];
      
      if( dir<NDIRS && dir == opp_dir(nu) ){
        //Two opposing arrows, easy to add
        //Note the factor of 4 from the dirac operator,
        //each dirac link has 0.5
        if( mersenne() < (4*m*m) ) {
          monomers_on(t,x,dir);
          diraclink[t][x] = NDIRS;
          diraclink[t2][x2] = NDIRS;
        }
      }
    }
  }
  return success;
}

int update_monomer()
{
  return update_monomers_at( (int) (mersenne()*VOLUME), (int) (mersenne()*ND)  );
}


/* Try to add a link at given coordinates */
int add_link_at(int t, int x)
{
  int success = 0;

  if( field[t][x] == 0 ) {

    int dir = diraclink[t][x];
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    int dir2 = diraclink[t2][x2];

    if( dir<NDIRS && dir2 == opp_dir(dir) ){
      //Two opposing arrows, easy to add
      //Note the factor of 4 from the dirac operator,
      //each dirac link has 0.5
      if( mersenne() < 4*U ) {
        link_on(t,x,dir);
        diraclink[t][x] = NDIRS;
        diraclink[t2][x2] = NDIRS;
        success = 1;
      }
    }
  }
  return success;
}


/* Try to remove a link at a given site */
int remove_link_at(int t, int x)
{
  int success = 0;

  if( field[t][x] >= LINK_TUP && field[t][x] <= LINK_XDN ) {
    /* Remove link at t,x */
    //Note the factor of 4 from the dirac operator,
    //each dirac link has 0.5
    int dir = field[t][x]-2;
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    if( mersenne() < 1./(4*U) ) {
      link_off(t,x,dir);
      /* Replace with 2 opposing arrows */
      diraclink[t][x] = dir;
      diraclink[t2][x2] = opp_dir(dir);
      success = 1;
    }
  }
  return success;
}


/* Get the sign of a given link */
int linksign( int t, int x, int dir ){
    int e = eta[t][x][dir%ND];
    if( dir == TDN ) e *= -1;
    if( dir == XDN ) e *= -1;
    if( t == NT-1 && dir == TUP ) e*= -1;
    if( t == 0 && dir == TDN ) e*= -1;
#ifdef ANTISYMMETRIC
    if( x == NX-1 && dir == XUP ) e*= -1;
    if( x == 0 && dir == XDN ) e*= -1;
#endif
    return e;
}

/* Calculate the sign of the propagator between the head and tail of a worm,
   the fermion propagator. Just follow the worm and count the signs. */
int propagator_sign(int t0, int x0, int t, int x){
  int dir = diraclink[t0][x0];
  int sign = linksign(t0,x0,dir);
  int t1 = tdir(t0, dir), x1= xdir(x0, dir);
  while( t1 != t || x1 != x ){
    dir = diraclink[t1][x1];
    sign *= linksign(t1,x1,dir);
    t1 = tdir(t1, dir), x1 = xdir(x1, dir);
  }
  return( sign );
}

/* Figure out the sign of the current
   configuration. In the most general
   case we go over all fermion loops. */
int negative_loops(){
  int sector = 0;

  // Mark the checked sites to avoid double counting
  int **checked = malloc( NT*sizeof(int *) );

  for(int t=0; t<NT; t++){
    checked[t] = malloc( NX*sizeof(int) );
    for(int x=0;x<NX;x++){
    checked[t][x] = 0;
    }
  }

  // Cycle trough all the sites and follow around each loop
  for(int t=0; t<NT; t++) for(int x=0;x<NX;x++) if(field[t][x]==0) if(checked[t][x] == 0){
    int dir, loop_sign = -1;
    int t1 = t, x1 = x;
    do {
      dir = diraclink[t1][x1];
      if( dir == NDIRS || dir == 10 ){
        loop_sign = 1;
        break;
      }
      loop_sign *= linksign(t1,x1,dir);
      checked[t1][x1] = 1;
      t1 = tdir(t1, dir), x1 = xdir(x1, dir);
    } while( t1 != t || x1 != x );

    if( loop_sign < 0 ){
      sector += 1;
    }
  }

  for(int t=0; t<NT; t++){
    free(checked[t]);
  }
  free(checked);
  
  return sector;
}

int configuration_sign(){
  int sign = 1;
  for(int t=0; t<NT; t++) for(int x=0;x<NX;x++) if(field[t][x]==0){
    int dir = diraclink[t][x];
    sign *= linksign(t,x,dir);
  }
  return sign;
}


/* Find a link pointing at a given site */
int find_link_pointing_at( int t, int x ){
  int t2,x2, dir;
  for( dir = 0; dir<NDIRS; dir++ ) {
    t2 = tdir(t, dir), x2 = xdir(x, dir);
    if( diraclink[t2][x2] == opp_dir(dir) ){
      //Found a site pointing to t,x
      break;
    }
  }
  return dir;
}


/* The propability of removing a link */
double remove_link_weight(int t, int x){
  double weight;
  int dir = diraclink[t][x];
  weight = 2;    //There is always a factor 0.5 for each link
  if( dir == TUP ) weight *= exp(-mu);
  if( dir == TDN ) weight *= exp(mu);
  return weight;
}

/* The propability of adding a link */
double add_link_weight(int t, int x, int dir){
  double weight;
  weight = 0.5;    //There is always a factor 0.5 for each link
  if( dir == TUP ) weight *= exp(mu);
  if( dir == TDN ) weight *= exp(-mu);
  return weight;
}

// A step in the worm update that updates monomers
void dirac_worm_add_monomer( int *t, int *x ){
  int t2, x2, dir;
  double p;
  dir = mersenne()*NDIRS;
  t2 = tdir(*t, dir), x2 = xdir(*x, dir);
  
  if( field[t2][x2] == 0 ){
    int removeddir = diraclink[t2][x2];
    if( removeddir == opp_dir(dir) ){
      p = m*remove_link_weight(t2, x2);
      if( mersenne() < p ){
        field[*t][*x] = MONOMER;
        diraclink[*t][*x] = NDIRS;
        diraclink[t2][x2] = 10;
        *t=t2; *x=x2;
      }
    }
  } else if( field[t2][x2] == MONOMER ){
    p = add_link_weight(*t, *x, dir) / m;
    if( mersenne() < p ){
      field[t2][x2] = 0;
      diraclink[*t][*x] = dir;
      diraclink[t2][x2] = 10;
      *t=t2; *x=x2;
    }
  }
}

//Update the Dirac background using a worm update
int update_dirac_background(){
  int t0, x0, t, x, dir;
  double p;
  int started = 0;

  //Pick a site
  t= mersenne()*NT, x=mersenne()*NX;
  
  if( field[t][x] == 0 ){
    //We break a link to create a fermion correlator.
    //The point the original link points to is the
    //starting point of the correlator

    p = remove_link_weight(t,x);
    if( mersenne() < p ){
      dir = diraclink[t][x];
      t0 = tdir(t, dir), x0 = xdir(x, dir);
      diraclink[t][x] = 10;
      started = 1;
    }
  } else if( field[t][x] == MONOMER ){
    // Selected a mass monomer. We just remove the monomer,
    // leaving an empty site. This is both the start and the
    // end point of the worm. Mass monomers allow local correlators.
    // The step with propability 0.5 corrects for the propability 
    // to try the reverse update
    if( mersenne() < 0.5 ){
      p = 1/m;
      if( mersenne() < p ){
        field[t][x] = 0;
        diraclink[t][x] = 10;
        t0 = t; x0 = x;
        started = 1;
      }
    }
  }
  
  if( started == 0 ){
    return 0;
  }
    
  for(int i=0;;i++) {
    //Also calculate the propagator between the two defects
    //This is the fermion propagator

    // First choose between adding propagating worm by
    // 1. Adding a Dirac link
    // 2. Adding/Removing a link
    // 3. Adding/Removing a monomer

    int choice = 3*mersenne();
    if( choice == 0 ){

      // Choose ending the worm or adding/removing monomers. The 0.5
      // propability here is corrected in the start worm update.
      if( mersenne() < 0.5 ){
        // Check if the start and end points overlap. This can
        // happen with mass monomers. If they do, try to close
        // the worm by adding a monomer.
        if( t0 == t &&  x0 == x ){
          // Match the other worm ending update, which picks a
          // random direction, propability 1/NDIRS
          if( mersenne() < 1.0/(NDIRS) ){
            p = m;
            if( mersenne() < p ){
              field[t][x] = MONOMER;
              diraclink[t][x] = NDIRS;
              break;
            }
          }
        }
      } else {
        // Add or remove a monomer
        dirac_worm_add_monomer( &t, &x );
      }

    } 
    if( choice == 1 ){
      // Add or remove a link
      // This is the only place in the worm where dimers
      // are updated
      int t2, x2, dir;
      dir = mersenne()*NDIRS;
      t2 = tdir(t, dir), x2 = xdir(x, dir);
      
      if( field[t2][x2] == 0 ){
        add_link_at( t2, x2 );
      } else {
        remove_link_at( t2, x2 ); 
      }

    } 
    if( choice == 2 ){
      // Pick a direction to propagate the worm
      int t2, x2;
  
      dir = mersenne()*NDIRS;
      t2 = tdir(t, dir), x2 = xdir(x, dir);

      // Check if a new link to this direction closes the worm
      // This needs to have the same weight as opening the worm
      if( t2 == t0 && x2 == x0 ) {
        // The weight for adding the link
        p = add_link_weight(t,x,dir);
        if( mersenne() < p ){
          diraclink[t][x] = dir;
          break;
        }
  
      } else {
      
        // Propagate the worm by adding a link
        if(field[t2][x2]==0) {
          // If the middle site is empty, we can point a link to it
          // There will then be a second link pointing at it, which we need
          // to remove. That will become the head of the worm
          int removeddir;
          int t3,x3;
          for( int dir2 = 0; dir2<NDIRS; dir2++ ) {
            t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
            if( diraclink[t3][x3] == opp_dir(dir2) ){
              //Found the site pointing to t2,x2
              removeddir = diraclink[t3][x3];
              break;
            }
          }
 
          //Calculate the propability and flip the link
          p = add_link_weight(t,x,dir);
          p *= remove_link_weight(t3,x3);
          if( mersenne() < p ){
            //flip the link
            diraclink[t][x] = dir;
            diraclink[t3][x3] = 10;
            t=t3; x=x3;
          }
        }
      }
    }
  }

  return 1;
}



/* A local plaquette update for the worldline configuration */
int plaquette_update(){
  int x1, t1, x2, t2, x3, t3, x4, t4, dir1, dir2, changes=0;
  // Pick a site
  t1= mersenne()*NT, x1=mersenne()*NX;

  if( field[t1][x1] == 0 ){
    dir1 = diraclink[t1][x1];
  
    // Pick orthogonal direction
    dir2 = (dir1+1)%ND;
    if(mersenne()<0.5){
      dir2 = opp_dir(dir2);
    }

    t2 = tdir(t1, dir1), x2 = xdir(x1, dir1);
    t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
    t4 = tdir(t1, dir2), x4 = xdir(x1, dir2);
    
    if(field[t4][x4] == MONOMER){
      // A monomer and a corner get flipped
      // >  v      v  o
      // o  >  ->  >  >
      if( diraclink[t2][x2] == dir2 ){
        diraclink[t1][x1] = dir2;
        field[t2][x2] = MONOMER;
        diraclink[t2][x2] = NDIRS;
        diraclink[t4][x4] = dir1;
        field[t4][x4] = 0;
        changes = 1;
      }
    } else {
      if( diraclink[t3][x3] == opp_dir(dir1) ){
        // Two links pointing to opposing directions are flipped
        // >  <      >  v
        // >  <  ->  ^  < 
        diraclink[t1][x1] = dir2;
        diraclink[t3][x3] = opp_dir(dir2);
        changes = 1;
      }
    }
  }
  return changes;
}

/* Try to flip the direction of a loop */
int flip_loop(){
  int x0, t0, t, x, t2, x2, dir, nt=0, changes=0;
  double weight;
  // Pick a site
  t0 = mersenne()*NT, x0 = mersenne()*NX;
  if( field[t0][x0] == 0 ){
    t = t0; x = x0;
    do {
      dir = diraclink[t][x];
      nt += dir == TUP;
      nt -= dir == TDN;
      t = tdir(t, dir), x = xdir(x, dir);
    } while( t0 != t || x0 != x );
  
    weight = exp(-nt*mu);
    if( mersenne() < weight ){
      t = t0; x = x0;
      dir = diraclink[t][x];
      do {
        int dir2;
        t = tdir(t, dir), x = xdir(x, dir);
        dir2 = dir;
        dir = diraclink[t][x];
        diraclink[t][x] = opp_dir(dir2);
      } while( t0 != t || x0 != x );
      changes = 1;
    }
  }
  return changes;
}


/* A number of thermalisation updates without accept/reject */
void thermalise( int nsteps ){
  for( int i=0; i<nsteps; i++ ){
    /* Worm update */
    update_dirac_background();
  }
}


/* A full update function. A single worm update followed by a number of random
   link and monomer updates */
int update_config( int nsteps )
{
  int changes=0;

  for( int i=0; i<nsteps; i++ ){
    /* Worm update */
    changes += update_dirac_background();
  }

  return changes;
}











void print_config()
{
  for (int t=NT; t--; ) {
    for (int x=0; x<NX; x++){
      int empty = 1;
      if(field[t][x]==MONOMER)  { empty = 0; printf(" o "); }
      if(field[t][x]==LINK_TUP) { empty = 0; printf(" | "); }
      if(field[t][x]==LINK_XUP) { empty = 0; printf(" --"); }
      if(field[t][x]==LINK_TDN) { empty = 0; printf(" | "); }
      if(field[t][x]==LINK_XDN) { empty = 0; printf("-- "); }
      if(field[t][x]==EMPTY) { empty = 0; printf("*"); }
      if(diraclink[t][x]==10) printf(" . ");
      if(diraclink[t][x]==TUP) printf(" ^ ");
      if(diraclink[t][x]==TDN) printf(" v ");
      if(diraclink[t][x]==XUP) printf(" > ");
      if(diraclink[t][x]==XDN) printf(" < ");
    }
    printf(" \n");
  }
  printf(" \n");
}



/* Measure the fermion and chiral charges */
void measure_charge(int *c, int *q){
  int _c = 0;
  int _q = 0;
  for (int x=0; x<NX; x++){
    if( diraclink[0][x] == TUP ){
      _c+= 1;
      _q+= ((x)%2==1?1:-1); 
    }
    if( diraclink[1][x] == TDN ){
      _c-= 1;
      _q+= ((x)%2==1?1:-1);
    }
    if( field[0][x] == LINK_TUP ){
      _q+= 2*((x)%2==1?1:-1);
    }
  }
  *c = _c;
  *q = _q;
}






//In the precense of (source) monomers, the sign
//can be negative. Assuming two monomers, find 
//the sign of the configuration
//Essentially count the surfaces crossed
//on a path between the monomers
int sign_with_monomers(int t0, int x0, int t1, int x1){
  int charge = 0;
  int ds, t = t0, x = x0;
  //First travel in t-direction
  if(t0<t1) {
    for(;t<t1-1;){
      t++;
      if(diraclink[t][x]==XUP) charge++;
      if(diraclink[t][xup[x]]==XDN) charge++;
    }
  }
  else {
    for(;t>t1;t--){
      if(diraclink[t][x]==XUP) charge++;
      if(diraclink[t][xup[x]]==XDN) charge++;
    }
  }

  //Then in the x direction
  if(x0<x1) {
    for(;x<x1-1;){
      x++;
      if(diraclink[t][x]==TUP) charge++;
      if(diraclink[tup[t]][x]==TDN) charge++;
    }
  }
  else {
    for(;x>x1;x--){
      if(diraclink[t][x]==TUP) charge++;
      if(diraclink[tup[t]][x]==TDN) charge++;
    }
  }
  return charge%2 == 0 ? 1:-1;
}




/* A worm for measuring the susceptibility in a wordline background
 * does not take the signs into account */

/* Propagate the worm by moving one of the sources to an unoccupied site.
 * This will always trigger a worldline worm update.
 * For detailed balance the worm can start in two ways:
 *   Moving the source and creating a defect in the background
 *   Just creating starting a worm update next to the source
 * Equivalently it can close in two way
 *   Hitting the source and moving it to close the worm
 *   Closing by hitting the defect
 */
int move_monomer(int monomer, int *ts0, int *xs0, int dir){

  int t = *ts0, x = *xs0;
  int t2 = tdir(t, dir), x2 = xdir(x, dir);
  int t0,x0;  //point where the propagator starts
  int ts,xs;  //the location of the other end of the worm (source)
  
  // Check correct use
  if( field[t][x] != monomer ){
    errormessage("Starting point of move_monomer is not a source monomer\n");
  }
  
  // This algorithm is only used when the neighbouring site
  // is unoccupied
  if( field[t2][x2] != 0 ){
    return 0;
  }
  
  // Start by finding the neighboring sites
  int done=1;
  int dir2 = find_link_pointing_at(t2,x2);
  int t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
  int dir3 = find_link_pointing_at(t3,x3);
  int t4 = tdir(t3, dir3), x4 = xdir(x3, dir3);


  // Branch between moving the source by two sites and keeping it at place
  // (the second one is required for detailed balance)
  if( mersenne() > 0.5) {

    double p = 2;
    if( dir == TUP ) p *= exp(mu);
    if( dir == TDN ) p *= exp(-mu);
    if( diraclink[t3][x3] == TUP ) p *= exp(-mu);
    if( diraclink[t3][x3] == TDN ) p *= exp(mu);
    if( diraclink[t4][x4] == TUP ) p *= exp(-mu);
    if( diraclink[t4][x4] == TDN ) p *= exp(mu);

    if( mersenne() < p ){
      // Remove the source monomer at the original point and move it to point 3
      field[t][x] = 0;
      field[t3][x3] = monomer;

      // This also means cutting the Worldline link at the site
      diraclink[t4][x4] = 10;
      diraclink[t3][x3] = NDIRS;
      diraclink[t][x] = dir; //this now points to site 2
      
      // Remember the new sites
      t0 = t; x0 = x;
      ts = t3; xs = x3;
      t=t4; x=x4;
      done = 0;
    }
    
  } else {

    // Don't move the source, just start a worm
    double p = 2;
    if( diraclink[t4][x4] == TUP ) p *= exp(-mu);
    if( diraclink[t4][x4] == TDN ) p *= exp(mu);

    if( mersenne() < p ){
      t0 = t3, x0 = x3;
      diraclink[t4][x4] = 10;
      ts = t; xs = x;
      t = t4; x = x4;
      done = 0;
    }
    
  }
  
  for(;done==0;){
    
      // Now keep updating the worm until it closes
      int dir = mersenne()*NDIRS;
      t2 = tdir(t, dir), x2 = xdir(x, dir);

      //If pointing at the source, move it back and close the worm
      if( t2==ts && x2 == xs ) {

        // Find the relevant directions
        int t3,x3;
        int dir2 = diraclink[t0][x0];
        t3 = tdir(t0, dir2), x3 = xdir(x0, dir2);
        for(dir2 = 0; dir2<NDIRS; dir2++){
            if( tdir(t2,dir2) == t3 && xdir(x2,dir2) == x3 ){
                break;
            }
        }
      
       // Calculate the propability of moving the source
       double p = 0.5;
       int dir0 = diraclink[t0][x0];
       if( dir == TUP ) p *= exp(mu);
       if( dir == TDN ) p *= exp(-mu);
       if( dir2 == TUP ) p *= exp(mu);
       if( dir2 == TDN ) p *= exp(-mu);
       if( diraclink[t0][x0] == TUP ) p *= exp(-mu);
       if( diraclink[t0][x0] == TDN ) p *= exp(mu);

       // Move and close the worm
       if( mersenne() < p) {
          field[t0][x0] = monomer;
          field[t2][x2] = 0;
          diraclink[t][x] = dir;
          diraclink[t0][x0] = NDIRS;
          diraclink[t2][x2] = dir2;
          *ts0 = t0; *xs0 = x0; 
          done = 1;
        }

      } else if(t2 == t0 && x2 == x0) {

        //If pointing at the defect, try to close the worm
        double p = 0.5;
        if( dir == TUP ) p *= exp(mu);
        if( dir == TDN ) p *= exp(-mu);
        if( mersenne() < p ){
          diraclink[t][x] = dir;
          *ts0 = ts; *xs0 = xs; 
          done = 1;
        }

      } else if(field[t2][x2]==0) {

        // This is the standard worm step
        // Find the directions and the propability of the move
        int removeddir;
        int t3,x3; 
        for( int dir2 = 0; dir2<NDIRS; dir2++ ) {
          t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
          if( diraclink[t3][x3] == opp_dir(dir2) ){
            break;
          }
        }
        double p=1;
        if( dir == TUP ) p *= exp(mu);
        if( dir == TDN ) p *= exp(-mu);
        if( diraclink[t3][x3] == TUP ) p *= exp(-mu);
        if( diraclink[t3][x3] == TDN ) p *= exp(mu);
        
        //flip and follow to a new site with propability p
        if( mersenne() < p ){
          //flip the link
          diraclink[t][x] = dir;
          diraclink[t3][x3] = 10;
          t=t3; x=x3;
        }
      }
  }
}



/* The actual worm for susceptibility. Does not take signs into account */
double measure_susceptibility_with_background( ){
  #ifdef OPENX
  int nsteps = 0;
  int n_attempts = 1;
  int n_links = 0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    if(field[t][x] >= LINK_TUP){
      n_links +=1;
    }
  }
  double scale_factor = n_links/(2.*VOLUME*n_attempts);
  
  for(int i=0; i<n_attempts;i++){
    //Pick a random link
    int t=0,x=0;
    if(n_links>0) do{
      t= mersenne()*NT, x=mersenne()*NX;
    } while(field[t][x] < LINK_TUP);
  
    if( field[t][x] >= LINK_TUP){
      //Replace the link with two source monomers
      int dir = field[t][x]-2;
      int t2 = tdir(t, dir), x2 = xdir(x, dir);

      field[t][x] = SOURCE_MONOMER;
      field[t2][x2] = SOURCE_MONOMER;
      n_links--;
    
      //Propagate the sources untill the worm closes
      for(;;){
        nsteps+=sign_with_monomers(t,x,t2,x2);
        int dir = mersenne()*NDIRS;
        int t3 = tdir(t2, dir), x3 = xdir(x2, dir);
      
        if( t3==t && x3==x ) {
           // Back at the original site, turn into a link 
           field[t][x]=0; field[t2][x2]=0;
           link_on(t2,x2,dir);
           n_links++;
           break;

         } else if( field[t3][x3] >= LINK_TUP ) {
           // found a link, flip it  
           int linkdir = field[t3][x3] - LINK_TUP;
           int t4 = tdir(t3,linkdir), x4 = xdir(x3,linkdir);

           field[t2][x2] = 0;
           link_off(t3,x3,linkdir); link_on(t2,x2,dir);
           field[t4][x4] = SOURCE_MONOMER;

           t2 = t4; x2 = x4;

         } else {
           // No link at found, run a worm that may move the
           // Source by two sites by updating the background field
           move_monomer_in_worldline(SOURCE_MONOMER, &t2, &x2, dir);
         }
      }
    }
  }

  return (double)nsteps*scale_factor;
  #else
  return NAN;
  #endif
}





void setup_lattice(long seed){
  /* "Warm up" the rng generator */
  seed_mersenne( seed );
  for (int i=0; i<543210; i++) mersenne();

  /* Allocate location and field arrays */
  field = malloc( NT*sizeof(int *) );
  diraclink = malloc( NT*sizeof(int *) );
  old_field = malloc( NT*sizeof(int *) );
  old_diraclink = malloc( NT*sizeof(int *) );
  eta = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++){
    field[t] = malloc( (NX+1)*sizeof(int) );
    diraclink[t] = malloc( (NX+1)*sizeof(int) );
    old_field[t] = malloc( (NX+1)*sizeof(int) );
    old_diraclink[t] = malloc( (NX+1)*sizeof(int) );
    eta[t] = malloc( (NX+1)*sizeof(int *) );
    for (int x=0; x<NX+1; x++) {
     eta[t][x] = malloc( 2*sizeof(int) );
    }
  }
  xup = malloc( (NX+1)*sizeof(int) );
  xdn = malloc( (NX+1)*sizeof(int) );
  tup = malloc( NT*sizeof(int) );
  tdn = malloc( NT*sizeof(int) );

  /* fill up the index array */
  for(int i=0; i<NT; i++) {
    tup[i] = (i+1) % NT;
    tdn[i] = (i-1+NT) % NT;
  }
  for(int i=0; i<NX+1; i++) {
    xup[i] = (i+1) % NX;
    xdn[i] = (i-1+NX) % NX;
  }
#ifdef OPENX  //Never match boundaries to actual sites
  xdn[0] = NX;
  xup[NX-1] = NX;
  xup[NX] = NX;
#endif

  /* fill the staggered eta array */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    eta[t][x][1] = 1;
    if( x%2 == 0 ){
      eta[t][x][0] = 1;
    } else {
      eta[t][x][0] = -1;
    }
#ifdef OPENX  
    eta[t][NX][0] = eta[t][NX][1] = 0;
#endif
  }

  int ** field_copy = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++) field_copy[t] = malloc( (NX+1)*sizeof(int) );

  #ifdef OPENX
  for (int t=0; t<NT; t++) {
    field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
    field_copy[t][NX] = EMPTY;
    diraclink[t][NX] = EMPTY;
  }
  #endif
}









