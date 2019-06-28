/****************************************************************
 * Simulate the 1+1D Thirring in the gauge representation using
 * a simple metropolis update
 *
 ****************************************************************/

#include "worldline.h"

/* Gauge field */
double ***A;

double gauge_action(){
    double S = 0;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
        S+= 1.-cos(A[t][x][dir]);
    }
    return U*S;
}



double complex determinant();
void update_gauge(){
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
      double s1 = -U*cos(A[t][x][dir]);
      double old = A[t][x][dir];
      A[t][x][dir] = 2*M_PI*mersenne()-M_PI;

      double s2 = -U*cos(A[t][x][dir]);
      double edS = exp(-s2+s1);
      
      if( mersenne() > edS )
        A[t][x][dir] = old;
    }
}








int main(void){
    
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  long seed;
  int n_loops, n_measure;
    
  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_loops);

  printf(" Updates / measurement : ");
  scanf("%d",&n_measure);
  
  printf(" Average over configurations : ");
  scanf("%d",&n_average);

  printf(" m : ");
  scanf("%lf",&m);

  printf(" g : ");
  scanf("%lf",&g);

  printf(" mu : ");
  scanf("%lf",&mu);

  printf(" Random number : ");
  scanf("%ld",&seed);
  seed_mersenne( seed );
  
  /* "Warm up" the rng generator */
  for ( int i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 2D quenched Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" g %f \n", g);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );
    
  eta = malloc( NT*sizeof(int *) );
  A = malloc( NT*sizeof(double *) );
  for (int t=0; t<NT; t++){
    A[t] = malloc( (NX+1)*sizeof(double *) );
    eta[t] = malloc( (NX+1)*sizeof(int *) );
    for (int x=0; x<NX+1; x++) {
      eta[t][x] = malloc( ND*sizeof(int) );
      A[t][x] = malloc( ND*sizeof(double) );
    }
  }
  xup = malloc( (NX+1)*sizeof(int) );
  xdn = malloc( (NX+1)*sizeof(int) );
  tup = malloc( NT*sizeof(int) );
  tdn = malloc( NT*sizeof(int) );
  
  /* fill up the index array */
  for ( int i=0; i<NT; i++ ) {
    tup[i] = (i+1) % NT;
    tdn[i] = (i-1+NT) % NT;
  }
  for (int i=0; i<NX+1; i++) {
    xup[i] = (i+1) % NX;
    xdn[i] = (i-1+NX) % NX;
  }

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    for( int dir=0; dir<ND; dir++ ) A[t][x][dir] = 0;
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

  for (int i=1; i<n_loops+1; i++) {
    update_gauge();
    
    if((i%n_measure)==0){
      /* Statistics */
      measure();
    }
  }
  
  
  return 0;
}















