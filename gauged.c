#define MAIN

#include "worldline.h"

/* Gauge field */
double ***A;
double g;

double gauge_action(){
    double S = 0;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
        S+= 1.-cos(A[t][x][dir]);
    }
    return S/(g*g);
}



void update_gauge(){
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
      double s1 = -cos(A[t][x][dir])/(g*g);
      double old = A[t][x][dir];
      A[t][x][dir] = 2*M_PI*mersenne()-M_PI;

      double s2 = -cos(A[t][x][dir])/(g*g);
      double edS = exp(-s2+s1);
      
      if( mersenne() > edS )
        A[t][x][dir] = old;
    }
}



double gauge_phase(){
  double phase = 0;
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) {
    int dir = diraclink[t][x];
    if( dir < ND ){
      phase += A[t][x][dir];
    } else {
      int tA = tdir(t,dir), xA = xdir(x,dir);
      dir = opp_dir(dir);
      phase -= A[tA][xA][dir];
    }
  }
  return phase;
}






int main(void){
    
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  int i,n_loops,n_measure,n_average;
  long seed;

  U=0; // Need a better way to disable dimers?
    
  /* Read in the input */
  get_int("Sites in the t direction", &NT);
  get_int("Sites in the x direction", &NX);
  get_int("Number of updates", &n_loops);
  get_int("Updates / measurement", &n_measure);
  get_int("Average over", &n_average);

  get_double("mass", &m);
  get_double("g", &g);
  get_double("mu", &mu);

  get_long("Random seed", &seed);
  get_char("Configuration filename ", configuration_filename);

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  //printf(" Git commit ID GIT_COMMIT_ID  \n");
  printf(" 2D Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" G %f \n", g);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );
  printf(" Configuration file %s\n", configuration_filename );

  
  // Set up lattice variables and fields
  setup_lattice(seed);
  read_configuration(configuration_filename);
    
  A = malloc( NT*sizeof(int **) );
  for (int t=0; t<NT; t++){
    A[t] = malloc( NX*sizeof(int*) );
    for (int x=0; x<NX+1; x++) {
     A[t][x] = malloc( 2*sizeof(int) );
    }
  }

  double gauge_action_sum = 0;

  for (int i=0; i<100; i++){
    update_gauge();
  }

  for (int i=1; i<n_loops+1; i++) {
    double phase;
    struct timeval start, end;
    double updatetime;

    for(int n=0; n<n_measure; n++){
      update_gauge();
      update_config(1);
    }

    gauge_action_sum += gauge_action();

    phase = gauge_phase();
    int sector = negative_loops();
    phase += M_PI*sector;
    printf("Sector %d\n", sector);
    printf("Phase %g\n", phase);
    
    if(i%n_average==0){
      printf("Gauge action %g \n", gauge_action_sum/n_average);
      gauge_action_sum = 0;
    }
  }
  
  
  return 0;
}




