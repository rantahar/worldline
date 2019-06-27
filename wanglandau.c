

#define MAIN

#include "worldline.h"
#include "SAD.h"

double active_sector_free_energy;


void init_sector_weights( double Weights[MAX_SECTOR], int max_init_steps ){
  for( int i=0; i<MAX_SECTOR; i++){
    Weights[i]=0.01;
  }
  for( int i=0; i<max_init_steps; i++){
    update(1);
    Weights[current_sector] += 1;
    if( Weights[current_sector] > 1000 ){
      // Effective stepsize is reduced to 0.001
      break;
    }
  }
}

void init_free_energy( int max_init_steps ){
  double Weights[MAX_SECTOR];

  printf(" Initialising free energy by direct measurement\n" );
  for( int i=0; i<MAX_SECTOR; i++){
    WangLaundau_F[i]=0;
  }

  /* Run a couple of Newton steps */
  for( int i=0; i<2; i++ ){
    init_sector_weights( Weights, max_init_steps );
    double average, sum = 0;
    for( int s=0; s<MAX_SECTOR; s++){
      sum += Weights[s];
    }
    average = sum/MAX_SECTOR;
    for( int s=0; s<MAX_SECTOR; s++){
      double logw = log(Weights[s]/average);
      logw = fmax( logw, -2);
      WangLaundau_F[s] += logw;
    }
  }
  for( int s=0; s<MAX_SECTOR; s++){
    printf("INIT SECTOR %d %g %g \n", s, WangLaundau_F[s], Weights[s]);
  }
}

void WangLaundau_setup( int max_init_steps ){
  FILE *config_file;
  config_file = fopen(init_parameter_filename, "rb");
  int initialized = 0;
  
  if(config_file) {
    printf(" Reading initial free energy\n" );
    for( int s=0; s<MAX_SECTOR; s++){
      fscanf(config_file, "%lf %ld %d\n", &WangLaundau_F[s], &WangLaundau_iteration[s], &WL_measure_sector[s]);
      if( WangLaundau_iteration[s] > 0 ){
        initialized = 1;
      }
    }
    fscanf(config_file, "%ld %d\n", &WL_nstep, &last_range_update);
    fclose(config_file);
  }

  if( initialized == 0 ) {
    for( int s=0; s<MAX_SECTOR; s++){
      WangLaundau_iteration[s] = 0;
      WL_measure_sector[s] = 0;
    }
    //init_free_energy( max_init_steps );
  }

  current_sector = count_negative_loops();
}

void WangLaundau_write_energy(){
  FILE * config_file;

  config_file = fopen(init_parameter_filename,"wb");
  if (config_file){
    for( int s=0; s<MAX_SECTOR; s++){
      fprintf(config_file, "%g %ld %d\n", WangLaundau_F[s], WangLaundau_iteration[s], WL_measure_sector[s]);
    }
    fprintf(config_file, "%d %d\n", WL_nstep, last_range_update);
    fclose(config_file);
  } else {
    printf("Could not write Wang Landau energy\n");
    exit(1);
  }
}

// Update the free energy in the Wang-Landau algorithm
void WangLaundau_update(sector){
#ifdef SAD
  int max_sector, max_histogram = 0;
  double Tmin = 100;
  double Ns;
  double step;
  double tolerance = -40;
  double maximum = -10000;

  if( WL_nstep==1 ){ // first call
    step = 1;
    WL_measure_sector[sector] = 1;
    WangLaundau_F[sector] += step;
    WangLaundau_iteration[sector]++;
    WL_nstep++;
    return;
  }

  for( int s=0; s<MAX_SECTOR; s++){
    if( WangLaundau_iteration[s] > max_histogram ){
      max_histogram = WangLaundau_iteration[s];
      max_sector = s;
    }
  }
  
  if( !WL_measure_sector[max_sector] ){
    WL_measure_sector[max_sector] = 1;
    last_range_update = WL_nstep;
    printf("New sector %d %d\n", max_sector, WL_nstep);
  }

  WangLaundau_iteration[sector]++;
  
  if( WL_measure_sector[sector] ){
    double Ns = 1;
    for( int s=0; s<MAX_SECTOR; s++ ){
      Ns += WL_measure_sector[sector];
    }
    double t0 = Ns/Tmin;
    double st2 = (double) WL_nstep * (double) WL_nstep;
    double l = last_range_update;
    double e = t0 + WL_nstep/l;
    double d = t0 + st2/(Ns*l);
    step = llr_alpha*e/d;
    //printf("SAD %g %d %d %g %g\n", step, last_range_update, WL_nstep, e, d);
    WangLaundau_F[sector] += step;
    if( step <= 0 ){
      printf("Negative step in SAD %g %g %g\n",step, e, d);
      exit(1);
    }

    for( int s=0; s<MAX_SECTOR; s++)
      if( WangLaundau_F[s] > maximum )
        maximum = WangLaundau_F[s];
    for( int s=0; s<MAX_SECTOR; s++){
      WangLaundau_F[s] -= maximum;
      if( WangLaundau_F[s] < tolerance ){
        WangLaundau_F[s] = tolerance;
      }
    }
  }
  WL_nstep++;

#else
  if( WL_measure_sector[sector] ){
    double step = llr_alpha/(WangLaundau_iteration[sector]+llr_constant_steps);
    WangLaundau_F[sector] += step;
    WangLaundau_iteration[sector]++;
  }
#endif
}

double WangLaundau_weight(new_sector,old_sector){
  return exp(WangLaundau_F[old_sector]-WangLaundau_F[new_sector]);
}


int llr_accept(){
  int accept = 1;
  int sector;
  double weight;
  sector = count_negative_loops();
  if( sector != current_sector ){
    weight = WangLaundau_weight(sector,current_sector);
    if( mersenne() < weight ){
      sector_changes += 1;
      accept = 1;
      current_sector = sector;
    } else {
      accept = 0;
    }
  }
  if( accept ){
    llr_accepted += 1;
  }
  return accept;
}

int negative_loops(){
  return current_sector;
}



/* A full update function. A single worm update followed by a number of random
   link and monomer updates */
int update( int nsteps )
{
  int changes=0;
  save_config();

  changes += update_config(nsteps);

  if( ! llr_accept() ){
    restore_config();
  }

  return changes;
}


/* Main function
 */
int main(int argc, char* argv[])
{
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  current_sector  = 0;
  llr_accepted =0;
  sector_changes =0;
  WL_nstep =1;
  last_range_update =1;

  int i,n_loops,n_measure,n_average,llr_update_every;
  long seed;

  /* Read in the input */
  get_int("Number of updates", &n_loops);
  get_int("Updates / measurement", &n_measure);
  get_int("Average over", &n_average);

  get_double("mass", &m);
  get_double("U", &U);
  get_double("mu", &mu);

  get_long("Random seed", &seed);

  get_char(" Configuration filename ", configuration_filename);


  get_double("Wang Landau step size", &llr_alpha);
  get_int("Wang Landau steps with dampened decay", &llr_constant_steps);
  get_char(" Initial values file ", init_parameter_filename);

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  //printf(" Git commit ID GIT_COMMIT_ID  \n");
  printf(" 2D Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Size of fluctuation matrix %d\n", max_changes );
  printf(" Random seed %ld\n", seed );
  printf(" Configuration file %s\n", configuration_filename );
  printf(" Wang Landau step size %g\n", llr_alpha );
  printf(" Wang Landau %d first steps with dampened decay\n", llr_constant_steps );
  printf(" Wang Landau weight file %s\n", init_parameter_filename );

  
  // Set up lattice variables and fields
  setup_lattice(seed);
  read_configuration(configuration_filename);
  
  int sum_sign=0;
  WangLaundau_setup( n_average*n_measure );

  struct timeval start, end;
  double updatetime=0, measuretime = 0;
  gettimeofday(&start,NULL);

  /* and the update/measure loop */
  for (i=1; i<n_loops+1; i++) {

    /* Update */
    update(n_measure);


    /* Time and report */
    gettimeofday(&end,NULL);
    updatetime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    gettimeofday(&start,NULL);

    int sector = negative_loops();
    int sign = 1-(sector%2)*2;
    sum_sign += sign;

    // Update the free energy in the sector
    WangLaundau_update(sector);

    gettimeofday(&end,NULL);
    measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    if(i%n_average==0){
      printf("\n%d, %d updates in %.3g seconds\n", i*n_measure, n_average*n_measure, 1e-6*updatetime);
      printf("%d, %d measurements in %.3g seconds\n", i*n_measure, n_average, 1e-6*measuretime);
      printf("%d, acceptance %.3g, sector change rate %.3g \n", i*n_measure, (double)llr_accepted/n_average, (double)sector_changes/n_average);
      llr_accepted = 0; sector_changes = 0;

      updatetime = 0; measuretime = 0;

      printf("SIGN %g\n", (double)sum_sign/n_average);

      WangLaundau_write_energy();
      write_configuration(configuration_filename);

      sum_sign = 0;
    }
      
    gettimeofday(&start,NULL);
  }

  printf(" ** simulation done\n");

  return(0);
}

