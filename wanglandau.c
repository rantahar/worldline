

#define MAIN
//#define SAD

#include "worldline.h"


double llr_gaussian_weight; // Used in thermalisation even with wall
int constant_steps; // Number of (roughly) constant steps at start
double stepsize;       // Step size
double WangLaundau_F[MAX_SECTOR];
long WangLaundau_iteration[MAX_SECTOR];
int WL_measure_sector[MAX_SECTOR];
char parameter_filename[100];
int current_sector;
int WL_accepted;
int sector_changes;
long WL_nstep;
int last_range_update;
double Tmin = 10;
double tolerance = -40;
int max_Ns=1;



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
  config_file = fopen(parameter_filename, "rb");
  int initialized = 0;
  
  if(config_file) {
    printf(" Reading initial free energy\n" );
    for( int s=0; s<MAX_SECTOR; s++){
      fscanf(config_file, "%lf %ld\n", &WangLaundau_F[s], &WangLaundau_iteration[s]);
      if( WangLaundau_iteration[s] > 0 ){
        initialized = 1;
      }
    }
    fscanf(config_file, "%ld %d %d\n", &WL_nstep, &last_range_update, &max_Ns);
    //printf("%d %d %d\n", WL_nstep, last_range_update, max_Ns);
    fclose(config_file);
  }

  if( initialized == 0 ) {
    for( int s=0; s<MAX_SECTOR; s++){
      WangLaundau_iteration[s] = 0;
      WL_measure_sector[s] = 0;
    }
    //init_free_energy( max_init_steps );
  }

  current_sector = negative_loops();
}

void WangLaundau_write_energy(){
  FILE * config_file;

  config_file = fopen(parameter_filename,"wb");
  if (config_file){
    for( int s=0; s<MAX_SECTOR; s++){
      fprintf(config_file, "%g %ld\n", WangLaundau_F[s], WangLaundau_iteration[s]);
    }
    fprintf(config_file, "%d %d %d\n", WL_nstep, last_range_update, max_Ns);
    fclose(config_file);
  } else {
    printf("Could not write Wang Landau energy\n");
    exit(1);
  }
}

// Update the free energy in the Wang-Landau algorithm
void WangLaundau_update(sector){
  double maximum = -10000;
  double  min_F = 10000;
  double step;
  int Ns = 1;

  WangLaundau_iteration[sector]++;

  for( int s=0; s<MAX_SECTOR; s++){
    if( WangLaundau_F[s] < min_F ){
      min_F = WangLaundau_F[s];
    }
  }

  for( int s=0; s<MAX_SECTOR; s++){
    if( WangLaundau_F[s] > (min_F+0.1) ){
      Ns+=1;
      WL_measure_sector[s] = 1;
    } 
  }
  if( Ns > max_Ns ){
    printf("New sectors %d %d\n", WL_nstep, Ns);
    last_range_update = WL_nstep;
    max_Ns = Ns;
  }
  
#ifdef SAD
  
  double t0 = max_Ns/Tmin;
  double st2 = (double) WL_nstep * (double) WL_nstep;
  double l = last_range_update;
  double e = t0 + WL_nstep/l;
  double d = t0 + st2/(max_Ns*l);
  step = e/d;
  printf("SAD %g %d %d %g %g %g\n", step, last_range_update, WL_nstep, t0, e, d);
  WangLaundau_F[sector] += step;
  if( step <= 0 ){
    printf("Negative step in SAD %g %g %g\n",step, e, d);
    exit(1);
  }

#else
  step = stepsize*max_Ns*constant_steps/(WL_nstep+max_Ns*constant_steps);
  WangLaundau_F[sector] += step;
#endif

  for( int s=0; s<MAX_SECTOR; s++)
    if( WangLaundau_F[s] > maximum )
      maximum = WangLaundau_F[s];
  for( int s=0; s<MAX_SECTOR; s++){
    WangLaundau_F[s] -= maximum;
    if( WangLaundau_F[s] < tolerance ){
      WangLaundau_F[s] = tolerance;
    }
  }
  WL_nstep++;
}

double WangLaundau_weight(new_sector,old_sector){
  return exp(WangLaundau_F[old_sector]-WangLaundau_F[new_sector]);
}


int WL_accept(){
  int accept = 1;
  int sector;
  double weight;
  sector = negative_loops();
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
    WL_accepted += 1;
  }
  return accept;
}


/* A full update function. A single worm update followed by a number of random
   link and monomer updates */
int update( int nsteps )
{
  int changes=0;
  save_config();

  changes += update_config(nsteps);

  if( ! WL_accept() ){
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
  WL_accepted =0;
  sector_changes =0;
  WL_nstep =1;
  last_range_update =1;

  int i,n_loops,n_measure,n_average,llr_update_every;

  setup_lattice();

  /* Read in the input */
  get_int("Number of updates", &n_loops);
  get_int("Updates / measurement", &n_measure);
  get_int("Updates between saves", &n_average);
  printf("\n %d updates per measurements\n", n_measure );

  read_thirring_parameters();

  #ifdef SAD
  get_double("Wang Landau Tmin", &Tmin);
  printf("\n Wang Landau Tmin %g\n", Tmin );
  #else
  get_int("Wang Landau t0", &constant_steps);
  get_double("Wang Landau stepsize", &stepsize);
  printf("\n Wang Landau t0 %d\n", constant_steps );
  printf("\n Wang Landau step size %g\n", stepsize );
  #endif
  get_double("Wang Landau tolerance", &tolerance);
  get_char("Result filename ", parameter_filename);

  printf("\n Wang Landau tolerance %g\n", tolerance );
  printf("\n Wang Landau weight file %s\n", parameter_filename );
  
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

    int sector = current_sector;
    int sign = 1-(sector%2)*2;
    sum_sign += sign;

    // Update the free energy in the sector
    WangLaundau_update(sector);

    gettimeofday(&end,NULL);
    measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    if(i%n_average==0){
      printf("\n%d, %d updates in %.3g seconds\n", i*n_measure, n_average*n_measure, 1e-6*updatetime);
      printf("%d, %d measurements in %.3g seconds\n", i*n_measure, n_average, 1e-6*measuretime);
      printf("%d, acceptance %.3g, sector change rate %.3g \n", i*n_measure, (double)WL_accepted/n_average, (double)sector_changes/n_average);
      WL_accepted = 0; sector_changes = 0;

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

