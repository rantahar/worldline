

#define MAIN

#include "worldline.h"


double Sampling_F[MAX_SECTOR];
char parameter_filename[100];
int current_sector;


void setup_free_energy( ){
  FILE *config_file;
  config_file = fopen(parameter_filename, "rb");
  int initialized = 0;
  
  if(config_file) {
    double weight_sum = 0;
    printf(" Reading initial free energy\n" );
    for( int s=0; s<MAX_SECTOR; s++){
      char tmp1[100];
      fscanf(config_file, "%lf %s\n", &Sampling_F[s], tmp1);
      weight_sum += exp(Sampling_F[s]);
    }
    fclose(config_file);
    for( int s=0; s<MAX_SECTOR; s++){
      Sampling_F[s] -= log(weight_sum);
    }
  } else {
    errormessage("Could not read weight file!\n");
  }

  current_sector = negative_loops();
}


double subtracted_weight(int s){
#ifdef SUBTRACT
  int even = ((int)s/2)*2;
  int odd = even+1;
  double weight = 0.5*log(exp(Sampling_F[even]) + exp(Sampling_F[odd]));
  //double relative = exp(Sampling_F[even]) - exp(Sampling_F[odd]);
  //double logexp = log(fabs(relative));
  //weight -= logexp;
  //printf("Sector %d, weight %g, relative %g, logexp %g\n", s, weight, relative, logexp);
#else
  double weight = Sampling_F[s];
#endif
  return weight;
}

int transition_accept(){
  int accept = 1;
  int sector;
  double weight;
  sector = negative_loops();
  if( sector != current_sector ){
    double minimum =-40;
    double e1 = subtracted_weight(current_sector);
    double e2 = subtracted_weight(sector);
    if(e1 < minimum)
      e1 = minimum;
    if(e2 < minimum)
      e2 = minimum;
    weight = exp(e1-e2);
    //printf("From %d to %d, %g %g weight %g\n", current_sector, sector, Sampling_F[current_sector],Sampling_F[sector], weight);
    if( mersenne() < weight ){
      current_sector = sector;
      accept = 1;
    } else {
      accept = 0;
    }
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

  if( ! transition_accept() ){
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

  int i,n_loops,n_measure,n_average;

  setup_lattice();

  /* Read in the input */
  get_int("Number of updates", &n_loops);
  get_int("Updates / measurement", &n_measure);
  get_int("Updates between saves", &n_average);
  printf("\n %d updates per measurements\n", n_measure );

  read_thirring_parameters();

  get_char("Initial values file ", parameter_filename);

  printf("\n Weight file %s\n", parameter_filename );
  
  int sum_sign=0;
  setup_free_energy( n_average*n_measure );
  double sectors[MAX_SECTOR];
  for(int i=0; i<MAX_SECTOR; i++){
    sectors[i] = 0;
    printf("Free energy %d %g\n", i, Sampling_F[i]);
  }

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

    sectors[sector]++;

    gettimeofday(&end,NULL);
    measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    if(i%n_average==0){
      printf("\n%d, %d updates in %.3g seconds\n", i*n_measure, n_average*n_measure, 1e-6*updatetime);
      printf("%d, %d measurements in %.3g seconds\n", i*n_measure, n_average, 1e-6*measuretime);

      updatetime = 0; measuretime = 0;

      printf("SIGN %g\n", (double)sum_sign/n_average);
      for(int s=0; s<MAX_SECTOR; s++){
        printf("SECTOR %d %g\n", s, sectors[s]/(double)n_average);
        sectors[s] = 0;
      }

      write_configuration(configuration_filename);

      sum_sign = 0;
    }
      
    gettimeofday(&start,NULL);
  }

  printf(" ** simulation done\n");

  return(0);
}

