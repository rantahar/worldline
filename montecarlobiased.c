

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
    printf(" Reading initial free energy\n" );
    for( int s=0; s<MAX_SECTOR; s++){
      long tmp1; int tmp2;
      fscanf(config_file, "%lf %ld %d\n", &Sampling_F[s], &tmp1, &tmp2);
    }
    fclose(config_file);
  } else {
    errormessage("Could not read weight file!\n");
  }

  current_sector = negative_loops();
}



int transition_accept(){
  int accept = 1;
  int sector;
  double weight;
  sector = negative_loops();
  if( sector != current_sector ){
    double minimum =-40;
    double e1 = Sampling_F[current_sector];
    double e2 = Sampling_F[sector];
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

  int i,n_loops,n_measure,n_average,llr_update_every;
  long seed;

  /* Read in the input */
  get_int("Sites in the t direction", &NT);
  get_int("Sites in the x direction", &NX);
  get_int("Number of updates", &n_loops);
  get_int("Updates / measurement", &n_measure);
  get_int("Average over", &n_average);
  get_double("mass", &m);
  get_double("U", &U);
  get_double("mu", &mu);
  get_long("Random seed", &seed);
  get_char(" Configuration filename ", configuration_filename);

  get_char("Initial values file ", parameter_filename);


  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  //printf(" Git commit ID GIT_COMMIT_ID  \n");
  printf(" 2D Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );
  printf(" Configuration file %s\n", configuration_filename );
  printf(" Weight file %s\n", parameter_filename );

  
  // Set up lattice variables and fields
  setup_lattice(seed);
  read_configuration(configuration_filename);
  
  int sum_sign=0;
  setup_free_energy( n_average*n_measure );
  double sectors[MAX_SECTOR];
  for(int i=0; i<MAX_SECTOR; i++){
    sectors[i] = 0;
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

