#define MAIN

#include "worldline.h"


int llr_target;
double llr_gaussian_weight = 5; // Used in thermalisation even with wall
double llr_a = 0;             // The measurable a in the llr method
int llr_constant_steps = 100; // Number of (roughly) constant steps at start
double llr_alpha = 0.1;       // Step size
int current_sector;
int llr_accepted;
int sector_changes;


// In LLR, modify the acceptance rate based on the
// number of negative loops
void LLR_update( double deltaS ){
  static int iter = 1;
  double step = llr_alpha*llr_constant_steps/(iter+llr_constant_steps);
  llr_a -= step*deltaS;
  iter ++;
}

// Get the modified weight of a sector
double LLR_weight( sector ){
  double distance, logweight, weight, a;
  distance = sector - llr_target-0.5;
  logweight = -(distance*distance-0.25)*llr_gaussian_weight;
  if( distance < 0 ){
    logweight += 0.5*llr_a;
  } else {
    logweight -= 0.5*llr_a;
  }
  weight = exp(logweight);
  return weight;
}

/* Check wether to accept a configuration */
int current_sector = 0;
int llr_accepted=0;
int sector_changes=0;
int llr_accept(){
  int accept = 1;
  int sector;
  double current_distance, previous_distance, weight;
  sector = negative_loops();
  if( sector != current_sector ){
    weight = LLR_weight(sector)/LLR_weight(current_sector);
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



/* Perform an update and accept/reject */
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

  int i,n_loops,n_measure,n_average,llr_update_every;

  setup_lattice();

  /* Read in the input */
  get_int("Number of updates", &n_loops);
  get_int("Updates / measurement", &n_measure);
  get_int("Updates between saves", &n_average);
  printf("\n %d updates per measurements\n", n_measure );

  read_thirring_parameters();

  get_double("LLR step size", &llr_alpha);
  get_int("LLR steps with dampened decay", &llr_constant_steps);
  get_int("Target LLR sector", &llr_target);
  get_double("LLR initial a", &llr_a);
  get_int("Updates / LLR update", &llr_update_every);

  printf("\n LLR target %d\n", llr_target );
  printf(" LLR updated every %ld updates\n", llr_update_every );
  printf(" LLR step size %g\n", llr_alpha );
  printf(" LLr %d first steps with dampened decay\n", llr_constant_steps );

  
  /* and the update/measure loop */
  int sum_sign=0;
  int sectors[MAX_SECTOR];
  for(i=0; i<MAX_SECTOR; i++)
    sectors[i] = 0;
  double sum_llr_a = 0;

  int sector=0;
  for (i=1;; i++) {
    // Wait for the target sector to be reached before
    // starting measurement runs
  
    update( 1 );
  
    sector = current_sector;
    if( sector == llr_target || sector == llr_target+1 ) {
      break;
    }
    if( i== n_loops ){
      printf( "Did not reach LLR target sector in %d updates\n", n_loops );
      exit(1);
    }
  }
  printf( "Reached LLR target sector in %d thermalisation updates\n", i );

  llr_accepted = 0;
  sector_changes = 0;

  struct timeval start, end;
  double updatetime=0, measuretime = 0;
  gettimeofday(&start,NULL);

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

    // Update the LLR transition propability
    if(i%llr_update_every==0){
      double llr_dS = (double)(sectors[llr_target]-sectors[llr_target+1])/(double)llr_update_every;
      LLR_update( llr_dS );
      sectors[llr_target] = 0;
      sectors[llr_target+1] = 0;
      sum_llr_a += llr_a;
    }

    gettimeofday(&end,NULL);
    measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    if(i%n_average==0){
      printf("\n%d, %d updates in %.3g seconds\n", i*n_measure, n_average*n_measure, 1e-6*updatetime);
      printf("%d, %d measurements in %.3g seconds\n", i*n_measure, n_average, 1e-6*measuretime);

      printf("%d, acceptance %.3g, sector change rate %.3g \n", i*n_measure, (double)llr_accepted/n_average, (double)sector_changes/n_average);
      llr_accepted = 0; sector_changes = 0;
        
      updatetime = 0; measuretime = 0;

      printf("SIGN %g\n", (double)sum_sign/n_average);

      double llr_a_ave = sum_llr_a/n_average*llr_update_every;
      printf("LLR a_%d = %g, exp(a) = %g\n", llr_target, llr_a_ave, exp(llr_a_ave));
      sum_llr_a = 0;

      write_configuration(configuration_filename);

      sum_sign = 0;
    }
      
    gettimeofday(&start,NULL);
  }

  printf(" ** simulation done\n");

  return(0);
}

