#define MAIN
#define MEASURE_SECTOR
#include "worldline.h"

int target_sector;
double gaussian_weight = 5; // Used in thermalisation even with wall
int current_sector = 0;


double sector_accept(){
  double logweight, weight;
  int distance, current_distance;
  int accept;
  int sector = negative_loops();
  distance = sector - target_sector;
  current_distance = current_sector - target_sector;
  logweight = -(distance*distance - current_distance*current_distance)*gaussian_weight;
  weight = exp(logweight);
  if( mersenne() < weight ){
    accept = 1;
    current_sector = sector;
  } else {
    accept = 0;
  }
  return accept;
}



/* Perform an update and accept/reject */
int update( int nsteps )
{
  int changes=0;

  for( int i=0; i<nsteps; i++ ){
    save_config();
    changes += update_config(1);
    if( ! sector_accept() ){
      restore_config();
    }
  }

  /* Only accept as real configuration if it matches the target sector */
  while(current_sector != target_sector){
    changes += update_config(1);
    current_sector = negative_loops();
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
  get_char("Configuration filename ", configuration_filename);

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  //printf(" Git commit ID GIT_COMMIT_ID  \n");
  printf(" 2D Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );
  printf(" Configuration file %s\n", configuration_filename );

  get_int("Target sector", &target_sector);
  printf(" Measuring expectation values in sector %d\n", target_sector );

  
  // Set up lattice variables and fields
  setup_lattice(seed);
  read_configuration(configuration_filename);
  
  /* and the update/measure loop */
  int sum_monomers = 0;
  int sum_links = 0;
  int sum_charge = 0;
  int sum_c2 = 0;
  int sum_q = 0;
  int sum_q2 = 0;
  double sum_susc = 0;
  double sum_susc_wb = 0;
  int sum_sign=0;
  int sectors[MAX_SECTOR];
  for(i=0; i<MAX_SECTOR; i++)
    sectors[i] = 0;

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

    int n_monomers = 0, n_links = 0;
    for (int t=0; t<NT; t++){
      for (int x=0; x<NX; x++) {
        if(field[t][x] == MONOMER){
          n_monomers +=1;
        }
        if(field[t][x] >= LINK_TUP){
          n_links +=1;
        }
      }
    }

    sum_monomers += sign*n_monomers;
    sum_links += sign*n_links;

    int c, q;
    measure_charge(&c, &q);
    sum_charge += sign*c;
    sum_c2 += sign*c*c;
    sum_q += sign*q;
    sum_q2 += sign*q*q; 

    gettimeofday(&end,NULL);
    measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    if(i%n_average==0){
      printf("\n%d, %d updates in %.3g seconds\n", i*n_measure, n_average*n_measure, 1e-6*updatetime);
      printf("%d, %d measurements in %.3g seconds\n", i*n_measure, n_average, 1e-6*measuretime);

      //measure_propagator( 1 ); //This includes an invertion and therefore takes time
        
      updatetime = 0; measuretime = 0;

      printf("MONOMERS %g \n", ((double)sum_monomers)/n_average);
      printf("LINKS %g \n", (double)sum_links/n_average);
      printf("CHARGE %g %g \n", (double)sum_charge/n_average, (double)sum_c2/n_average);
      printf("QCHI %g %g \n", (double)sum_q/n_average, (double)sum_q2/n_average);
      printf("SIGN %g\n", (double)sum_sign/n_average);

      write_configuration(configuration_filename);

      sum_monomers = 0; sum_links = 0; sum_charge = 0;
      sum_c2 = 0; sum_q = 0; sum_q2 = 0; sum_susc = 0;
      sum_susc_wb = 0; sum_sign = 0;
    }
      
    gettimeofday(&start,NULL);
  }

  printf(" ** simulation done\n");

  return(0);
}
