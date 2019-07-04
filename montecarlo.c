#define MAIN
#include "worldline.h"


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
  double sum_llr_a = 0;

  struct timeval start, end;
  double updatetime=0, measuretime = 0;
  gettimeofday(&start,NULL);

  for (i=1; i<n_loops+1; i++) {

    /* Update */
    update_config(n_measure);

    /* Time and report */
    gettimeofday(&end,NULL);
    updatetime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    gettimeofday(&start,NULL);

    int sector = negative_loops();
    int sign = 1-(sector%2)*2;
    sum_sign += sign;

    // Just count hits to each sector
    sectors[sector] += 1;
    if( m == 0 )
      sum_susc_wb += sign*measure_susceptibility_with_background();

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
      if( m == 0 )
        printf("SUSCEPTIBILITY %g \n", (double)sum_susc_wb/n_average);
      printf("SIGN %g\n", (double)sum_sign/n_average);

      for(int s=0; s<MAX_SECTOR; s++){
        printf("SECTOR %d %g \n", s, (double)sectors[s]/n_average);
        sectors[s] = 0;
      }


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
