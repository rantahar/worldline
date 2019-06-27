#ifndef SAD_H
#define SAD_H
#define SAD

#include "worldline.h"

/* LLR parameters */
EXTERN int llr_target;
EXTERN double llr_gaussian_weight; // Used in thermalisation even with wall
EXTERN double llr_a;             // The measurable a in the llr method
EXTERN int llr_constant_steps; // Number of (roughly) constant steps at start
EXTERN double llr_alpha;       // Step size


EXTERN double WangLaundau_F[MAX_SECTOR];
EXTERN long WangLaundau_iteration[MAX_SECTOR];
EXTERN int WL_measure_sector[MAX_SECTOR];
EXTERN char init_parameter_filename[100];
EXTERN int current_sector;
EXTERN int llr_accepted;
EXTERN int sector_changes;
EXTERN long WL_nstep;
EXTERN int last_range_update;



#endif