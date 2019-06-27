#ifndef LLR_H
#define LLR_H
/* LLR parameters */
#include "worldline.h"

EXTERN int llr_target;
EXTERN double llr_gaussian_weight = 5; // Used in thermalisation even with wall
EXTERN double llr_a = 0;             // The measurable a in the llr method
EXTERN int llr_constant_steps = 100; // Number of (roughly) constant steps at start
EXTERN double llr_alpha = 0.1;       // Step size

EXTERN int current_sector;
EXTERN int llr_accepted;
EXTERN int sector_changes;


#endif