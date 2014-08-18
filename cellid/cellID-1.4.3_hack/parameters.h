//V1.2a The double variable were float in the previous version. 
//Using double for compatibility with GOptions library

extern double max_d_over_s_cut;
extern int max_pixels_per_cell;
extern int min_pixels_per_cell;
extern double background_reject_factor;

extern double max_split_d_over_minor;
extern double I_over_U_for_match;
#ifndef _parameters_
#define _parameters_
double max_d_over_s_cut;
int max_pixels_per_cell;
int min_pixels_per_cell;

double background_reject_factor;
double max_split_d_over_minor;
double I_over_U_for_match;
#endif
