
#include "mex.h"
#include "math.h"

#define N_DATA 41
#define N_INPUT 0
#define N_PARS 28
#define N_STATES 25
#define N_OBSPARS 1
#define N_ICPARS 0
#define N_OBS 42
#define N_TIME 41
#define SENSDIM 1325

void objectiveFn(double *obj, double *s, double *p, double *d, double *u );

