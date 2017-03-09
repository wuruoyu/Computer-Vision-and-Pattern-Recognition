/*
  C version. Based on Java code by Tom Gibara
  Rewritten by RuoyuWu
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ffabs(x) ( (x) >= 0 ? (x) : -(x) ) 
#define GAUSSIAN_CUT_OFF 0.005f
#define MAGNITUDE_SCALE 100.0f
#define MAGNITUDE_LIMIT 1000.0f
#define MAGNITUDE_MAX ((int) (MAGNITUDE_SCALE * MAGNITUDE_LIMIT))
