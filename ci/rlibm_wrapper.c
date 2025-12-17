#include <fenv.h>
#include "rlibm.h"

float log2f (float x) {
  int rnd = fegetround ();
  fesetround (FE_TONEAREST);
  double temp  = rlibm_log2f (x);
  fesetround (rnd);
  return (float) temp;
}
