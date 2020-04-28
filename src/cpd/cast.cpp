#include "cast.h"

uintptr_t
warthog::cpd::wt_to_label(double& b)
{
    uintptr_t ret;
    memcpy(&ret, &b, sizeof(double));
    return ret;
}

double
warthog::cpd::label_to_wt(uintptr_t &b)
{
  double ret;
  memcpy(&ret, &b, sizeof(uintptr_t));
  return ret;
}
