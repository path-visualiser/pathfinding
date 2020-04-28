#ifndef __CAST_H_
#define __CAST_H_

#include <cstdint>
#include <cstring>

namespace warthog
{

namespace cpd
{

// Hack of reinterpret_cast
uintptr_t
wt_to_label(double& b);

double
label_to_wt(uintptr_t &b);

} // namespace cpd

} // namespace warthog

#endif // __CAST_H_
