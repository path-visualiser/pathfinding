#ifndef __GEOGRAPHY_H_
#define __GEOGRAPHY_H_

#include <bits/stdint-intn.h>

namespace warthog
{

namespace geo
{

// Get the bearing from pos1 to pos2
double
get_bearing(double lng1, double lat1, double lng2, double lat2);

// Get the bearing of a location wrt North
double
true_bearing(double lng, double lat);

// Get angle between (AB) and (AC) given their point coordinates
double
get_angle(double xa, double ya, double xb, double yb, double xc, double yc);

}
}

#endif // __GEOGRAPHY_H_
