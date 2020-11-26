// Module to handle geographic coordinates (lat/long)
#ifndef __GEOGRAPHY_H_
#define __GEOGRAPHY_H_

#include <bits/stdint-uintn.h>

namespace warthog
{

namespace geo
{

static const double DIMACS_RATIO = 1e6;

// Get the bearing from pos1 to pos2
double
get_bearing(double lat1, double lng1, double lat2, double lng2);

// Get the bearing of a location wrt North
double
true_bearing(double lat, double lng);

// Get angle between (ABC) given the lat/long of the three points
double
get_angle(double lat_a, double lng_a, double lat_b, double lng_b, double lat_c,
          double lng_c);

// Check whether, for three vectors rooted at the same point s, (sb) is between
// (sa) and (sc).
bool
between(double lat_a, double lng_a, double lat_b, double lng_b, double lat_c,
        double lng_c);

// Check whether vector (sb) is between (sa) and (sc)
bool
between(double lat_s, double lng_s, double lat_a, double lng_a, double lat_b,
        double lng_b, double lat_c, double lng_c);

}
}

#endif // __GEOGRAPHY_H_
