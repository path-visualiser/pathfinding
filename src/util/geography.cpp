#include "geography.h"
#include <cmath>

double
deg_to_rad(double deg)
{
    return deg * M_PI / 180;
}

double
rad_to_deg(double rad)
{
    return fmod(rad * 180 / M_PI + 360, 360);
}

// We calculate bearing with the following formula:
// θ = atan2( sin Δλ ⋅ cos φ2 , cos φ1 ⋅ sin φ2 − sin φ1 ⋅ cos φ2 ⋅ cos Δλ )
// where:
//   φ1,λ1 is the start point,
//   φ2,λ2 the end point
//   Δλ is the difference in longitude
double
warthog::geo::get_bearing(double lat1, double lng1, double lat2, double lng2)
{
   double l1 = deg_to_rad(lng1);
   double l2 = deg_to_rad(lng2);
   double p1 = deg_to_rad(lat1);
   double p2 = deg_to_rad(lat2);
   double y = sin(l2 - l1) * cos(p2);
   double x = cos(p1) * sin(p2) - sin(p1) * cos(p2) * cos(l2 - l1);

   return rad_to_deg(std::atan2(y, x));
}

// We use the magnetic North (86°26′52.8″N 175°20′45.06″E) as to avoid angles
// being 180° for all points.
double
warthog::geo::true_bearing(double lat, double lng)
{
    return warthog::geo::get_bearing(
        deg_to_rad(86.448), deg_to_rad(175.5968),
        deg_to_rad(lng), deg_to_rad(lat));
}

// The angle (ABC) is defined as the bearing from C to A when B is the origin.
double
warthog::geo::get_angle(double lat_a, double lng_a, double lat_b, double lng_b,
                        double lat_c, double lng_c)
{
    return warthog::geo::get_bearing(
        lat_a - lat_b, lng_a - lng_b, lat_c - lat_b, lng_c - lng_b);
}

// This works by checking that the longitude of A and C differ in sign wrt B.
bool
warthog::geo::between(double lat_a, double lng_a, double lat_b, double lng_b,
                      double lat_c, double lng_c)
{
    return (lng_a - lng_b) * (lng_c - lng_b) < 0;
}

// Variant of 'between' where we pass the origin.
bool
warthog::geo::between(double lat_o, double lng_o, double lat_a, double lng_a,
                      double lat_b, double lng_b, double lat_c, double lng_c)
{
    double lng_sa = lng_o - lng_a;
    double lng_sb = lng_o - lng_b;
    double lng_sc = lng_o - lng_c;
    // TODO Offset latitudes too?
    return warthog::geo::between(lat_a, lng_sa, lat_b, lng_sb, lat_c, lng_sc);
}
