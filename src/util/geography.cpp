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

// Latitude of North pole is 90, longitude is 0
//
// In radians, this means that φ2 = π / 2 and λ2 = 0
//
// Further: cos φ2 = 0 and sin φ2 = 1 and Δλ = -λ1
//
// So: θ = atan2(0 , cos φ1)
//
// We use the magnetic North (86°26′52.8″N 175°20′45.06″E) as to avoid angles
// being 180° for all points.
double
warthog::geo::true_bearing(double lat, double lng)
{
    // double phi = deg_to_rad(lat);

    // return rad_to_deg(std::atan2(0, cos(phi)));
    return warthog::geo::get_bearing(
        deg_to_rad(86.44778), deg_to_rad(175.35483),
        deg_to_rad(lng), deg_to_rad(lat));
}

// Taken from: https://stackoverflow.com/questions/3486172/angle-between-3-points
double
warthog::geo::get_angle(double xa, double ya, double xb, double yb,
                        double xc, double yc)
{
    double xab = xb - xa;
    double yab = yb - ya;
    double xcb = xb - xc;
    double ycb = yb - yc;

    double dot = xab * xcb + yab * ycb; // dot product
    double cross = xab * ycb - yab * xcb; // cross product
    double alpha = atan2(cross, dot);

    return rad_to_deg(alpha);
}
