#include <Eigen/Dense>
#include <iostream>
#include "math/math.h"

double math::deg_to_rad(double degrees)
{
    return degrees * 0.0174532925;
}

double math::rad_to_deg(double radians)
{
    return radians * 57.2957795;
}

double math::something(double a)
{
    return a * a;
}