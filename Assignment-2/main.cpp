#include <iostream>

#include <Eigen/Dense>

#include "math/math.h"

int main()
{
    std::cout << math:: deg_to_rad(180) << std::endl;
    std::cout << math:: something(2) ;
    return 0;
}

//Task 1: Functions and algorithms
//T.1 a) Calculate an Euler ZYX vector from a rotation matrix.

//Define a function to see if two numbers are equal, preventing floating point errors
bool floatEquals(double a, double b)
{
    return std::abs(a - b) < 1e-6;
}
 //Complete Task 1
Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    if(floatEquals(r(2,0),1.0))
    {
        b = -EIGEN_PI / 2.0;
        a = 0.0;
        c = -std::atan2(r(0,1), r(1,1));
    }
    else if(floatEquals(r(2,0),-1.0))
    {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0,1), r(1,1));
    }
    else
    {
        b = std::atan2(r(2,0), std::sqrt(r(0,0)*r(0,0)+r(1,0)*r(1,0)));
        a = std::atan2(r(1,0), r(0,0));
        c = std::atan2(r(2,1), r(2,2));
    }

    return Eigen::Vector3d(a, b, c);
}

//T.1 b) Create a twist vector from velocity vectors.

Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    return Eigen::VectorXd(w(0), w(1), w(2), v(0), v(1), v(2));
}

//T.1 c) Create a screw axis.

Eigen::Vector3d vector_cross_product

Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    Eigen::Vector3d v =  + s * h;
    return Eigen::VectorXd(s(0), s(1), s(2)), v(0), v(1), v(2);
}
