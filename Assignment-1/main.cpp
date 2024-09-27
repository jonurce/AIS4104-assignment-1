#include <iostream>

#include <Eigen/Dense>

#include "math/math.h"

//Task 1: Skew symmetric matrix
//1.a) Create a function that calculates the skew-symmetric matrix representation of a vector v ∈ R3
Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v)
{
    Eigen::Matrix3d matrix;
    matrix << 0, -v(2), v(1),
    v(2), 0, -v(0),
    -v(1), v(0), 0;
    return matrix;
}

//1.b) Verify the function for constructing the skew-symmetric matrix
void skew_symmetric_test()
{
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}


//Task 2: Rotation Matrices
//2.a)Create a rotation matrix from reference frame axes
Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
{
    Eigen::Matrix3d matrix;
    matrix << x(0), y(0), z(0),
        x(1), y(1), z(1),
        x(2), y(2), z(2)
        ;
    return matrix;
}

//2.b)Create a rotation matrix from rotating θ degrees about the principal axis x
Eigen::Matrix3d rotate_x(double degrees)
{
    double radians = math::deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix << 1, 0, 0,
    0, std::cos(radians), -std::sin(radians),
    0, std::sin(radians), std::cos(radians);
    return matrix;
}

//2.c)Create a rotation matrix from rotating θ degrees about the principal axis y
Eigen::Matrix3d rotate_y(double degrees)
{
    double radians = math::deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix << std::cos(radians), 0, std::sin(radians),
    0, 1, 0,
    -std::sin(radians), 0, std::cos(radians);
    return matrix;
}

//2.d) Create a rotation matrix from rotating θ degrees about the principal axis z
Eigen::Matrix3d rotate_z(double degrees)
{
    double radians = math::deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    matrix << std::cos(radians), -std::sin(radians), 0,
    std::sin(radians), std::cos(radians), 0,
    0, 0, 1;
    return matrix;
}

//2.e)Create a rotation matrix from rotating θ degrees about an arbitrary axis ω
Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees)
{
    double radians = math::deg_to_rad(degrees);
    Eigen::Matrix3d matrix;
    double c = std::cos(radians);
    double s = std::sin(radians);
    double w1 = axis(0);
    double w2 = axis(1);
    double w3 = axis(2);
    matrix << c+(w1*w1)*(1-c), w1*w2*(1-c)-w3*s, w1*w3*(1-c)+w2*s,
    w1*w2*(1-c)+w3*s, c+(w2*w2)*(1-c), w2*w3*(1-c)-w1*s,
    w1*w3*(1-c)-w2*s, w2*w3*(1-c)+w1*s, c+(w3*w3)*(1-c);
    return matrix;
}

//2.f)Create a rotation matrix from Euler angles (ZYX)
Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = rotate_z(e(0));
    Eigen::Matrix3d rot_ybeta = rotate_y(e(1));
    Eigen::Matrix3d rot_xgamma = rotate_x(e(2));
    matrix=rot_zalpha*rot_ybeta*rot_xgamma;
    return matrix;
    }

//2.g)Verify the functions for constructing rotation matrices
void rotation_matrix_test()
{
    Eigen::Matrix3d rot =
    rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa =
    rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa =
    rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
    Eigen::Vector3d{-0.5, -0.5, 0.707107},
    Eigen::Vector3d{0.707107, -0.707107, 0.0});
    std::cout << "Rotation matrix from Euler: " << std::endl;
    std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl;
    std::cout << "Rotation matrix from frame axes: " << std::endl;
    std::cout << rot_fa << std::endl << std::endl;
}


//Task 3: Transformation matrices
//3.a)Create a transformation matrix from a rotation matrix and translation vector
Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d matrix;
    matrix << r(0,0), r(0,1), r(0,2), p(0),
    r(1,0), r(1,1), r(1,2), p(1),
    r(2,0), r(2,1), r(2,2), p(2),
    0,0,0,1;
return matrix;
}

//3.b)Verify the function for creating a transformation matrix
void transformation_matrix_test()
{
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl;
}

//3.c)Transform vector from body-frame to fixed-frame coordinates
void transform_vector()
{
    Eigen::Vector3d e_rad{math::deg_to_rad(60.0), math::deg_to_rad(45.0), math::deg_to_rad(0.0)};
    Eigen::Vector3d p_wa{0.0, 0.0, 1.0};
    Eigen::Vector3d v_a{2.5, 3.0, -10};

    //As we are transforming a vector, we do not need to take into account the translation of the frames
    Eigen::Matrix3d rot_wa = rotation_matrix_from_euler_zyx(e_rad);
    Eigen::Vector3d v_w = rot_wa*v_a;

    /*If the vector was representing a point in space, we would need to take into account both
    the translation a rotation. This would be the code:
    Eigen::Matrix4d T_wa = transformation_matrix(rot_wa, p_wa);
    Eigen::Vector4d v4d_a{2.5, 3.0, -10, 1.0};
    Eigen::Vector4d v4d_w = T_wa*v4d_a;
    Eigen::Vector3d v3d_w{v4d_w(0), v4d_w(1), v4d_w(2)};
    */

    std::cout << "Vector v expressed in w coordinates: " << std::endl;
    std::cout << v_w << std::endl;
    std::cout << std::endl << std::endl;
}

int main()
{
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    transform_vector();
    return 0;
}
