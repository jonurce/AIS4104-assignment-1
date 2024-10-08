#include <iostream>

#include <Eigen/Dense>

#include "math/math.h"

//Task 1: Skew symmetric matrix
//1.b) Verify the function for constructing the skew-symmetric matrix
void skew_symmetric_test()
{
    Eigen::Matrix3d skew_matrix = math::skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

//Task 2: Rotation Matrices
//2.g)Verify the functions for constructing rotation matrices
void rotation_matrix_test()
{
    Eigen::Matrix3d rot = math::rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa = math::rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa = math::rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
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
//3.b)Verify the function for creating a transformation matrix
void transformation_matrix_test()
{
    Eigen::Matrix3d r = math::rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << math::transformation_matrix(r, v) << std::endl;
}
//3.c)Transform vector from body-frame to fixed-frame coordinates
void transform_vector()
{
    Eigen::Vector3d e_rad{math::deg_to_rad*60.0, math::deg_to_rad*45.0, math::deg_to_rad*0.0};
    Eigen::Vector3d p_wa{0.0, 0.0, 1.0};
    Eigen::Vector3d v_a{2.5, 3.0, -10};

    //As we are transforming a vector, we do not need to take into account the translation of the frames
    Eigen::Matrix3d rot_wa = math::rotation_matrix_from_euler_zyx(e_rad);
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
    //Task 1: Skew symmetric matrix
    std::cout << std::endl << "TASK 1" << std::endl;
    skew_symmetric_test();

    //Task 2: Rotation Matrices
    std::cout << std::endl << "TASK 2" << std::endl;
    rotation_matrix_test();

    //Task 3: Transformation matrices
    std::cout << std::endl << "TASK 3" << std::endl;
    transformation_matrix_test();
    transform_vector();
    return 0;
}
