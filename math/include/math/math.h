#ifndef MATH_H
#define MATH_H

#include <Eigen/Dense>

namespace math
{
    //Assignment 1:
    //Assignment 1. Task 1: Skew symmetric matrix
    constexpr double deg_to_rad = 0.0174532925;
    constexpr double rad_to_deg = 57.2957795;
    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);
    Eigen::Vector3d vector_from_skew_symmetric(Eigen::Matrix3d skew);
    bool floatEquals(double a, double b);
    //Assignment 1. Task 2: Rotation Matrices
    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z);
    Eigen::Matrix3d rotate_x(double degrees);
    Eigen::Matrix3d rotate_y(double degrees);
    Eigen::Matrix3d rotate_z(double degrees);
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees);
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_zxy(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_yxz(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_xyz(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_xzy(const Eigen::Vector3d &e);
    //Assignment 1. Task 3: Transformation matrices
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);

    //Assignment 2:
    //Assignment 2. Task 1: Functions and algorithms
    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::Vector3d vector_cross_product(const Eigen::Vector3d &a, const Eigen::Vector3d &b);
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
    double cot(double deg);
    //Assignment 2. Task 3: Matrix exponentials and logarithms
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double t_deg);
    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double t_deg);
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t);
    //Assignment 2. Task 4: Forward kinematics: 3R planar open chain
    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);
    Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &link_L,const std::vector<double> &joint_positions);
    Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &link_L,const std::vector<double> &joint_positions);
    //Assignment 2. Task 5: Forward kinematics: UR3e 6R open chain
    Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions);
    Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions);

    //Assignment 3:
    //Assignment 3. Task 1: Functions and algorithms
    Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v);
    bool is_average_below_eps(const std::vector<double> &values, double eps, uint8_t n_values);
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_space_chain();
    Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions);
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain();
    Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions);
    //Assignment 3. Task 2: Numerical optimization
    std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f, double x_0, double dx_0, double eps);
    std::pair<uint32_t, double> gradient_descent_root_find(const std::function<double(double)> &f, double x_0, double gamma, double dx_0, double eps);
    //Assignment 3. Task 3: Velocity kinematics: UR3e 6R open chain
    Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions);
    Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions);
    //Assignment 3. Task 4: Inverse kinematics: UR3e 6R open chain
    std::pair<int, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &current_joint_positions, double gamma, double v_e, double w_e);

}

#endif


