#ifndef MATH_H
#define MATH_H

#include <Eigen/Dense>

namespace math
{
    double deg_to_rad(double degrees);
    double rad_to_deg(double radians);
    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);
    Eigen::Vector3d vector_from_skew_symmetric(Eigen::Matrix3d skew);
    bool floatEquals(double a, double b);
    Eigen::Matrix3d rotate_x(double degrees);
    Eigen::Matrix3d rotate_y(double degrees);
    Eigen::Matrix3d rotate_z(double degrees);
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_zxy(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_yxz(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_xyz(const Eigen::Vector3d &e);
    Eigen::Matrix3d rotation_matrix_from_euler_xzy(const Eigen::Vector3d &e);
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);
    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::Vector3d vector_cross_product(const Eigen::Vector3d &a, const Eigen::Vector3d &b);
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
    double cot(double deg);
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double t_deg);
    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double t_deg);
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t);
    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);
    Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &link_L,const std::vector<double> &joint_positions);
    Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &link_L,const std::vector<double> &joint_positions);
    Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions);
    Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions);


}

#endif


