#include <iostream>
#include <utility>
#include <Eigen/Dense>
#include <math/math.h>
#include <functional>

//Assignment 3
//Task 1: Functions and algorithms
void ur3e_test_fk()
{
    std::cout << "Forward kinematics tests" << std::endl;
    math::print_pose("Joint Configuration 0 SPACE",math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
    math::print_pose("Joint Configuration 0 BODY",math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
    std::cout << std::endl;

    math::print_pose("Joint Configuration 1 SPACE",math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0}) * math::deg_to_rad));
    math::print_pose("Joint Configuration 1 BODY",math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0}) * math::deg_to_rad));
    std::cout << std::endl;

    math::print_pose("Joint Configuration 2 SPACE",math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
    math::print_pose("Joint Configuration 2 BODY",math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
    std::cout << std::endl;

    math::print_pose("Joint Configuration 3 SPACE",math::ur3e_space_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
    math::print_pose("Joint Configuration 3 BODY",math::ur3e_body_fk(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0}) * math::deg_to_rad));
}

//Task 2: Numerical optimization
void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0) {
    double eps = 10e-7;
    double dx = 0.5;
    auto [iterations, x_hat] = math::newton_raphson_root_find(f, x0, dx, eps);
    std::cout << "NR root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" <<
    f(x_hat) << std::endl;
}
void test_gradient_descent_root_find(const std::function<double(double)> &f, double x0) {
    double eps = 10e-7;
    double dx = 0.5;
    double gamma = 0.01;
    auto [iterations, x_hat] = math::gradient_descent_root_find(f, x0, gamma, dx, eps);
    std::cout << "GD root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" <<
    f(x_hat) << std::endl;
}
void test_root_find()
{
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x){return (x - 3.0) * (x - 3.0) - 1.0;};
    test_newton_raphson_root_find(f1, -20);
    test_gradient_descent_root_find(f1, -20);
}

//Task 3: Velocity kinematics: UR3e 6R open chain
void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions){
    Eigen::Matrix4d tsb = math::ur3e_body_fk(joint_positions);
    auto [m, space_screws] = math::ur3e_space_chain();
    Eigen::MatrixXd jb = math::ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = math::ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = math::adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = math::adjoint_matrix(tsb.inverse());
    std::cout << "Jb: " << std::endl << jb << std::endl << "Ad_tbs*Js:" << std::endl << ad_tbs * js <<
    std::endl << std::endl;
    std::cout << "Js: " << std::endl << js << std::endl << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb <<
    std::endl << std::endl;
    std::cout << "d Jb: " << std::endl << jb - ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js - ad_tsb * jb << std::endl << std::endl;
}
void ur3e_test_jacobian(){
    std::cout << "Jacobian matrix tests" << std::endl;
    ur3e_test_jacobian(math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad);
    ur3e_test_jacobian(math::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) * math::deg_to_rad);
}

//Task 4: Inverse kinematics: UR3e 6R open chain
void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0)
{
    std::cout << "Test from pose" << std::endl;
    Eigen::Matrix4d t_sd = math::transformation_matrix(math::rotation_matrix_from_euler_zyx(zyx), pos);
    //std::cout << "t_sd: " << t_sd << std::endl;
    double gamma = 1e-2; double v_e = 4e-3; double w_e = 4e-3;
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0, gamma, v_e, w_e);
    Eigen::Matrix4d t_ik = math::ur3e_space_fk(j_ik);
    math::print_pose(" IK pose",t_ik);
    math::print_pose("Desired pose",t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg << std::endl << std::endl;
}
void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0)
{
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = math::ur3e_space_fk(joint_positions);
    double gamma = 1e0; double v_e = 4e-3; double w_e = 4e-3;
    auto [iterations, j_ik] = math::ur3e_ik_body(t_sd, j0, gamma, v_e, w_e);
    Eigen::Matrix4d t_ik = math::ur3e_space_fk(j_ik);
    math::print_pose(" IK pose",t_ik);
    math::print_pose("Desired pose",t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg << std::endl << std::endl;
}
void ur3e_ik_test()
{
    Eigen::VectorXd j_t0 = math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad;
    Eigen::VectorXd j_t1 = math::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -89.0, 0.0, 0.0, 0.0}) * math::deg_to_rad;
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505},Eigen::Vector3d{0.0, 90.0, -90.0}, j_t0);
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505},Eigen::Vector3d{0.0, 90.0, -90.0}, j_t1);

    Eigen::VectorXd j_t2 = math::std_vector_to_eigen(std::vector<double>{50.0, -30.0, 20.0, 0.0, -30.0, 50.0}) * math::deg_to_rad;
    Eigen::VectorXd j_d1 = math::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0,-50.0}) * math::deg_to_rad;
    ur3e_ik_test_configuration(j_d1, j_t0);
    ur3e_ik_test_configuration(j_d1, j_t2);
}

int main() {
    //Task 1: Functions and algorithms
    std::cout << std::endl << "TASK 1" << std::endl;
    ur3e_test_fk();

    //Task 2: Numerical optimization
    std::cout << std::endl << "TASK 2" << std::endl;
    test_root_find();

    //Task 3: Velocity kinematics: UR3e 6R open chain
    std::cout << std::endl << "TASK 3" << std::endl;
    ur3e_test_jacobian();

    //Task 4: Inverse kinematics: UR3e 6R open chain
    std::cout << std::endl << "TASK 4" << std::endl;
    ur3e_ik_test();

    return 0;
}