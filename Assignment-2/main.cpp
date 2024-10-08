#include <iostream>

#include <Eigen/Dense>

#include "math/math.h"

//Task 2: Wrenches
//T.2 a) Change reference frame of a wrench as measured by a force/torque sensor.
void task2a()
{
    Eigen::Vector3d f_w{-30, 0, 0};
    Eigen::Vector3d m_s{0, 0, 2};
    Eigen::Vector3d e_ws_yzx{60, -60, 0};

    Eigen::Matrix3d rot_ws = math::rotation_matrix_from_euler_yzx(e_ws_yzx);

    Eigen::Vector3d m_w = rot_ws * m_s;
    Eigen::Vector3d f_s = rot_ws.transpose() * f_w;

    std::cout << "f_w:";
    std::cout << f_w.transpose() << std::endl;
    std::cout << "m_w:";
    std::cout << m_w.transpose() << std::endl << std::endl;
    std::cout << "f_s:";
    std::cout << f_s.transpose() << std::endl;
    std::cout << "m_s:";
    std::cout << m_s.transpose() << std::endl << std::endl;
}
//T.2 b) Calculate the sum of wrenches expressed in different reference frames.
void task2b()
{
    double m_a = 0.1; //kg
    double g = 10; //(m/s^2)
    double m_h = 0.5; //kg
    double L_1 = 0.1; //meters
    double L_2 = 0.15; //meters

    Eigen::Vector3d w_h{0,0,0};
    Eigen::Vector3d f_h{0,-g*m_h,0};

    Eigen::Vector3d w_a{0,0,0};
    Eigen::Vector3d f_a{0,0,g*m_a};

    Eigen::VectorXd F_h = math::twist(w_h,f_h);
    Eigen::VectorXd F_a = math::twist(w_a,f_a);

    Eigen::Matrix3d rot_hf = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d rot_af;
    rot_af << 1, 0, 0,
    0, 0, 1,
    0, -1, 0;

    Eigen::Vector3d p_hf {-L_1, 0, 0};
    Eigen::Vector3d p_af {-(L_1+L_2), 0, 0};

    Eigen::Matrix4d T_hf = math::transformation_matrix(rot_hf,p_hf);
    Eigen::Matrix4d T_af = math::transformation_matrix(rot_af,p_af);

    Eigen::MatrixXd Ad_hf = math::adjoint_matrix(T_hf);
    Eigen::MatrixXd Ad_af = math::adjoint_matrix(T_af);

    Eigen::VectorXd F_f = Ad_hf.transpose()*F_h + Ad_af.transpose()*F_a;

    std::cout << "F_f:";
    std::cout << F_f.transpose() << std::endl << std::endl;

}

int main() {
    //Task 2 a)
    std::cout << std::endl << "TASK 2A" << std::endl;
    task2a();

    //Task 2 b)
    std::cout << std::endl << "TASK 2B" << std::endl;
    task2b();

    //Task 4: Forward kinematics: 3R planar open chain
    std::cout << std::endl << "TASK 4" << std::endl;
    std::vector<double> L = {10,10,10};

    std::vector<double> j4_1 = {0,0,0};
    Eigen::Matrix4d T4_1 = math::planar_3r_fk_transform(L, j4_1);
    std::wcout << "Configuration 1: ";
    for (double i: j4_1){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 1 Using D-H", T4_1);
    Eigen::Matrix4d T4s_1 = math::planar_3r_fk_screw(L, j4_1);
    math::print_pose("Configuration 1 Using PoE", T4s_1);

    std::vector<double> j4_2 = {90,0,0};
    Eigen::Matrix4d T4_2 = math::planar_3r_fk_transform(L, j4_2);
    std::wcout << "Configuration 2: ";
    for (double i: j4_2){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 2 Using D-H", T4_2);
    Eigen::Matrix4d T4s_2 = math::planar_3r_fk_screw(L, j4_2);
    math::print_pose("Configuration 2 Using PoE", T4s_2);

    std::vector<double> j4_3 = {0,90,0};
    Eigen::Matrix4d T4_3 = math::planar_3r_fk_transform(L, j4_3);
    std::wcout << "Configuration 3: ";
    for (double i: j4_3){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 3 Using D-H", T4_3);
    Eigen::Matrix4d T4s_3 = math::planar_3r_fk_screw(L, j4_3);
    math::print_pose("Configuration 3 Using PoE", T4s_3);

    std::vector<double> j4_4 = {0,0,90};
    Eigen::Matrix4d T4_4 = math::planar_3r_fk_transform(L, j4_4);
    std::wcout << "Configuration 4: ";
    for (double i: j4_4){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 4 Using D-H", T4_4);
    Eigen::Matrix4d T4s_4 = math::planar_3r_fk_screw(L, j4_4);
    math::print_pose("Configuration 4 Using PoE", T4s_4);

    std::vector<double> j4_5 = {10,-15,2.75};
    Eigen::Matrix4d T4_5 = math::planar_3r_fk_transform(L, j4_5);
    std::wcout << "Configuration 5: ";
    for (double i: j4_5){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 5 Using D-H", T4_5);
    Eigen::Matrix4d T4s_5 = math::planar_3r_fk_screw(L, j4_5);
    math::print_pose("Configuration 5 Using PoE", T4s_5);


    //Task 5: Forward kinematics: UR3e 6R open chain
    std::cout << std::endl << "TASK 5" << std::endl;
    std::vector<double> j_example = {0,-90,0,0,90,0};
    Eigen::Matrix4d T_example = math::ur3e_fk_transform(j_example);
    std::wcout << "Configuration Example: ";
    for (double i: j_example){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration Example Using D-H", T_example);
    Eigen::Matrix4d Ts_example = math::ur3e_fk_screw(j_example);
    math::print_pose("Configuration Example Using PoE", Ts_example);

    std::vector<double> j_1 = {0,0,0,-90,0,0};
    Eigen::Matrix4d T_1 = math::ur3e_fk_transform(j_1);
    std::wcout << "Configuration 1: ";
    for (double i: j_1){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 1 Using D-H", T_1);
    Eigen::Matrix4d Ts_1 = math::ur3e_fk_screw(j_1);
    math::print_pose("Configuration 1 Using PoE", Ts_1);

    std::vector<double> j_2 = {0,-180,0,0,0,0};
    Eigen::Matrix4d T_2 = math::ur3e_fk_transform(j_2);
    std::wcout << "Configuration 2: ";
    for (double i: j_2){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 2 Using D-H", T_2);
    Eigen::Matrix4d Ts_2 = math::ur3e_fk_screw(j_2);
    math::print_pose("Configuration 2 Using PoE", Ts_2);

    std::vector<double> j_3 = {0,-90,0,0,0,0};
    Eigen::Matrix4d T_3 = math::ur3e_fk_transform(j_3);
    std::wcout << "Configuration 3: ";
    for (double i: j_3){std::cout << i;std::cout << ",";}std::cout << std::endl;
    math::print_pose("Configuration 3 Using D-H", T_3);
    Eigen::Matrix4d Ts_3 = math::ur3e_fk_screw(j_3);
    math::print_pose("Configuration 3 Using PoE", Ts_3);

    return 0;
}