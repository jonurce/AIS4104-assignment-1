#include <Eigen/Dense>
#include <iostream>
#include "math/math.h"

//Assignment 1:
//Assignment 1. Task 1: Skew symmetric matrix
Eigen::Matrix3d math::skew_symmetric(Eigen::Vector3d v)
{
    Eigen::Matrix3d matrix;
    matrix << 0, -v(2), v(1),
    v(2), 0, -v(0),
    -v(1), v(0), 0;
    return matrix;
}
Eigen::Vector3d math::vector_from_skew_symmetric(Eigen::Matrix3d skew)
{
    Eigen::Vector3d vector{skew(2,1), skew(0,2), skew(1,0)};
    return vector;
}
bool math::floatEquals(double a, double b)
{
    return std::abs(a - b) < 1e-3;
}
//Assignment 1. Task 2: Rotation Matrices
Eigen::Matrix3d math::rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
{
    Eigen::Matrix3d matrix;
    matrix << x(0), y(0), z(0),
        x(1), y(1), z(1),
        x(2), y(2), z(2)
        ;
    return matrix;
}
Eigen::Matrix3d math::rotate_x(double degrees)
{
    double radians = math::deg_to_rad*degrees;
    Eigen::Matrix3d matrix;
    matrix << 1, 0, 0,
    0, std::cos(radians), -std::sin(radians),
    0, std::sin(radians), std::cos(radians);
    return matrix;
}
Eigen::Matrix3d math::rotate_y(double degrees)
{
    double radians = math::deg_to_rad*degrees;
    Eigen::Matrix3d matrix;
    matrix << std::cos(radians), 0, std::sin(radians),
    0, 1, 0,
    -std::sin(radians), 0, std::cos(radians);
    return matrix;
}
Eigen::Matrix3d math::rotate_z(double degrees)
{
    double radians = math::deg_to_rad*degrees;
    Eigen::Matrix3d matrix;
    matrix << std::cos(radians), -std::sin(radians), 0,
    std::sin(radians), std::cos(radians), 0,
    0, 0, 1;
    return matrix;
}
Eigen::Matrix3d math::rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees)
{
    double radians = math::deg_to_rad*degrees;
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
Eigen::Matrix3d math::rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = math::rotate_z(e(0));
    Eigen::Matrix3d rot_ybeta = math::rotate_y(e(1));
    Eigen::Matrix3d rot_xgamma = math::rotate_x(e(2));
    matrix=rot_zalpha*rot_ybeta*rot_xgamma;
    return matrix;
}
Eigen::Matrix3d math::rotation_matrix_from_euler_zxy(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = math::rotate_z(e(0));
    Eigen::Matrix3d rot_ybeta = math::rotate_y(e(2));
    Eigen::Matrix3d rot_xgamma = math::rotate_x(e(1));
    matrix=rot_zalpha*rot_xgamma*rot_ybeta;
    return matrix;
}
Eigen::Matrix3d math::rotation_matrix_from_euler_yzx(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = math::rotate_z(e(1));
    Eigen::Matrix3d rot_ybeta = math::rotate_y(e(0));
    Eigen::Matrix3d rot_xgamma = math::rotate_x(e(2));
    matrix=rot_ybeta*rot_zalpha*rot_xgamma;
    return matrix;
}
Eigen::Matrix3d math::rotation_matrix_from_euler_yxz(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = math::rotate_z(e(2));
    Eigen::Matrix3d rot_ybeta = math::rotate_y(e(0));
    Eigen::Matrix3d rot_xgamma = math::rotate_x(e(1));
    matrix=rot_ybeta*rot_xgamma*rot_zalpha;
    return matrix;
}
Eigen::Matrix3d math::rotation_matrix_from_euler_xyz(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = math::rotate_z(e(2));
    Eigen::Matrix3d rot_ybeta = math::rotate_y(e(1));
    Eigen::Matrix3d rot_xgamma = math::rotate_x(e(0));
    matrix=rot_xgamma*rot_ybeta*rot_zalpha;
    return matrix;
}
Eigen::Matrix3d math::rotation_matrix_from_euler_xzy(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d rot_zalpha = math::rotate_z(e(1));
    Eigen::Matrix3d rot_ybeta = math::rotate_y(e(2));
    Eigen::Matrix3d rot_xgamma = math::rotate_x(e(0));
    matrix=rot_xgamma*rot_zalpha*rot_ybeta;
    return matrix;
}
//Assignment 1. Task 3: Transformation matrices
Eigen::Matrix4d math::transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d matrix;
    matrix << r(0,0), r(0,1), r(0,2), p(0),
    r(1,0), r(1,1), r(1,2), p(1),
    r(2,0), r(2,1), r(2,2), p(2),
    0,0,0,1;
    return matrix;
}

//Assignment 2:
//Assignment 2. Task 1: Functions and algorithms
Eigen::Vector3d math::euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
{
    double a_rad = 0.0; //Z axis Alpha
    double b_rad = 0.0; //Y axis Beta
    double c_rad = 0.0; //X axis Gamma

    if(math::floatEquals(r(2,0),1.0))
    {
        b_rad = -EIGEN_PI / 2.0;
        a_rad = 0.0;
        c_rad = -std::atan2(r(0,1), r(1,1));
    }
    else if(math::floatEquals(r(2,0),-1.0))
    {
        b_rad = EIGEN_PI / 2.0;
        a_rad = 0.0;
        c_rad = std::atan2(r(0,1), r(1,1));
    }
    else
    {
        b_rad = std::atan2(-r(2,0), std::sqrt(r(0,0)*r(0,0)+r(1,0)*r(1,0)));
        a_rad = std::atan2(r(1,0), r(0,0));
        c_rad = std::atan2(r(2,1), r(2,2));
    }

    double a_deg = math::rad_to_deg*a_rad; //Z axis Alpha
    double b_deg = math::rad_to_deg*b_rad; //Y axis Beta
    double c_deg = math::rad_to_deg*c_rad; //X axis Gamma

    return Eigen::Vector3d{a_deg, b_deg, c_deg};
}
Eigen::VectorXd math::twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::VectorXd tw(6);
    tw << w(0), w(1), w(2), v(0), v(1), v(2);
    return tw;
}
Eigen::Vector3d math::vector_cross_product(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
{
    return Eigen::Vector3d{a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)};
}
Eigen::VectorXd math::screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    Eigen::Vector3d v = math::vector_cross_product(-s,q) + h * s;
    Eigen::VectorXd screw(6);
    screw << s(0), s(1), s(2), v(0), v(1), v(2);
    return screw;
}
Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d R;
    R << tf(0,0), tf(0,1), tf(0,2),
    tf(1,0), tf(1,1), tf(1,2),
    tf(2,0), tf(2,1), tf(2,2);

    Eigen::Matrix3d skew_p = math::skew_symmetric(Eigen::Vector3d(tf(0,3), tf(1,3), tf(2,3)));
    Eigen::MatrixXd ad(6,6);

    Eigen::Matrix3d a = skew_p * R;

    ad << R(0,0), R(0,1), R(0,2), 0, 0, 0,
    R(1,0), R(1,1), R(1,2), 0, 0, 0,
    R(2,0), R(2,1), R(2,2), 0, 0, 0,
    a(0,0), a(0,1), a(0,2), R(0,0), R(0,1), R(0,2),
    a(1,0), a(1,1), a(1,2), R(1,0), R(1,1), R(1,2),
    a(2,0), a(2,1), a(2,2), R(2,0), R(2,1), R(2,2);
    return ad;;
}
double math::cot(double deg)
{
    double rad = math::deg_to_rad*deg;
    return std::cos(rad)/std::sin(rad);
}
//Assignment 2. Task 3: Matrix exponentials and logarithms
Eigen::Matrix3d math::matrix_exponential(const Eigen::Vector3d &w, double t_deg)
{
    double t_rad = math::deg_to_rad*t_deg;
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + std::sin(t_rad)*math::skew_symmetric(w) + (1-std::cos(t_rad))*math::skew_symmetric(w)*math::skew_symmetric(w);
    return R;
}
std::pair<Eigen::Vector3d, double> math::matrix_logarithm(const Eigen::Matrix3d &r)
{
    double t_rad = 0;
    double t_deg = 0;
    Eigen::Vector3d w{0, 0, 0};

    if (r == Eigen::Matrix3d::Identity())
    {}
    else if (r(0,0)+r(1,1)+r(2,2) == -1)
    {
        t_rad = EIGEN_PI;
        t_deg = math::rad_to_deg*t_rad;

        Eigen::Vector3d w_1 {r(0,2), r(1,2), 1 + r(2,2)};
        w_1 = w_1*(1/sqrt(2*(1+r(2,2))));

        Eigen::Vector3d w_2 {r(0,1), 1+ r(1,1),r(2,1)};
        w_2 = w_2*(1/sqrt(2*(1+r(1,1))));

        Eigen::Vector3d w_3 {1 + r(0,0),r(1,0),r(2,0)};
        w_3 = w_3*(1/sqrt(2*(1+r(0,0))));

        if (matrix_exponential(w_1,t_deg) == r)
        {
            w = w_1;
        }
        else if (matrix_exponential(w_2,t_deg) == r)
        {
            w = w_2;
        }
        else if (matrix_exponential(w_3,t_deg) == r)
        {
            w = w_3;
        }
        else
        {
            t_rad = 0;
            t_deg = 0;
        }
    }
    else
    {
        t_rad = 1/std::cos(0.5*(r(0,0)+r(1,1)+r(2,2)-1));
        t_deg = math::rad_to_deg*t_rad;
        Eigen::Matrix3d skew_w = (r-r.transpose())/(2*std::sin(t_rad));
        w = math::vector_from_skew_symmetric(skew_w);
    }

    while (t_deg>=180) {t_deg-=360; t_rad-=2*EIGEN_PI;}
    while (t_deg<=-180) {t_deg+=360; t_rad+=2*EIGEN_PI;}

    return {w, t_deg};
}
Eigen::Matrix4d math::matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double t_deg)
{
    double t_rad = math::deg_to_rad*t_deg;
    Eigen::Matrix3d R = matrix_exponential(w, t_deg);
    Eigen::Vector3d Gv = (Eigen::Matrix3d::Identity()*t_rad + (1-std::cos(t_rad))*math::skew_symmetric(w) + (t_rad-std::sin(t_rad))*math::skew_symmetric(w)*math::skew_symmetric(w))*v;

    Eigen::Matrix4d T = math::transformation_matrix(R,Gv);
    return T;
}
std::pair<Eigen::VectorXd, double> math::matrix_logarithm(const Eigen::Matrix4d &t)
{
    Eigen::Matrix3d R;
    R << t(0,0), t(0,1), t(0,2),
    t(1,0), t(1,1), t(1,2),
    t(2,0), t(2,1), t(2,2);

    Eigen::Vector3d p{t(0,3), t(1,3), t(2,3)};

    Eigen::Vector3d w{0,0,0};
    Eigen::Vector3d v{0, 0, 0};
    double t_rad = 0;
    double t_deg = 0;

    if (R == Eigen::Matrix3d::Identity())
    {
        v = p/sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
        t_rad = sqrt(p(0)*p(0)+p(1)*p(1)+p(2)*p(2));
        t_deg = math::rad_to_deg*t_rad;
    }
    else
    {
        std::pair<Eigen::VectorXd, double> a = math::matrix_logarithm(R);
        w=a.first;
        t_deg=a.second;
        t_rad = math::deg_to_rad*t_deg;
        v = (Eigen::Matrix3d::Identity()/t_rad - 0.5*math::skew_symmetric(w) + (1/t_rad - 0.5*cot(math::rad_to_deg*t_rad/2))*math::skew_symmetric(w)*math::skew_symmetric(w))*p;
    }

    while (t_deg>=180) {t_deg-=360; t_rad-=2*EIGEN_PI;}
    while (t_deg<=-180) {t_deg+=360; t_rad+=2*EIGEN_PI;}

    Eigen::VectorXd s = twist(w,v);
    return {s, t_deg};
}
//Assignment 2. Task 4: Forward kinematics: 3R planar open chain
void math::print_pose(const std::string &label, const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d R;
    R << tf(0,0), tf(0,1), tf(0,2),
    tf(1,0), tf(1,1), tf(1,2),
    tf(2,0), tf(2,1), tf(2,2);

    Eigen::Vector3d p{tf(0,3), tf(1,3), tf(2,3)};

    Eigen::Vector3d e_zyx = math::euler_zyx_from_rotation_matrix(R);

    std::cout << "Euler ZYX Angles and Position for ";
    std::cout << label;
    std::cout << " pose:"<< std::endl << std::endl;
    std::cout << "Euler ZYX Angles:" << std::endl;
    std::cout << "Alpha = ";
    std::cout << e_zyx(0) << std::endl;
    std::cout << "Beta = ";
    std::cout << e_zyx(1) << std::endl;
    std::cout << "Gamma = ";
    std::cout << e_zyx(2) << std::endl << std::endl;
    std::cout << "Linear position vector:" << std::endl;
    std::cout << "p = ";
    std::cout << p.transpose() << std::endl << std::endl;

}
Eigen::Matrix4d math::planar_3r_fk_transform(const std::vector<double> &link_L,const std::vector<double> &joint_positions)
{
    Eigen::Vector3d w{0,0,1};
    Eigen::Vector3d p_01{0,0,0};
    Eigen::Vector3d p_12{link_L[0],0,0};
    Eigen::Vector3d p_23{link_L[1],0,0};
    Eigen::Vector3d p_34{link_L[2],0,0};

    Eigen::Matrix3d R_01 = math::matrix_exponential(w,joint_positions[0]);
    Eigen::Matrix3d R_12 = math::matrix_exponential(w,joint_positions[1]);
    Eigen::Matrix3d R_23 = math::matrix_exponential(w,joint_positions[2]);
    Eigen::Matrix3d R_34 = math::matrix_exponential(w,0);

    Eigen::Matrix4d T_01 = math::transformation_matrix(R_01,p_01);
    Eigen::Matrix4d T_12 = math::transformation_matrix(R_12,p_12);
    Eigen::Matrix4d T_23 = math::transformation_matrix(R_23,p_23);
    Eigen::Matrix4d T_34 = math::transformation_matrix(R_34,p_34);

    Eigen::Matrix4d T = T_01*T_12*T_23*T_34;
    return T;
}
Eigen::Matrix4d math::planar_3r_fk_screw(const std::vector<double> &link_L,const std::vector<double> &joint_positions)
{
    Eigen::Vector3d w{0,0,1};
    Eigen::Vector3d q_01{0,0,0};
    Eigen::Vector3d q_02{link_L[0],0,0};
    Eigen::Vector3d q_03{link_L[0]+link_L[1],0,0};
    Eigen::Vector3d q_04{link_L[0]+link_L[1]+link_L[2],0,0};

    Eigen::Matrix4d M = math::transformation_matrix(Eigen::Matrix3d::Identity(),q_04);

    Eigen::VectorXd s_03 = math::screw_axis(q_03,w,0);
    Eigen::Vector3d v_03{s_03(3), s_03(4), s_03(5)};
    Eigen::Matrix4d S_03 = math::matrix_exponential(w,v_03,joint_positions[2]);

    Eigen::VectorXd s_02 = math::screw_axis(q_02,w,0);
    Eigen::Vector3d v_02{s_02(3), s_02(4), s_02(5)};
    Eigen::Matrix4d S_02 = math::matrix_exponential(w,v_02,joint_positions[1]);

    Eigen::VectorXd s_01 = math::screw_axis(q_01,w,0);
    Eigen::Vector3d v_01{s_01(3), s_01(4), s_01(5)};
    Eigen::Matrix4d S_01 = math::matrix_exponential(w,v_01,joint_positions[0]);

    Eigen::Matrix4d T = S_01*S_02*S_03*M;
    return T;

}
//Assignment 2. Task 5: Forward kinematics: UR3e 6R open chain
Eigen::Matrix4d math::ur3e_fk_screw(const std::vector<double> &joint_positions)
{
    //Parameters for UR3E
    //To consider the 180 degree rotation of the {s} frame of the robot in class,
    //the following values and vectors should be changed in the sign (put a minus):
    //values: L_1, L_2, W_1, W_2; vectors: w_2, w_3, w_4, w_6
    double L_1 = 0.24355; double W_1 = 0.13105; double H_1 = 0.15185; //meters
    double L_2 = 0.2132; double W_2 = 0.0921; double H_2 = 0.08535; //meters

    double t_01 = joint_positions[0]; double t_12 = joint_positions[1]; double t_23 = joint_positions[2];
    double t_34 = joint_positions[3]; double t_45 = joint_positions[4]; double t_56 = joint_positions[5];

    Eigen::Vector3d w_1{0,0,1};
    Eigen::Vector3d w_2{0,1,0};
    Eigen::Vector3d w_3{0,1,0};
    Eigen::Vector3d w_4{0,1,0};
    Eigen::Vector3d w_5{0,0,-1};
    Eigen::Vector3d w_6{0,1,0};

    Eigen::Vector3d q_01{0,0,H_1};
    Eigen::Vector3d q_02{0,W_1,H_1};
    Eigen::Vector3d q_03{L_1,W_1,H_1};
    Eigen::Vector3d q_04{L_1+L_2,0,H_1};
    Eigen::Vector3d q_05{L_1+L_2,W_1,H_1};
    Eigen::Vector3d q_06{L_1+L_2,W_1+W_2,H_1-H_2};

    Eigen::Vector3d z_sb{0,1,0};
    Eigen::Vector3d y_sb{0,0,1};
    Eigen::Vector3d x_sb{-1,0,0};
    Eigen::Matrix3d R_sb;
    R_sb << x_sb(0), y_sb(0), z_sb(0),
    x_sb(1), y_sb(1), z_sb(1),
    x_sb(2), y_sb(2), z_sb(2);
    Eigen::Vector3d q_sb{L_1+L_2,W_1+W_2,H_1-H_2};
    Eigen::Matrix4d M = math::transformation_matrix(R_sb,q_sb);

    Eigen::VectorXd s_06 = math::screw_axis(q_06,w_6,0);
    Eigen::Vector3d v_06{s_06(3), s_06(4), s_06(5)};
    Eigen::Matrix4d S_06 = math::matrix_exponential(w_6,v_06,t_56);

    Eigen::VectorXd s_05 = math::screw_axis(q_05,w_5,0);
    Eigen::Vector3d v_05{s_05(3), s_05(4), s_05(5)};
    Eigen::Matrix4d S_05 = math::matrix_exponential(w_5,v_05,t_45);

    Eigen::VectorXd s_04 = math::screw_axis(q_04,w_4,0);
    Eigen::Vector3d v_04{s_04(3), s_04(4), s_04(5)};
    Eigen::Matrix4d S_04 = math::matrix_exponential(w_4,v_04,t_34);

    Eigen::VectorXd s_03 = math::screw_axis(q_03,w_3,0);
    Eigen::Vector3d v_03{s_03(3), s_03(4), s_03(5)};
    Eigen::Matrix4d S_03 = math::matrix_exponential(w_3,v_03,t_23);

    Eigen::VectorXd s_02 = math::screw_axis(q_02,w_2,0);
    Eigen::Vector3d v_02{s_02(3), s_02(4), s_02(5)};
    Eigen::Matrix4d S_02 = math::matrix_exponential(w_2,v_02,t_12);

    Eigen::VectorXd s_01 = math::screw_axis(q_01,w_1,0);
    Eigen::Vector3d v_01{s_01(3), s_01(4), s_01(5)};
    Eigen::Matrix4d S_01 = math::matrix_exponential(w_1,v_01,t_01);

    Eigen::Matrix4d T = S_01*S_02*S_03*S_04*S_05*S_06*M;
    return T;
}
Eigen::Matrix4d math::ur3e_fk_transform(const std::vector<double> &joint_positions)
{

    //Parameters from UR3E
    double L_1 = 0.24355; double W_1 = 0.13105; double H_1 = 0.15185; //meters
    double L_2 = 0.2132; double W_2 = 0.0921; double H_2 = 0.08535; //meters

    double t_00 = 180; //180 degree rotated axes for the robot that we have in class
    //I do not use t_00 & it's relatives for the exercise
    double t_01 = joint_positions[0]; double t_12 = joint_positions[1]; double t_23 = joint_positions[2];
    double t_34 = joint_positions[3]; double t_45 = joint_positions[4]; double t_56 = joint_positions[5];

    Eigen::Vector3d w_0{0,0,1}; //Rotation axis for t_00
    Eigen::Vector3d w_1{0,0,1};
    Eigen::Vector3d w_2{0,1,0};
    Eigen::Vector3d w_3{0,1,0};
    Eigen::Vector3d w_4{0,1,0};
    Eigen::Vector3d w_5{0,0,-1};
    Eigen::Vector3d w_6{0,1,0};

    Eigen::Vector3d q_00{0,0,0}; //Position for t_00
    Eigen::Vector3d q_01{0,0,H_1};
    Eigen::Vector3d q_12{0,0,0};
    Eigen::Vector3d q_23{L_1,0,0};
    Eigen::Vector3d q_34{L_2,0,0};
    Eigen::Vector3d q_45{0,W_1,0};
    Eigen::Vector3d q_56{0,W_2,-H_2};
    Eigen::Vector3d q_6b{0,0,0};

    Eigen::Matrix3d R_00 = math::matrix_exponential(w_0,t_00); //Rotation matrix for t_00
    Eigen::Matrix3d R_01 = math::matrix_exponential(w_1,t_01);
    Eigen::Matrix3d R_12 = math::matrix_exponential(w_2,t_12);
    Eigen::Matrix3d R_23 = math::matrix_exponential(w_3,t_23);
    Eigen::Matrix3d R_34 = math::matrix_exponential(w_4,t_34);
    Eigen::Matrix3d R_45 = math::matrix_exponential(w_5,t_45);
    Eigen::Matrix3d R_56 = math::matrix_exponential(w_6,t_56);
    Eigen::Matrix3d R_6b;
    R_6b << -1, 0, 0,
    0, 0, 1,
    0, 1, 0;

    Eigen::Matrix4d T_00 = math::transformation_matrix(R_00,q_00); //Transformation matrix for t_00
    Eigen::Matrix4d T_01 = math::transformation_matrix(R_01,q_01);
    Eigen::Matrix4d T_12 = math::transformation_matrix(R_12,q_12);
    Eigen::Matrix4d T_23 = math::transformation_matrix(R_23,q_23);
    Eigen::Matrix4d T_34 = math::transformation_matrix(R_34,q_34);
    Eigen::Matrix4d T_45 = math::transformation_matrix(R_45,q_45);
    Eigen::Matrix4d T_56 = math::transformation_matrix(R_56,q_56);
    Eigen::Matrix4d T_6b = math::transformation_matrix(R_6b,q_6b);

    //To consider the initial 180 degree rotation for the robot in class,
    //T_00 has to be added at the begining of the multiplication:
    //Eigen::Matrix4d T = T_00*T_01*T_12*T_23*T_34*T_45*T_56*T_6b; //With the initial 180 degree rotation
    Eigen::Matrix4d T = T_01*T_12*T_23*T_34*T_45*T_56*T_6b;
    return T;
}

//Assignment 3:
//Assignment 3. Task 1: Functions and algorithms
Eigen::VectorXd math::std_vector_to_eigen(const std::vector<double> &v) {
    Eigen::VectorXd r(v.size());
    for (int i = 0; i < v.size(); i++)
        r[i]=v[i];
    return r;
}
bool math::is_average_below_eps(const std::vector<double> &values, double eps, uint8_t n_values) {
    if (values.size() < n_values)
        return false;
    else {
        double average = 0.0;
        for (int i = (values.size()-n_values); i < values.size(); i++) {average=average+values[i];}
        average = average/n_values;
        if (average > eps)
            return false;
        else
            return true;
    }
}
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> math::ur3e_space_chain() {
    double L_1 = 0.24355; double L_2 = 0.2132; //meters
    double W_1 = 0.13105 ; double W_2 = 0.0921; //meters
    double H_1 = 0.15185; double H_2 = 0.08535; //meters

    Eigen::Vector3d w_1{0,0,1};
    Eigen::Vector3d w_2{0,1,0};
    Eigen::Vector3d w_3{0,1,0};
    Eigen::Vector3d w_4{0,1,0};
    Eigen::Vector3d w_5{0,0,-1};
    Eigen::Vector3d w_6{0,1,0};

    Eigen::Vector3d q_01{0,0,H_1};
    Eigen::Vector3d q_02{0,W_1,H_1};
    Eigen::Vector3d q_03{L_1,W_1,H_1};
    Eigen::Vector3d q_04{L_1+L_2,0,H_1};
    Eigen::Vector3d q_05{L_1+L_2,W_1,H_1};
    Eigen::Vector3d q_06{L_1+L_2,W_1+W_2,H_1-H_2};

    Eigen::Vector3d z_sb{0,1,0};
    Eigen::Vector3d y_sb{0,0,1};
    Eigen::Vector3d x_sb{-1,0,0};
    Eigen::Matrix3d R_sb;
    R_sb << x_sb(0), y_sb(0), z_sb(0),
    x_sb(1), y_sb(1), z_sb(1),
    x_sb(2), y_sb(2), z_sb(2);
    Eigen::Vector3d q_sb{L_1+L_2,W_1+W_2,H_1-H_2};
    Eigen::Matrix4d M = math::transformation_matrix(R_sb,q_sb);

    Eigen::VectorXd s_01 = math::screw_axis(q_01,w_1,0);
    Eigen::VectorXd s_02 = math::screw_axis(q_02,w_2,0);
    Eigen::VectorXd s_03 = math::screw_axis(q_03,w_3,0);
    Eigen::VectorXd s_04 = math::screw_axis(q_04,w_4,0);
    Eigen::VectorXd s_05 = math::screw_axis(q_05,w_5,0);
    Eigen::VectorXd s_06 = math::screw_axis(q_06,w_6,0);

    std::vector<Eigen::VectorXd> s={s_01,s_02,s_03,s_04,s_05,s_06};

    return {M,s};
}
Eigen::Matrix4d math::ur3e_space_fk(const Eigen::VectorXd &joint_positions) {
    Eigen::Matrix4d M = math::ur3e_space_chain().first;
    std::vector<Eigen::VectorXd> s = math::ur3e_space_chain().second;
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    for (int i=0; i<s.size(); i++) {
        double t = math::rad_to_deg*joint_positions[i];
        Eigen::Vector3d w{s[i](0), s[i](1), s[i](2)};
        Eigen::Vector3d v{s[i](3), s[i](4), s[i](5)};
        Eigen::Matrix4d S = math::matrix_exponential(w,v,t);
        T=T*S;
    }
    T=T*M;
    return T;
}
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> math::ur3e_body_chain() {
    Eigen::Matrix4d M_s = math::ur3e_space_chain().first;
    std::vector<Eigen::VectorXd> s_s = math::ur3e_space_chain().second;
    Eigen::Matrix4d M_b=M_s.inverse();
    Eigen::MatrixXd A = math::adjoint_matrix(M_b);
    std::vector<Eigen::VectorXd> s_b(s_s.size());
    for (int i=0; i<s_s.size(); i++) {s_b[i]= A*s_s[i];}
    return {M_b,s_b};
}
Eigen::Matrix4d math::ur3e_body_fk(const Eigen::VectorXd &joint_positions) {
    Eigen::Matrix4d M_s = math::ur3e_space_chain().first;
    std::vector<Eigen::VectorXd> s_b = math::ur3e_body_chain().second;
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    for (int i=0; i<s_b.size(); i++) {
        double t = math::rad_to_deg*joint_positions[i];
        Eigen::Vector3d w{s_b[i](0), s_b[i](1), s_b[i](2)};
        Eigen::Vector3d v{s_b[i](3), s_b[i](4), s_b[i](5)};
        Eigen::Matrix4d S = math::matrix_exponential(w,v,t);
        T=T*S;
    }
    T=M_s*T;
    return T;
}
//Assignment 3. Task 2: Numerical optimization
std::pair<uint32_t, double> math::newton_raphson_root_find(const std::function<double(double)> &f,
    double x_0, double dx_0, double eps) {
    std::vector<double> x = {x_0};
    uint32_t i=0;
    double error=eps+1;
    while (error>eps) {
        x.push_back(x[i] - ((2*dx_0) / (f(x[i]+dx_0) - f(x[i]-dx_0))) * f(x[i]));
        i++;
        error = abs(f(x[i])-f(x[i-1]));
        //error = abs(f(x[i])-0);
    }
    return {i, x[i]};
}
std::pair<uint32_t, double> math::gradient_descent_root_find(const std::function<double(double)> &f,
    double x_0, double gamma, double dx_0, double eps) {
    std::vector<double> x = {x_0};
    uint32_t i=0;
    x.push_back(x[i] - gamma*(f(x[i]+dx_0) - f(x[i]))/(dx_0));
    std::vector<double> error= {abs(f(x[i])-0),abs(f(x[i+1])-0)};
    i++;
    while (error[i]>eps) {
        x.push_back(x[i] - gamma*(error[i] - error[i-1])/(x[i]-x[i-1]));
        i++;
        //error = abs(f(x[i])-f(x[i-1]));
        error.push_back(abs(f(x[i])-0));
    }
    return {i, x[i]};
}
//Assignment 3. Task 3: Velocity kinematics: UR3e 6R open chain
Eigen::MatrixXd math::ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions) {
    auto [M,s] = math::ur3e_space_chain();
    int n = current_joint_positions.size();
    Eigen::MatrixXd J(n,n);
    Eigen::Matrix4d E = Eigen::Matrix4d::Identity();
    J.col(0)=s[0];
    for (int i=1; i<n; i++) {
        E = E*math::matrix_exponential({s[i-1](0),s[i-1](1),s[i-1](2)},
                {s[i-1](3),s[i-1](4),s[i-1](5)},math::rad_to_deg*current_joint_positions[i-1]);
        J.col(i) = math::adjoint_matrix(E)*s[i];
    }
    return J;
}
Eigen::MatrixXd math::ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions) {
    auto [M,s] = math::ur3e_body_chain();
    int n = current_joint_positions.size();
    Eigen::MatrixXd J(n,n);
    Eigen::Matrix4d E = Eigen::Matrix4d::Identity();
    J.col(n-1)=s[n-1];
    for (int i=n-2; i>=0; i--) {
        E = E*(math::matrix_exponential({s[i+1](0),s[i+1](1),s[i+1](2)},
        {s[i+1](3),s[i+1](4),s[i+1](5)},
            math::rad_to_deg*current_joint_positions(i+1)).inverse());
        J.col(i) = math::adjoint_matrix(E)*s[i];
    }
    return J;
}
//Assignment 3. Task 4: Inverse kinematics: UR3e 6R open chain
std::pair<int, Eigen::VectorXd> math::ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd
&current_joint_positions, double gamma, double v_e, double w_e) {
    std::pair<Eigen::VectorXd,double> pair;
    Eigen::Matrix4d V;
    Eigen::VectorXd theta = current_joint_positions;
    Eigen::VectorXd twist(6);
    int i=1;

    while ((twist.head(3).norm() > w_e or twist.tail(3).norm() > v_e) and i<1e4) {
        V = math::ur3e_space_fk(theta).inverse()*t_sd;
        pair=math::matrix_logarithm(V);
        twist = pair.first*pair.second*math::deg_to_rad;
        theta += gamma * math::ur3e_body_jacobian(theta).completeOrthogonalDecomposition().pseudoInverse() * twist;
        i++;
    }

    for (int j=0; j<theta.size(); j++){
        while (theta(j)>=EIGEN_PI) {theta(j)-=2*EIGEN_PI;}
        while (theta(j)<=-EIGEN_PI) {theta(j)+=2*EIGEN_PI;}
    }

    return {i, theta};
}
