#include <iostream>

#include <Eigen/Dense>

#include "math/math.h"

//Task 1: Functions and algorithms
//T.1 a) Calculate an Euler ZYX vector from a rotation matrix.
Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
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

    double a_deg = math::rad_to_deg(a_rad); //Z axis Alpha
    double b_deg = math::rad_to_deg(b_rad); //Y axis Beta
    double c_deg = math::rad_to_deg(c_rad); //X axis Gamma

    return Eigen::Vector3d{a_deg, b_deg, c_deg};
}

//T.1 b) Create a twist vector from velocity vectors.
Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::VectorXd tw(6);
    tw << w(0), w(1), w(2), v(0), v(1), v(2);
    return tw;
}

//T.1 c) Create a screw axis.
//Define cross product for two vectors
Eigen::Vector3d vector_cross_product(const Eigen::Vector3d &a, const Eigen::Vector3d &b)
{
    return Eigen::Vector3d{a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)};
}
//Complete task 1 c)
Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    Eigen::Vector3d v = vector_cross_product(-s,q) + h * s;
    Eigen::VectorXd screw(6);
    screw << s(0), s(1), s(2), v(0), v(1), v(2);
    return screw;
}

//T.1 d) Create the adjoint representation of a homogeneous transformation matrix.
Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf)
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

//T.1 e) Calculate the cotangent.
double cot(double deg)
{
    double rad = math::deg_to_rad(deg);
    return std::cos(rad)/std::sin(rad);
}

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

    Eigen::VectorXd F_h = twist(w_h,f_h);
    Eigen::VectorXd F_a = twist(w_a,f_a);

    Eigen::Matrix3d rot_hf = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d rot_af;
    rot_af << 1, 0, 0,
    0, 0, 1,
    0, -1, 0;

    Eigen::Vector3d p_hf {-L_1, 0, 0};
    Eigen::Vector3d p_af {-(L_1+L_2), 0, 0};

    Eigen::Matrix4d T_hf = math::transformation_matrix(rot_hf,p_hf);
    Eigen::Matrix4d T_af = math::transformation_matrix(rot_af,p_af);

    Eigen::MatrixXd Ad_hf = adjoint_matrix(T_hf);
    Eigen::MatrixXd Ad_af = adjoint_matrix(T_af);

    Eigen::VectorXd F_f = Ad_hf.transpose()*F_h + Ad_af.transpose()*F_a;

    std::cout << "F_f:";
    std::cout << F_f.transpose() << std::endl << std::endl;

}

//Task 3: Matrix exponentials and logarithms
//T.3 a) Implement the matrix exponential for rotation matrices.

Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double t_deg)
{
    double t_rad = math::deg_to_rad(t_deg);
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + std::sin(t_rad)*math::skew_symmetric(w) + (1-std::cos(t_rad))*math::skew_symmetric(w)*math::skew_symmetric(w);
    return R;
}

//T.3 b) Implement the matrix logarithm for rotation matrices.
std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r)
{
    double t_rad = 0;
    double t_deg = 0;
    Eigen::Vector3d w{0, 0, 0};

    if (r == Eigen::Matrix3d::Identity())
    {}
    else if (r(0,0)+r(1,1)+r(2,2) == -1)
    {
        t_rad = EIGEN_PI;
        t_deg = math::rad_to_deg(t_rad);

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
        t_deg = math::rad_to_deg(t_rad);
        Eigen::Matrix3d skew_w = (r-r.transpose())/(2*std::sin(t_rad));
        w = math::vector_from_skew_symmetric(skew_w);

    }
    return std::pair(w, t_deg);
}

//T.3 c) Implement the matrix exponential for homogeneous transformation matrices.
Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double t_deg)
{
    double t_rad = math::deg_to_rad(t_deg);
    Eigen::Matrix3d R = matrix_exponential(w, t_deg);
    Eigen::Vector3d Gv = (Eigen::Matrix3d::Identity()*t_rad + (1-std::cos(t_rad))*math::skew_symmetric(w) + (t_rad-std::sin(t_rad))*math::skew_symmetric(w)*math::skew_symmetric(w))*v;

    Eigen::Matrix4d T = math::transformation_matrix(R,Gv);
    return T;
}

//T.3 d) Implement the matrix logarithm for homogeneous transformation matrices.
std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t)
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
        t_deg = math::rad_to_deg(t_rad);
    }
    else
    {
        std::pair (w,t_deg) = matrix_logarithm(R);
        t_rad = math::deg_to_rad(t_deg);
        v = (Eigen::Matrix3d::Identity()/t_rad - 0.5*math::skew_symmetric(w) + (1/t_rad - 0.5*cot(math::rad_to_deg(t_rad/2)))*math::skew_symmetric(w)*math::skew_symmetric(w))*p;
    }

    Eigen::VectorXd s = twist(w,v);
    return std::pair(s, t_deg);
}

//Task 4: Forward kinematics: 3R planar open chain
//T.4 a) Print the Euler ZYX angles and linear position of a homogeneous transformation matrix.
void print_pose(const std::string &label, const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d R;
    R << tf(0,0), tf(0,1), tf(0,2),
    tf(1,0), tf(1,1), tf(1,2),
    tf(2,0), tf(2,1), tf(2,2);

    Eigen::Vector3d p{tf(0,3), tf(1,3), tf(2,3)};

    Eigen::Vector3d e_zyx = euler_zyx_from_rotation_matrix(R);

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

//T.4 b) Calculate forward kinematics using transformation matrices.
Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions)
{
    double t_01 = joint_positions[0];
    double t_12 = joint_positions[1];
    double t_23 = joint_positions[2];
    double t_34 = 0;
    double L_1 = 10;
    double L_2 = 10;
    double L_3 = 10;
    Eigen::Vector3d w{0,0,1};
    Eigen::Vector3d p_01{0,0,0};
    Eigen::Vector3d p_12{L_1,0,0};
    Eigen::Vector3d p_23{L_2,0,0};
    Eigen::Vector3d p_34{L_3,0,0};

    Eigen::Matrix3d R_01 = matrix_exponential(w,t_01);
    Eigen::Matrix3d R_12 = matrix_exponential(w,t_12);
    Eigen::Matrix3d R_23 = matrix_exponential(w,t_23);
    Eigen::Matrix3d R_34 = matrix_exponential(w,t_34);

    Eigen::Matrix4d T_01 = math::transformation_matrix(R_01,p_01);
    Eigen::Matrix4d T_12 = math::transformation_matrix(R_12,p_12);
    Eigen::Matrix4d T_23 = math::transformation_matrix(R_23,p_23);
    Eigen::Matrix4d T_34 = math::transformation_matrix(R_34,p_34);

    Eigen::Matrix4d T = T_01*T_12*T_23*T_34;
    return T;
}

//T.4 c) Calculate forward kinematics using the product of exponentials.
Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions)
{
    double t_01 = joint_positions[0];
    double t_12 = joint_positions[1];
    double t_23 = joint_positions[2];
    double t_34 = 0;
    double L_1 = 10;
    double L_2 = 10;
    double L_3 = 10;

    Eigen::Vector3d w{0,0,1};
    Eigen::Vector3d q_01{0,0,0};
    Eigen::Vector3d q_02{L_1,0,0};
    Eigen::Vector3d q_03{L_1+L_2,0,0};
    Eigen::Vector3d q_04{L_1+L_2+L_3,0,0};

    Eigen::Matrix4d M = math::transformation_matrix(Eigen::Matrix3d::Identity(),q_04);

    Eigen::VectorXd s_03 = screw_axis(q_03,w,0);
    Eigen::Vector3d v_03{s_03(3), s_03(4), s_03(5)};
    Eigen::Matrix4d S_03 = matrix_exponential(w,v_03,t_23);

    Eigen::VectorXd s_02 = screw_axis(q_02,w,0);
    Eigen::Vector3d v_02{s_02(3), s_02(4), s_02(5)};
    Eigen::Matrix4d S_02 = matrix_exponential(w,v_02,t_12);

    Eigen::VectorXd s_01 = screw_axis(q_01,w,0);
    Eigen::Vector3d v_01{s_01(3), s_01(4), s_01(5)};
    Eigen::Matrix4d S_01 = matrix_exponential(w,v_01,t_01);

    Eigen::Matrix4d T = S_01*S_02*S_03*M;
    return T;

}

//Task 5: Forward kinematics: UR3e 6R open chain
//T.5 a) Calculate forward kinematics using the product of exponentials.
Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions)
{
    //Parameters from the book. Example 4.5
    /*
    double L_1 = 0.425; //meters
    double L_2 = 0.392; //meters
    double W_1 = 0.109; //meters
    double W_2 = 0.082; //meters
    double H_1 = 0.089; //meters
    double H_2 = 0.095; //meters
    */

    //Parameters for UR3E
    //To consider the 180 degree rotation of the {s} frame of the robot in class,
    //the following values and vectors should be changed in the sign (put a minus):
    //values: L_1, L_2, W_1, W_2; vectors: w_2, w_3, w_4, w_6
    double L_1 = 0.24355; //meters
    double L_2 = 0.2132; //meters
    double W_1 = 0.13105 ; //meters
    double W_2 = 0.0921; //meters
    double H_1 = 0.15185; //meters
    double H_2 = 0.08535; //meters

    double t_01 = joint_positions[0];
    double t_12 = joint_positions[1];
    double t_23 = joint_positions[2];
    double t_34 = joint_positions[3];
    double t_45 = joint_positions[4];
    double t_56 = joint_positions[5];

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

    Eigen::VectorXd s_06 = screw_axis(q_06,w_6,0);
    Eigen::Vector3d v_06{s_06(3), s_06(4), s_06(5)};
    Eigen::Matrix4d S_06 = matrix_exponential(w_6,v_06,t_56);

    Eigen::VectorXd s_05 = screw_axis(q_05,w_5,0);
    Eigen::Vector3d v_05{s_05(3), s_05(4), s_05(5)};
    Eigen::Matrix4d S_05 = matrix_exponential(w_5,v_05,t_45);

    Eigen::VectorXd s_04 = screw_axis(q_04,w_4,0);
    Eigen::Vector3d v_04{s_04(3), s_04(4), s_04(5)};
    Eigen::Matrix4d S_04 = matrix_exponential(w_4,v_04,t_34);

    Eigen::VectorXd s_03 = screw_axis(q_03,w_3,0);
    Eigen::Vector3d v_03{s_03(3), s_03(4), s_03(5)};
    Eigen::Matrix4d S_03 = matrix_exponential(w_3,v_03,t_23);

    Eigen::VectorXd s_02 = screw_axis(q_02,w_2,0);
    Eigen::Vector3d v_02{s_02(3), s_02(4), s_02(5)};
    Eigen::Matrix4d S_02 = matrix_exponential(w_2,v_02,t_12);

    Eigen::VectorXd s_01 = screw_axis(q_01,w_1,0);
    Eigen::Vector3d v_01{s_01(3), s_01(4), s_01(5)};
    Eigen::Matrix4d S_01 = matrix_exponential(w_1,v_01,t_01);

    Eigen::Matrix4d T = S_01*S_02*S_03*S_04*S_05*S_06*M;
    return T;
}

//T.5 b) Calculate forward kinematics using homogeneous transformation matrices.
Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions)
{
    //Parameters from the book. Example 4.5
    /*double L_1 = 0.425; //meters
    double L_2 = 0.392; //meters
    double W_1 = 0.109; //meters
    double W_2 = 0.082; //meters
    double H_1 = 0.089; //meters
    double H_2 = 0.095; //meters
    */

    //Parameters from UR3E
    double L_1 = 0.24355; //meters
    double L_2 = 0.2132; //meters
    double W_1 = 0.13105 ; //meters
    double W_2 = 0.0921; //meters
    double H_1 = 0.15185; //meters
    double H_2 = 0.08535; //meters

    double t_00 = 180; //180 degree rotated axes for the robot that we have in class
    //I do not use t_00 & it's relatives for the exercise
    double t_01 = joint_positions[0];
    double t_12 = joint_positions[1];
    double t_23 = joint_positions[2];
    double t_34 = joint_positions[3];
    double t_45 = joint_positions[4];
    double t_56 = joint_positions[5];

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

    Eigen::Matrix3d R_00 = matrix_exponential(w_0,t_00); //Rotation matrix for t_00
    Eigen::Matrix3d R_01 = matrix_exponential(w_1,t_01);
    Eigen::Matrix3d R_12 = matrix_exponential(w_2,t_12);
    Eigen::Matrix3d R_23 = matrix_exponential(w_3,t_23);
    Eigen::Matrix3d R_34 = matrix_exponential(w_4,t_34);
    Eigen::Matrix3d R_45 = matrix_exponential(w_5,t_45);
    Eigen::Matrix3d R_56 = matrix_exponential(w_6,t_56);
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

int main() {
    //Task 2 a)
    task2a();

    //Task 2 b)
    task2b();

    //Task 4 a) & b)
    std::vector<double> j4_1 = {0,0,0};
    Eigen::Matrix4d T4_1 = planar_3r_fk_transform(j4_1);
    std::wcout << "Configuration 1: ";
    for (double i: j4_1){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 1 Using D-H", T4_1);
    Eigen::Matrix4d T4s_1 = planar_3r_fk_screw(j4_1);
    print_pose("Configuration 1 Using PoE", T4s_1);

    std::vector<double> j4_2 = {90,0,0};
    Eigen::Matrix4d T4_2 = planar_3r_fk_transform(j4_2);
    std::wcout << "Configuration 2: ";
    for (double i: j4_2){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 2 Using D-H", T4_2);
    Eigen::Matrix4d T4s_2 = planar_3r_fk_screw(j4_2);
    print_pose("Configuration 2 Using PoE", T4s_2);

    std::vector<double> j4_3 = {0,90,0};
    Eigen::Matrix4d T4_3 = planar_3r_fk_transform(j4_3);
    std::wcout << "Configuration 3: ";
    for (double i: j4_3){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 3 Using D-H", T4_3);
    Eigen::Matrix4d T4s_3 = planar_3r_fk_screw(j4_3);
    print_pose("Configuration 3 Using PoE", T4s_3);

    std::vector<double> j4_4 = {0,0,90};
    Eigen::Matrix4d T4_4 = planar_3r_fk_transform(j4_4);
    std::wcout << "Configuration 4: ";
    for (double i: j4_4){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 4 Using D-H", T4_4);
    Eigen::Matrix4d T4s_4 = planar_3r_fk_screw(j4_4);
    print_pose("Configuration 4 Using PoE", T4s_4);

    std::vector<double> j4_5 = {10,-15,2.75};
    Eigen::Matrix4d T4_5 = planar_3r_fk_transform(j4_5);
    std::wcout << "Configuration 5: ";
    for (double i: j4_5){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 5 Using D-H", T4_5);
    Eigen::Matrix4d T4s_5 = planar_3r_fk_screw(j4_5);
    print_pose("Configuration 5 Using PoE", T4s_5);


    //Task 5 a) & b)
    std::vector<double> j_example = {0,-90,0,0,90,0};
    Eigen::Matrix4d T_example = ur3e_fk_transform(j_example);
    std::wcout << "Configuration Example: ";
    for (double i: j_example){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration Example Using D-H", T_example);
    Eigen::Matrix4d Ts_example = ur3e_fk_screw(j_example);
    print_pose("Configuration Example Using PoE", Ts_example);

    std::vector<double> j_1 = {0,0,0,-90,0,0};
    Eigen::Matrix4d T_1 = ur3e_fk_transform(j_1);
    std::wcout << "Configuration 1: ";
    for (double i: j_1){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 1 Using D-H", T_1);
    Eigen::Matrix4d Ts_1 = ur3e_fk_screw(j_1);
    print_pose("Configuration 1 Using PoE", Ts_1);

    std::vector<double> j_2 = {0,-180,0,0,0,0};
    Eigen::Matrix4d T_2 = ur3e_fk_transform(j_2);
    std::wcout << "Configuration 2: ";
    for (double i: j_2){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 2 Using D-H", T_2);
    Eigen::Matrix4d Ts_2 = ur3e_fk_screw(j_2);
    print_pose("Configuration 2 Using PoE", Ts_2);

    std::vector<double> j_3 = {0,-90,0,0,0,0};
    Eigen::Matrix4d T_3 = ur3e_fk_transform(j_3);
    std::wcout << "Configuration 3: ";
    for (double i: j_3){std::cout << i;std::cout << ",";}std::cout << std::endl;
    print_pose("Configuration 3 Using D-H", T_3);
    Eigen::Matrix4d Ts_3 = ur3e_fk_screw(j_3);
    print_pose("Configuration 3 Using PoE", Ts_3);

    return 0;
}