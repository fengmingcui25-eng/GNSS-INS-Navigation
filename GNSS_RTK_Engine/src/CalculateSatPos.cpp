#include"RTK.h"
#include <Eigen\Dense>
using namespace Eigen;
using namespace std;

// 求解开普勒方程
double solveKeplerEquation(double M_t, double e)
{
    double EPSILON = 1e-12;
    double E_t = M_t;
    double E_t_next;
    do
    {
        E_t_next = M_t + e * sin(E_t);
        if (fabs(E_t_next - E_t) < EPSILON)
        {
            break;
        }
        E_t = E_t_next;
    } while (true);

    return E_t;
}

bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATMIDRES* Mid)
{
    // 根据卫星系统选择星历指针
    GPSEPHREC* eph = nullptr;
    if (Sys == GPS) {
        eph = GPSEph;
    }
    else if (Sys == BDS) {
        eph = BDSEph;
    }
    else {
        return false; 
    }
    // 检查星历有效性及PRN匹配
    if (!eph || eph->PRN != Prn || eph->System != Sys) {
        return false;
    }
    // 检查卫星健康状态
    if (eph->SVHealth != 0) {
        return false;
    }
    // 计算当前时间绝对秒数（GPS周秒）
    const double t_gps_sec = t->Week * 604800.0 + t->SecOfWeek;
    double t_sys_sec = t_gps_sec;
    // BDS系统需转换到BDT时间
    if (Sys == BDS) t_sys_sec -= (14.0 + 1356 * 604800);
    // 计算星历参考时间（TOE）的绝对秒数
    const double toe_sec = eph->TOE.Week * 604800.0 + eph->TOE.SecOfWeek;
    // 判断星历是否过期
    const double max_valid_interval = (Sys == GPS) ? 7500.0 : 3900.0; // 对更新时间进行设置GPS 2h, BDS 1h
    if (fabs(t_sys_sec - toe_sec) > max_valid_interval) return false;
    
    // 计算钟差参考时间（TOC）的绝对秒数
    const double toc_sec = eph->TOC.Week * 604800.0 + eph->TOC.SecOfWeek;
    const double dt = t_sys_sec - toc_sec; // 相对于钟差参考时间的时间差
    // 计算钟差和钟速
    Mid->SatClkOft = eph->ClkBias + eph->ClkDrift * dt + eph->ClkDriftRate * dt * dt;
    Mid->SatClkSft = eph->ClkDrift + 2 * eph->ClkDriftRate * dt ;

    // 群延差
    Mid->Tgd1 = eph->TGD1;
    Mid->Tgd2 = eph->TGD2;

    // 标记计算结果有效
    Mid->Valid = true;
    return true;
}

int CompGPSSatPVT(const int Prn,  GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid)
{
    double A = pow(Eph->SqrtA, 2);
    double n0;
    n0 = sqrt(u_GPS) / (pow(A, 1.5));
    GPSTIME t_k;
    t_k.Week = t->Week - Eph->TOE.Week;
    t_k.SecOfWeek = t->SecOfWeek - Eph->TOE.SecOfWeek;
    double tk = t_k.Week * 604800.0 + t_k.SecOfWeek;
    /*if (tk > 302400.0) tk -= 604800.0;
    else if (tk < -302400.0) tk += 604800.0;*/
    double n;
    n = n0 + Eph->DetlaN;

    double M_t; // 观测瞬间t时刻卫星的平近点角
    M_t = Eph->M0 + n * tk;
    double E_t;
    E_t = solveKeplerEquation(M_t, Eph->e);
    double f_t; // 真近点角
    f_t = atan2(sqrt(1 - Eph->e * Eph->e) * sin(E_t), (cos(E_t) - Eph->e));

    double u_t_dot;
    u_t_dot = Eph->omega + f_t;
    double r_t_dot;
    r_t_dot = A * (1 - Eph->e * cos(E_t));
    double u_t, r_t, i_t;
    u_t = u_t_dot + Eph->Cuc * cos(2 * u_t_dot) + Eph->Cus * sin(2 * u_t_dot);
    r_t = r_t_dot + Eph->Crc * cos(2 * u_t_dot) + Eph->Crs * sin(2 * u_t_dot);
    i_t = Eph->i0 + Eph->iDot * tk + Eph->Cic * cos(2 * u_t_dot) + Eph->Cis * sin(2 * u_t_dot);
    double x_t, y_t;
    x_t = r_t * cos(u_t);
    y_t = r_t * sin(u_t);
    double L_t;
    L_t = Eph->OMEGA + (Eph->OMEGADot - omega_e) * tk - omega_e * Eph->TOE.SecOfWeek;//Eph->OMEGADot
    
    Mid->SatPos[0] = x_t * cos(L_t) - y_t * cos(i_t) * sin(L_t);
    Mid->SatPos[1] = x_t * sin(L_t) + y_t * cos(i_t) * cos(L_t);
    Mid->SatPos[2] = y_t * sin(i_t);
    double EkDOT = n / (1 - Eph->e * cos(E_t));
    double faikDOT = sqrt(1 - Eph->e * Eph->e) * EkDOT / (1 - Eph->e * cos(E_t));
    double ukDOT = 2 * (Eph->Cus * cos(2 * u_t_dot) - Eph->Cuc * sin(2 * u_t_dot)) * faikDOT + faikDOT;
    double rkDOT = A * Eph->e * sin(E_t) * EkDOT + 2 * (Eph->Crs * cos(2 * u_t_dot) - Eph->Crc * sin(2 * u_t_dot)) * faikDOT;
    double IkDOT=Eph->iDot+ 2 * (Eph->Cis * cos(2 * u_t_dot) - Eph->Cic * sin(2 * u_t_dot)) * faikDOT;
    double OMEGAkDOT = Eph->OMEGADot - omega_e;
    MatrixXd RDOT(3, 4);
    RDOT << cos(L_t), -sin(L_t) * cos(i_t), -(x_t * sin(L_t) + y_t * cos(L_t) * cos(i_t)), y_t* sin(L_t)* sin(i_t),
        sin(L_t), cos(L_t)* cos(i_t), (x_t * cos(L_t) - y_t * sin(L_t) * cos(i_t)), -y_t * cos(L_t) * sin(i_t),
        0, sin(i_t), 0, y_t* cos(i_t);
    double xkDOT = rkDOT * cos(u_t) - r_t * ukDOT * sin(u_t);
    double ykDOT = rkDOT * sin(u_t) + r_t * ukDOT * cos(u_t);
    MatrixXd v(3, 1);
    MatrixXd mat(4, 1);
    mat << xkDOT,
        ykDOT,
        OMEGAkDOT,
        IkDOT;
    v = RDOT * mat;
    Mid->SatVel[0] = v(0, 0);
    Mid->SatVel[1] = v(1, 0);
    Mid->SatVel[2] = v(2, 0);

    Mid->SatClkOft = F * Eph->e * Eph->SqrtA * sin(E_t);
    Mid->SatClkSft = F * Eph->e * Eph->SqrtA * cos(E_t) * EkDOT;
    return 0;
}

int CompBDSSatPVT(const int Prn, GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid)
{
    t->Week = t->Week - 1356;
    t->SecOfWeek = t->SecOfWeek - 14.0;
    double A = pow(Eph->SqrtA, 2);
    double n0;
    n0 = sqrt(u_BDS) / (pow(A, 1.5));
    GPSTIME t_k;
    t_k.Week = t->Week - Eph->TOE.Week;
    t_k.SecOfWeek = t->SecOfWeek - Eph->TOE.SecOfWeek;
    double tk = t_k.Week * 604800.0 + t_k.SecOfWeek;
    double n;
    n = n0 + Eph->DetlaN;

    double M_t; // 观测瞬间t时刻卫星的平近点角
    M_t = Eph->M0 + n * tk;
    double E_t;
    E_t = solveKeplerEquation(M_t, Eph->e);
    double f_t; // 真近点角
    f_t = atan2(sqrt(1 - Eph->e * Eph->e) * sin(E_t), (cos(E_t) - Eph->e));

    double u_t_dot;
    u_t_dot = Eph->omega + f_t;
    double r_t_dot;
    r_t_dot = A * (1 - Eph->e * cos(E_t));
    double u_t, r_t, i_t;
    u_t = u_t_dot + Eph->Cuc * cos(2 * u_t_dot) + Eph->Cus * sin(2 * u_t_dot);
    r_t = r_t_dot + Eph->Crc * cos(2 * u_t_dot) + Eph->Crs * sin(2 * u_t_dot);
    i_t = Eph->i0 + Eph->iDot * tk + Eph->Cic * cos(2 * u_t_dot) + Eph->Cis * sin(2 * u_t_dot);
    double x_t, y_t;
    x_t = r_t * cos(u_t);
    y_t = r_t * sin(u_t);
    double L_t;
    MatrixXd Rx(3, 3), Rz(3, 3), X(3, 1), x(3, 1);
    Rx << 1, 0, 0,//GEO矩阵
        0, cos(-5.0 / 180.0 * pi), sin(-5.0 / 180.0 * pi),
        0, -sin(-5.0 / 180.0 * pi), cos(-5.0 / 180.0 * pi);
    Rz << cos(Omega_e_dot_BDS * tk), sin(Omega_e_dot_BDS * tk), 0,
        -sin(Omega_e_dot_BDS * tk), cos(Omega_e_dot_BDS * tk), 0,
        0, 0, 1;
    double EkDOT = n / (1 - Eph->e * cos(E_t));
    if ((Prn >= 1 && Prn <= 5) || (Prn >= 59 && Prn <= 63)/*Eph->i0 < 30.0 / 180.0 * pi*/)//GEO卫星
    {
        L_t = Eph->OMEGA + Eph->OMEGADot * tk - Omega_e_dot_BDS * Eph->TOE.SecOfWeek;
        X(0, 0) = x_t * cos(L_t) - y_t * cos(i_t) * sin(L_t);
        X(1, 0) = x_t * sin(L_t) + y_t * cos(i_t) * cos(L_t);
        X(2, 0) = y_t * sin(i_t);
        
        x = Rz * Rx * X;
        Mid->SatPos[0] = x(0, 0);
        Mid->SatPos[1] = x(1, 0);
        Mid->SatPos[2] = x(2, 0);

        double faikDOT = sqrt(1 - Eph->e * Eph->e) * EkDOT / (1 - Eph->e * cos(E_t));
        double ukDOT = 2 * (Eph->Cus * cos(2 * u_t_dot) - Eph->Cuc * sin(2 * u_t_dot)) * faikDOT + faikDOT;
        double rkDOT = A * Eph->e * sin(E_t) * EkDOT + 2 * (Eph->Crs * cos(2 * u_t_dot) - Eph->Crc * sin(2 * u_t_dot)) * faikDOT;
        double IkDOT = Eph->iDot + 2 * (Eph->Cis * cos(2 * u_t_dot) - Eph->Cic * sin(2 * u_t_dot)) * faikDOT;
        double OMEGAkDOT = Eph->OMEGADot;
        MatrixXd RDOT(3, 4);
        RDOT << cos(L_t), -sin(L_t) * cos(i_t), -(x_t * sin(L_t) + y_t * cos(L_t) * cos(i_t)), y_t* sin(L_t)* sin(i_t),
            sin(L_t), cos(L_t)* cos(i_t), (x_t * cos(L_t) - y_t * sin(L_t) * cos(i_t)), -y_t * cos(L_t) * sin(i_t),
            0, sin(i_t), 0, y_t* cos(i_t);
        double xkDOT = rkDOT * cos(u_t) - r_t * ukDOT * sin(u_t);
        double ykDOT = rkDOT * sin(u_t) + r_t * ukDOT * cos(u_t);
        MatrixXd v(3,1),v1(3, 1), v2(3, 1);
        MatrixXd mat(4, 1);
        mat << xkDOT,
            ykDOT,
            OMEGAkDOT,
            IkDOT;
        v1 = Rz * Rx * RDOT * mat;
        MatrixXd RzDOT(3, 3);
        RzDOT << -Omega_e_dot_BDS *sin(Omega_e_dot_BDS * tk), Omega_e_dot_BDS* cos(Omega_e_dot_BDS * tk), 0,
            -Omega_e_dot_BDS * cos(Omega_e_dot_BDS * tk), -Omega_e_dot_BDS * sin(Omega_e_dot_BDS * tk), 0,
            0, 0, 0;
        v2 = RzDOT * Rx * X;
        v = v1 + v2;
        Mid->SatVel[0] = v(0, 0);
        Mid->SatVel[1] = v(1, 0);
        Mid->SatVel[2] = v(2, 0);
    }
    else if(Prn > 5 && Prn < 59)//MEO/IGSO
    {
        L_t = Eph->OMEGA + (Eph->OMEGADot - Omega_e_dot_BDS) * tk - Omega_e_dot_BDS * Eph->TOE.SecOfWeek;
        Mid->SatPos[0] = x_t * cos(L_t) - y_t * cos(i_t) * sin(L_t);
        Mid->SatPos[1] = x_t * sin(L_t) + y_t * cos(i_t) * cos(L_t);
        Mid->SatPos[2] = y_t * sin(i_t);

        double faikDOT = sqrt(1 - Eph->e * Eph->e) * EkDOT / (1 - Eph->e * cos(E_t));
        double ukDOT = 2 * (Eph->Cus * cos(2 * u_t_dot) - Eph->Cuc * sin(2 * u_t_dot)) * faikDOT + faikDOT;
        double rkDOT = A * Eph->e * sin(E_t) * EkDOT + 2 * (Eph->Crs * cos(2 * u_t_dot) - Eph->Crc * sin(2 * u_t_dot)) * faikDOT;
        double IkDOT = Eph->iDot + 2 * (Eph->Cis * cos(2 * u_t_dot) - Eph->Cic * sin(2 * u_t_dot)) * faikDOT;
        double OMEGAkDOT = Eph->OMEGADot - Omega_e_dot_BDS;
        MatrixXd RDOT(3, 4);
        RDOT << cos(L_t), -sin(L_t) * cos(i_t), -(x_t * sin(L_t) + y_t * cos(L_t) * cos(i_t)), y_t* sin(L_t)* sin(i_t),
            sin(L_t), cos(L_t)* cos(i_t), (x_t * cos(L_t) - y_t * sin(L_t) * cos(i_t)), -y_t * cos(L_t) * sin(i_t),
            0, sin(i_t), 0, y_t* cos(i_t);
        double xkDOT = rkDOT * cos(u_t) - r_t * ukDOT * sin(u_t);
        double ykDOT = rkDOT * sin(u_t) + r_t * ukDOT * cos(u_t);
        MatrixXd v(3, 1);
        MatrixXd mat(4, 1);
        mat << xkDOT,
            ykDOT,
            OMEGAkDOT,
            IkDOT;
        v = RDOT * mat;
        Mid->SatVel[0] = v(0, 0);
        Mid->SatVel[1] = v(1, 0);
        Mid->SatVel[2] = v(2, 0);
    }
    Mid->SatClkOft = F * Eph->e * Eph->SqrtA * sin(E_t);
    Mid->SatClkSft = F * Eph->e * Eph->SqrtA * cos(E_t) * EkDOT;

    return 0;
}

int CorrectEarthRotation(double Xr[], SATMIDRES* Mid)
{
    MatrixXd X(3, 1), XDOT(3, 1), Rz(3, 3);
    double L, delta_t;

    //for (int i = 0; i < 10; i++)
    //{
    //    L = sqrt(pow(Xr[0] - Mid->SatPos[0], 2) + pow(Xr[1] - Mid->SatPos[1], 2) + pow(Xr[2] - Mid->SatPos[2], 2));
    //    delta_t = L / c;


    //    X << Mid->SatPos[0], Mid->SatPos[1], Mid->SatPos[2];
    //    XDOT << Mid->SatVel[0], Mid->SatVel[1], Mid->SatVel[2];
    //    Rz << cos(omega_e * delta_t), sin(omega_e * delta_t), 0,
    //        -sin(omega_e * delta_t), cos(omega_e * delta_t), 0,
    //        0, 0, 1;
    //    X = Rz * X;
    //    XDOT = Rz * XDOT;
    //}
    L = sqrt(pow(Xr[0] - Mid->SatPos[0], 2) + pow(Xr[1] - Mid->SatPos[1], 2) + pow(Xr[2] - Mid->SatPos[2], 2));
    delta_t = L / C_Light;


    X << Mid->SatPos[0], Mid->SatPos[1], Mid->SatPos[2];
    XDOT << Mid->SatVel[0], Mid->SatVel[1], Mid->SatVel[2];
    Rz << cos(omega_e * delta_t), sin(omega_e * delta_t), 0,
        -sin(omega_e * delta_t), cos(omega_e * delta_t), 0,
        0, 0, 1;
    X = Rz * X;
    XDOT = Rz * XDOT;

    Mid->SatPos[0] = X(0, 0);
    Mid->SatPos[1] = X(1, 0);
    Mid->SatPos[2] = X(2, 0);

    Mid->SatVel[0] = XDOT(0, 0);
    Mid->SatVel[1] = XDOT(1, 0);
    Mid->SatVel[2] = XDOT(2, 0);

    return 0;
}