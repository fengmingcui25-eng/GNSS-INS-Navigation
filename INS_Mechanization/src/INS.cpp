#include"INS.h"


Matrix3d Q2DCM(const Eigen::Quaterniond& q)
{
    double q0 = q.w(), q1 = q.x(), q2 = q.y(), q3 = q.z();
    Matrix3d C;
    C << q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2),
        2 * (q1 * q2 + q0 * q3), q0* q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1),
        2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0* q0 - q1 * q1 - q2 * q2 + q3 * q3;
    return C;
}

Quaterniond Euler2Q(const Eigen::Vector3d& a_)
{
    double hf = 0.5 * a_.x(), ht = 0.5 * a_.y(), hp = 0.5 * a_.z();
    double cf = cos(hf), sf = sin(hf), ct = cos(ht), st = sin(ht), cp = cos(hp), sp = sin(hp);
    return Quaterniond(cf * ct * cp + sf * st * sp, sf * ct * cp - cf * st * sp,
        cf * st * cp + sf * ct * sp, cf * ct * sp - sf * st * cp);
}

Vector3d DCM2Euler(const Eigen::Matrix3d& C)
{
    return Vector3d(atan2(C(2, 1), C(2, 2)),
        atan2(-C(2, 0), std::sqrt(C(2, 1) * C(2, 1) + C(2, 2) * C(2, 2))),
        atan2(C(1, 0), C(0, 0)));
}

// 辅助函数：计算向量的斜对称矩阵
Matrix3d skewSymmetric(const Vector3d& v) {
    Matrix3d S;
    S << 0, -v.z(), v.y(),
        v.z(), 0, -v.x(),
        -v.y(), v.x(), 0;
    return S;
}

//正常重力计算
double GetGravity(double lat, double h)
{
    double sin2lat = sin(lat) * sin(lat);
    double sin4lat = sin2lat * sin2lat;

    double g0 = GRAVITY * (1 + 0.0052790414 * sin2lat + 0.0000232718 * sin4lat);
    double g = g0 - (3.087691089e-6 - 4.397731e-9 * sin2lat) * h + 0.721e-12 * h * h;

    return g;
}

//计算有害加速度
static Vector3d GetAgc(const Vector3d& vel, const Vector3d& pos)
{
    double lat = pos.x(), h = pos.z();
    double sl = std::sin(lat), cl = std::cos(lat), tl = std::tan(lat);
    double sq = 1.0 - E2_WGS84 * sl * sl;
    double RN = R_WGS84 / std::sqrt(sq);
    double RM = R_WGS84 * (1.0 - E2_WGS84) / std::pow(sq, 1.5);

    // 重力向量
    Vector3d gn(0, 0, GetGravity(lat, h));

    // 自转角速度
    Vector3d wie_n(w_e * cl, 0, -w_e * sl);

    // 运输角速度
    Vector3d wen_n(vel.y() / (RN + h), -vel.x() / (RM + h), -vel.y() * tl / (RN + h));

    // agc
    return gn - (2.0 * wie_n + wen_n).cross(vel);
}
/// <summary>
/// /////注意输入数据！输入为速度增量（单位m/s）和角度增量（单位rad）
/// </summary>
/// <param name="state"></param>
/// <param name="gyro_last"></param>
/// <param name="gyro"></param>
void AttitudeUpdate(InsState* state,double gyro_last[],double gyro[])
{
    // STEP 1: 等效旋转矢量法更新 b 系（二子样算法）
   // φ_k = Δθ_k + (1/12) * Δθ_{k-1} × Δθ_k
    Vector3d d_theta;
    Vector3d d_theta_old;
    for (int i = 0; i < 3; i++)
    {
        d_theta(i) = gyro[i]; d_theta_old(i) = gyro_last[i];
    }
    Vector3d phi_k = d_theta + (1.0 / 12.0) * d_theta_old.cross(d_theta);

    // 旋转矢量转四元数：q = [cos(‖φ‖/2), sin(‖φ‖/2)*φ/‖φ‖]
    double phi_norm = phi_k.norm();
    Quaterniond q_b_step;
    if (phi_norm < 1e-12) {
        q_b_step = Quaterniond::Identity();
    }
    else {
        double half_phi = 0.5 * phi_norm;
        q_b_step = Quaterniond(cos(half_phi),
            sin(half_phi) * phi_k.x() / phi_norm,
            sin(half_phi) * phi_k.y() / phi_norm,
            sin(half_phi) * phi_k.z() / phi_norm);
    }

    // STEP 2: 等效旋转矢量法更新 n 系
    // 计算导航系旋转角速度：ω = ω_ie_n + ω_en_n
    
    double sl = sin(state->pos.x()), cl = cos(state->pos.x()), tl = tan(state->pos.x());
    double sq = 1.0 - E2_WGS84 * sl * sl;
    double RN = R_WGS84 / sqrt(sq);
    double RM = R_WGS84 * (1.0 - E2_WGS84) / pow(sq, 1.5);

    Vector3d w_ie_n(w_e * cl, 0.0, -w_e * sl);
    Vector3d w_en_n(state->vel.y() / (RN + state->pos.z()), -state->vel.x() / (RM + state->pos.z()), -state->vel.y() * tl / (RN + state->pos.z()));
    Vector3d w_in_n = w_ie_n + w_en_n;

    // ζ_k = ω * Δt
    Vector3d zeta_k = w_in_n * dt;// dt=1 / SampleRate;

    // 旋转矢量转四元数（注意负号，因为导航系旋转方向与载体相反）
    double zeta_norm = zeta_k.norm();
    Quaterniond q_n_step;
    if (zeta_norm < 1e-12) {
        q_n_step = Quaterniond::Identity();
    }
    else {
        double half_zeta = 0.5 * zeta_norm;
        Vector3d neg_zeta = -zeta_k; // 取负
        q_n_step = Quaterniond(cos(half_zeta),
            sin(half_zeta) * neg_zeta.x() / zeta_norm,
            sin(half_zeta) * neg_zeta.y() / zeta_norm,
            sin(half_zeta) * neg_zeta.z() / zeta_norm);
    }

    // STEP 3: 四元数连乘更新：q_b^n(k) = q_n(k-1)^n(k) * q_b^n(k-1) * q_b(k)^b(k-1)
    state->q = q_n_step * state->q * q_b_step;

    // STEP 4: 归一化处理
    state->q.normalize();

    // 更新方向余弦矩阵和欧拉角
    state->C = Q2DCM(state->q);
    state->att = DCM2Euler(state->C);
}


void VelUpdate(InsState* state, double gyro_last[], double gyro[], double acc_last[], double acc[])
{
    // 1. 比力积分项（包含划桨效应补偿）
   // dv_fb = dv_k + 0.5 * dθ_k × dv_k + (1/12) * (dθ_{k-1} × dv_k + dv_{k-1} × dθ_k)
    
    Vector3d d_theta;
    Vector3d d_theta_old;
    Vector3d d_vel;
    Vector3d d_vel_old;
    for (int i = 0; i < 3; i++)
    {
        d_theta(i) = gyro[i]; d_theta_old(i) = gyro_last[i];
        d_vel(i) = acc[i]; d_vel_old(i) = acc_last[i];
    }

    Vector3d dv_scull = d_vel + 0.5 * d_theta.cross(d_vel) +
        (1.0 / 12.0) * (d_theta_old.cross(d_vel) + d_vel_old.cross(d_theta));

    // 2. 导航系旋转补偿
    double lat = state->pos.x(), h = state->pos.z();
    double sl = sin(lat), cl = cos(lat);
    double sq = 1.0 - E2_WGS84 * sl * sl;
    double RN = R_WGS84 / sqrt(sq);
    double RM = R_WGS84 * (1.0 - E2_WGS84) / pow(sq, 1.5);

    Vector3d w_ie_n(w_e * cl, 0.0, -w_e * sl);
    Vector3d w_en_n(state->vel.y() / (RN + h), -state->vel.x() / (RM + h),
        -state->vel.y() * tan(lat) / (RN + h));
    Vector3d zeta = (w_ie_n + w_en_n) * dt;// dt=1 / SampleRate;

    // 旋转补偿：[I - 0.5 * ζ×] * C * dv_fb
    Matrix3d I_minus_half_zeta = Matrix3d::Identity() - 0.5 * skewSymmetric(zeta);
    Vector3d dv_f_n = I_minus_half_zeta * state->C * dv_scull;

    // 3. 有害加速度补偿（梯形积分法）
    // 使用k-1时刻的有害加速度
    Vector3d agc_old = GetAgc(state->vel, state->pos);

    // 预测k时刻的速度（用于计算k时刻的有害加速度）
    Vector3d vel_pred = state->vel + dv_f_n + agc_old * dt;// dt=1 / SampleRate;
    Vector3d agc_new = GetAgc(vel_pred, state->pos);

    // 梯形积分：dv_gc = 0.5 * (agc_old + agc_new) * dt
    Vector3d dv_gc = 0.5 * (agc_old + agc_new) * dt;// dt=1 / SampleRate;

    // 4. 速度更新：v_k = v_{k-1} + dv_f_n + dv_gc
    state->vel += dv_f_n + dv_gc;
}

void PosUpdate(InsState* state,double vel_last[])
{
    Vector3d vel_old;
    for (int i = 0; i < 3; i++)
        vel_old(i) = vel_last[i];
    // 提取上一时刻的位置和速度
    double lat_old = state->pos.x();
    double lon_old = state->pos.y();
    double h_old = state->pos.z();

    Vector3d v_curr = state->vel;

    // 1. 高程更新（梯形积分法）
    // h_k = h_{k-1} - 0.5 * (v_D,k-1 + v_D,k) * Δt
    double h_new = h_old - 0.5 * (vel_old.z() + v_curr.z()) * dt;
    double h_bar = 0.5 * (h_old + h_new);  // 平均高程

    // 2. 纬度更新（梯形积分法）
    // 计算k-1时刻的子午圈半径
    double sl_old = sin(lat_old);
    double sq_old = 1.0 - E2_WGS84 * sl_old * sl_old;
    double RM_old = R_WGS84 * (1.0 - E2_WGS84) / pow(sq_old, 1.5);

    // φ_k = φ_{k-1} + (v_N,k-1 + v_N,k) / [2*(R_M + h̄)] * Δt
    double lat_new = lat_old + (0.5 * (vel_old.x() + v_curr.x()) * dt) / (RM_old + h_bar);

    // 计算平均纬度
    double lat_bar = 0.5 * (lat_old + lat_new);

    // 3. 经度更新（梯形积分法）
    // 计算平均纬度处的卯酉圈半径
    double sl_bar = sin(lat_bar);
    double sq_bar = 1.0 - E2_WGS84 * sl_bar * sl_bar;
    double RN_bar = R_WGS84 / sqrt(sq_bar);

    // λ_k = λ_{k-1} + (v_E,k-1 + v_E,k) / [2*(R_N + h̄)*cos(φ̄)] * Δt
    double lon_new = lon_old + (0.5 * (vel_old.y() + v_curr.y()) * dt) /
        ((RN_bar + h_bar) * cos(lat_bar));

    // 4. 更新位置向量
    //state->pos << lat_new, lon_new, h_new;
    state->pos(0) = lat_new;
    state->pos(1) = lon_new;
    state->pos(2) = h_new;
}

//零速修正实现
void ZuptCheck(double acc[], double gyro[], InsState* state)
{
    // 1. 调用 GetGravity 计算当地精确重力值
    Vector3d d_vel(acc[0],acc[1],acc[2]);
    Vector3d d_theta(gyro[0],gyro[1],gyro[2]);
    double g_local = GetGravity(state->pos.x(), state->pos.z());

    // ZUPT 阈值设置
    const double TH_ACC = 0.2;          // 加速度模长容差 (m/s^2)
    const double TH_GYRO = 0.05;         // 角速度模长容差 (rad/s)
    const double TIME_TO_LOCK = 10;     // 判定静止所需的最短持续时间 (s)

    // 2. 计算传感器数据的模长
    double acc_norm = d_vel.norm() / dt;//(m/s^2)
    double gyro_norm = d_theta.norm() / dt;

    // 3. 判定当前时刻是否满足静止条件
    bool is_static = (abs(acc_norm - g_local) < TH_ACC) && (gyro_norm < TH_GYRO);

    // 4. 计时器逻辑
    if (is_static)
    {
        state->static_time += dt;
    }
    else
    {
        state->static_time = 0.0;
    }

    // 5. 若静止时间达标，执行零速修正
    if (state->static_time >= TIME_TO_LOCK)
    {
        state->vel.setZero();
    }
}