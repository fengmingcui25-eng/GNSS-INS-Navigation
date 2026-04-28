#include <iostream>
#include <cstring> 
#include<fstream>
#include <iomanip>
#include "INS.h"

using namespace std;

int main() {
    FILE* FIMUData = nullptr;
    int lenR_imu;        // 实际接收长度
    int lenD_imu = 0;    // 缓冲区中剩余数据长度
    unsigned char buff_imu[MAXRAWLEN]; // 数据缓冲区
    IMUDATA obs;
    IMUDATA obs_last;
    Vector3d acc_bias_incr(0.0, 0.0, 0.0);  // 加速度计零偏增量
    Vector3d gyro_bias_incr(0.0, 0.0, 0.0);
    InsState state;
    double acc[3] = { 0.0 }, gyro[3] = { 0.0 };//前五分钟数据均值
    double Init_Atti[3] = { 0.0 };//roll;pitch;yaw,对准结果
    string outFile = "results_no_test.txt";
    ofstream fout(outFile);
    if (!fout.is_open())
    {
        cerr << "Error: Failed to create " << outFile << endl;
        return -1;
    }
    fout << fixed << setprecision(9);
    fout << "Time\tLat\tLon\tH\tVn\tVe\tVd\tRoll\tPitch\tYaw" << endl;

    cout << "开始惯导解算..." << endl;
    int idx_time = 0;

    //初始化
#if WorkMode == 1
    //**********************************采集数据***************************************
    double epochs = 442.0 * SampleRate;//前五分钟历元总数
    Vector3d start_Pos(30.5278839567 * DEG_TO_RAD, 114.3557289981 * DEG_TO_RAD, 19.701);/*初始位置30.5278839567;114.3557289981;19.701*/
    Vector3d start_Vel(0.0, 0.0, 0.0);/*初始速度*/
    Vector3d start_Att(0.6075357981, 0.2310624196, 5.0069750392);/*初始姿态5.0069750392   0.2310624196   0.6075357981,单位：°head;pitch;roll*/
    string filename = "group5.ASC";
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << endl;
        return 0;
    }
    string line;

    while (getline(file, line)) {
        // 只处理以 %RAWIMUSA 开头的行
        if (line.find("%RAWIMUSA") != 0) {
            continue; // 跳过其他类型的数据
        }

        //解码
        obs = parseIMULine(line);
        /******************************************************************************/
        //staticCoarseAlignmentENU(obs.accel, obs.gyro, Latitude, Init_Atti[0], Init_Atti[1], Init_Atti[2]);//得到结果为角度
        //cout << "" << Init_Atti[0] << "\t" << Init_Atti[1] << "\t" << Init_Atti[2] << endl;

        /******************************************************************************/
        //开始五分钟用于对准、误差补偿
        if (idx_time < epochs)
        {
            for (int i = 0; i < 3; i++)
            {
                acc[i] += obs.accel[i];
                gyro[i] += obs.gyro[i];
            }
            idx_time++;
            continue;
        }
        //计算前五分钟数据均值
        if (fabs(idx_time - epochs) < 0.001)
        {
            for (int i = 0; i < 3; i++)
            {
                acc[i] = acc[i] / epochs ;//(m/s)
                gyro[i] = gyro[i] / epochs ;//(m/s)
            }

            //零偏估计方法
            Vector3d acc_mean_val(acc[0]/dt, acc[1]/dt, acc[2]/dt);  // 加速度计测量均值 (m/s^2)

            // 计算当地重力
            double g_local = GetGravity(start_Pos.x(), start_Pos.z());/******************(m/s^2)*******************/

            // 利用已知初始姿态计算理论比力
            Vector3d init_att_rad = start_Att * DEG_TO_RAD;
            Quaterniond q_init = Euler2Q(init_att_rad);
            Matrix3d C_b_n = Q2DCM(q_init);
            Matrix3d C_n_b = C_b_n.transpose();

            // 导航系下的理论比力
            Vector3d f_n_theory(0.0, 0.0, -g_local);

            // 投影到载体系
            Vector3d f_b_theory = C_n_b * f_n_theory;

            // 计算全轴加速度零偏
            Vector3d acc_bias_val = acc_mean_val - f_b_theory;

            // 转换为增量形式（假设采样间隔为Dt，需要根据实际设置）

            acc_bias_incr = acc_bias_val * dt;
            cout << "零偏估计结果-加速度计零偏(m /s^2) : " << acc_bias_val.transpose() << endl;

            //初始化
            state.pos = start_Pos;
            state.vel = start_Vel;
            // 直接使用已知初始姿态
            state.att = start_Att * DEG_TO_RAD;//单位：rad
            //使用对准结果
            //for (int i = 0; i < 3; i++)state.att(i) = Init_Atti[i];


            state.q = q_init;
            state.C = C_b_n;
            state.static_time = 0.0;
            obs_last = obs;

            idx_time++;
            continue;
        }

        Vector3d current_acc(obs.accel[0], obs.accel[1], obs.accel[2]);
        Vector3d current_gyro(obs.gyro[0], obs.gyro[1], obs.gyro[2]);

        // 加速度计零偏补偿
        current_acc = current_acc - acc_bias_incr;

        // 将补偿后的值赋回obs
        for (int i = 0; i < 3; i++) {
            obs.accel[i] = current_acc(i);
        }
        //输入数据
        //陀螺
        Vector3d d_theta;
        Vector3d d_theta_old;
        //加速度计
        Vector3d d_vel;
        Vector3d d_vel_old;
        for (int i = 0; i < 3; i++)
        {
            d_theta(i) = obs.gyro[i];
            d_vel(i) = obs.accel[i];
            d_theta_old(i) = obs_last.gyro[i];
            d_vel_old(i) = obs_last.accel[i];
        }
        // 机械编排


        bool zupt_triggered = false;

        // 保存上一时刻数据
        double vel_last[3];
        for (int i = 0; i < 3; i++) {
            vel_last[i] = state.vel(i);
        }

        // 1：速度更新（使用k-1时刻的姿态）
        VelUpdate(&state, obs_last.gyro, obs.gyro, obs_last.accel, obs.accel);

        // 2：位置更新
        PosUpdate(&state, vel_last);

        // 3：姿态更新
        AttitudeUpdate(&state, obs_last.gyro, obs.gyro);

        // 4：零速修正
        //ZuptCheck(obs.accel, obs.gyro, &state);

        obs_last = obs;

        double out_lat = RAD_TO_DEG * (state.pos.x());
        double out_lon = RAD_TO_DEG * (state.pos.y());
        double out_h = state.pos.z();

        double out_vn = state.vel.x();
        double out_ve = state.vel.y();
        double out_vd = state.vel.z();

        double out_roll = RAD_TO_DEG * (state.att.x());
        double out_pitch = RAD_TO_DEG * (state.att.y());
        double out_yaw = RAD_TO_DEG * (state.att.z());

        fout << obs.time << "\t"
            << out_lat << "\t" << out_lon << "\t" << out_h << "\t"
            << out_vn << "\t" << out_ve << "\t" << out_vd << "\t"
            << out_roll << "\t" << out_pitch << "\t" << out_yaw << endl;
    }
    file.close();

#else
    ///************************************读取示例数据**********************************/
    Vector3d start_Pos(23.1373950708 * DEG_TO_RAD, 113.3713651222 * DEG_TO_RAD, 2.175);/*初始位置*/
    Vector3d start_Vel(0.0, 0.0, 0.0);/*初始速度*/
    Vector3d start_Att(0.0107951084511778, -2.14251290749072, -75.7498049314083);/*初始姿态*/



    if ((FIMUData = fopen("IMU.bin", "rb")) == NULL) { // 确保以二进制模式打开[6,8](@ref)
        printf("Cannot open IMU data file.\n");
        return -1;
    }

    // 利用已知初始姿态计算理论比力
    Vector3d init_att_rad = start_Att * DEG_TO_RAD;
    Quaterniond q_init = Euler2Q(init_att_rad);
    Matrix3d C_b_n = Q2DCM(q_init);
    Matrix3d C_n_b = C_b_n.transpose();


    //初始化
    state.pos = start_Pos;
    state.vel = start_Vel;
    // 直接使用已知初始姿态
    state.att = start_Att * DEG_TO_RAD;//单位：rad
    //使用对准结果
    //for (int i = 0; i < 3; i++)state.att(i) = Init_Atti[i];


    state.q = q_init;
    state.C = C_b_n;
    state.static_time = 0.0;

    while (1)
    {
        // 读取数据到缓冲区空闲部分
        lenR_imu = fread(buff_imu + lenD_imu, sizeof(unsigned char), MAXRAWLEN - lenD_imu, FIMUData);

        if (lenR_imu <= 0) {
            if (feof(FIMUData)) {
                printf("Reached end of IMU data file.\n");
            }
            else {
                perror("Error reading IMU data file");
            }
            break; // 退出循环
        }

        lenD_imu += lenR_imu; // 更新缓冲区中数据总长度


        // 循环处理缓冲区中所有完整的数据帧


        int processed_len = 0;
        while (lenD_imu >= 56)
        {
            if (Decode_IMUDATA(buff_imu + processed_len, &obs) == 0) { // 假设解码成功返回0
                // 解码成功，处理obs数据
                processed_len += 56; // 更新已处理数据长度
                lenD_imu -= 56;    // 更新剩余数据长度

                /////////////误差补偿
                if (obs.time - 91620.0 < 0.0001)
                {
                    obs_last = obs;
                    continue;
                }

                Vector3d current_acc(obs.accel[0], obs.accel[1], obs.accel[2]);
                Vector3d current_gyro(obs.gyro[0], obs.gyro[1], obs.gyro[2]);

                Vector3d last_acc(obs_last.accel[0], obs_last.accel[1], obs_last.accel[2]);
                Vector3d last_gyro(obs_last.gyro[0], obs_last.gyro[1], obs_last.gyro[2]);

                //输入数据
                //陀螺
                Vector3d d_theta;
                Vector3d d_theta_old;
                //加速度计
                Vector3d d_vel;
                Vector3d d_vel_old;
                for (int i = 0; i < 3; i++)
                {
                    d_theta(i) = obs.gyro[i];
                    d_vel(i) = obs.accel[i];
                    d_theta_old(i) = obs_last.gyro[i];
                    d_vel_old(i) = obs_last.accel[i];
                }
                // 机械编排

                bool zupt_triggered = false;

                // 保存上一时刻数据
                double vel_last[3];
                for (int i = 0; i < 3; i++)
                    vel_last[i] = state.vel(i);

                // ✅ 第一步：速度更新（使用k-1时刻的姿态）
                VelUpdate(&state, obs_last.gyro, obs.gyro, obs_last.accel, obs.accel);

                // ✅ 第二步：位置更新
                PosUpdate(&state, vel_last);

                // ✅ 第三步：姿态更新
                AttitudeUpdate(&state, obs_last.gyro, obs.gyro);

                obs_last = obs;

                double out_lat = RAD_TO_DEG * (state.pos.x());
                double out_lon = RAD_TO_DEG * (state.pos.y());
                double out_h = state.pos.z();

                double out_vn = state.vel.x();
                double out_ve = state.vel.y();
                double out_vd = state.vel.z();

                double out_roll = RAD_TO_DEG * (state.att.x());
                double out_pitch = RAD_TO_DEG * (state.att.y());
                double out_yaw = RAD_TO_DEG * (state.att.z());

                fout << obs.time << "\t"
                    << out_lat << "\t" << out_lon << "\t" << out_h << "\t"
                    << out_vn << "\t" << out_ve << "\t" << out_vd << "\t"
                    << out_roll << "\t" << out_pitch << "\t" << out_yaw << endl;

            }
            else {
                // 解码失败，可能数据不完整，保留剩余数据到下一轮读取
                break;
            }
        }


        // 将缓冲区中未处理的数据移动到缓冲区开头
        if (processed_len > 0) {
            memmove(buff_imu, buff_imu + processed_len, lenD_imu);
        }

        // 防止缓冲区溢出
        if (lenD_imu >= MAXRAWLEN) {
            fprintf(stderr, "Buffer overflow risk. No valid frame found in %d bytes.\n", lenD_imu);
            break;
        }
    }

    if (FIMUData) fclose(FIMUData);

#endif
    fout.close();
    cout << "解算完成! 结果已保存至 " << outFile << endl;
    return 0;
}