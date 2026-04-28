#include"INS.h"
#include <fstream>


using namespace std;
using namespace Eigen;

// 常量定义



//const double acc_scale = 1.0E-9;
//const double gyo_scale = 2.0E-8;
//const int SampleRate = 200;


// 矢量归一化函数
Vector3d vectorNormalize(const Vector3d& vec) {
    double norm = vec.norm();
    if (norm > 1e-10) {
        return vec / norm;
    }
    return vec;
}

void CalAverage(const vector<IMUDATA>& imu_data,double avg_acc[],double ave_gyro[])
{
    double sum_acc[3] = { 0.0 }, sum_gyro[3] = { 0.0 };
    int count = imu_data.size();

    for (const auto& data : imu_data) {
        for (int i = 0; i < 3; i++)
        {
            sum_acc[i] += data.accel[i];
            sum_gyro[i] += data.gyro[i];
        }
    }
    for (int i = 0; i < 3; i++)
    {
        avg_acc[i] = sum_acc[i] / count;
        ave_gyro[i] = sum_gyro[i] / count;
    }
}


// 静态解析粗对准函数NED
void staticCoarseAlignmentNED(double acc[],double gyro[], double latitude_deg, double& roll, double& pitch, double& yaw)
{
    // 转换为弧度
    double latitude = latitude_deg * DEG_TO_RAD;
    // 北东地坐标系下的理论矢量
    Vector3d g_n(0, 0, -GRAVITY);  // 重力向下

    // 地球自转矢量：北向和天向分量
    Vector3d wie_n(w_e * cos(latitude),
        0,
        -w_e * sin(latitude));

    // 载体系下的测量矢量
    Vector3d g_b = Vector3d::Zero(); //avg_acc;      // 加速度计测量
    Vector3d wie_b = Vector3d::Zero(); //= avg_gyro;   // 陀螺仪测量

    for (int i = 0; i < 3; i++)
    {
        g_b(i) = acc[i];
        wie_b(i) = gyro[i];
    }

   /* cout << "理论重力矢量: " << g_n.transpose() << endl;
    cout << "理论地球自转矢量: " << wie_n.transpose() << endl;*/

    // 导航系下的正交归一化基向量
    Vector3d v_g, v_omega, v_gomega;

    // v_g = g_n / |g_n|
    v_g = vectorNormalize(g_n);

    // v_omega = (g_n × wie_n) / |g_n × wie_n|
    Vector3d temp1 = g_n.cross(wie_n);
    v_omega = vectorNormalize(temp1);

    // v_gomega = (g_n × wie_n × g_n) / |g_n × wie_n × g_n|
    Vector3d temp2 = temp1.cross(g_n);
    v_gomega = vectorNormalize(temp2);

    //载体系下的正交归一化基向量
    Vector3d w_g, w_omega, w_gomega;

    // w_g = g_b / |g_b|
    w_g = vectorNormalize(g_b);

    // w_omega = (g_b × wie_b) / |g_b × wie_b|
    Vector3d temp3 = g_b.cross(wie_b);
    w_omega = vectorNormalize(temp3);

    // w_gomega = (g_b × wie_b × g_b) / |g_b × wie_b × g_b|
    Vector3d temp4 = temp3.cross(g_b);
    w_gomega = vectorNormalize(temp4);

    // 构造导航系和载体系的基向量矩阵
    Matrix3d V, W;
    V << v_g, v_omega, v_gomega;
    W << w_g, w_omega, w_gomega;

    // 计算姿态矩阵 C_bn = V * W^T
    Matrix3d C_bn = V * W.transpose();

    //return C_bn;
     // 北东地坐标系的欧拉角提取公式
    roll = atan2(-C_bn(2, 1), C_bn(2, 2));  // 横滚角
    pitch = -asin(C_bn(2, 0));               // 俯仰角
    yaw = atan2(C_bn(1, 0), C_bn(0, 0));    // 航向角

    // 转换为角度
    roll *= RAD_TO_DEG;
    pitch *= RAD_TO_DEG;
    yaw *= RAD_TO_DEG;

    // 确保航向角在0-360度范围内
    if (yaw < 0) yaw += 360.0;
}


// 静态解析粗对准函数ENU
void staticCoarseAlignmentENU(double acc[], double gyro[], double latitude_deg, double& roll, double& pitch, double& yaw)
{
    // 转换为弧度
    double latitude = latitude_deg * DEG_TO_RAD;
    // 北东地坐标系下的理论矢量
    Vector3d g_n(0, 0, GRAVITY);  // 重力向下

    // 地球自转矢量：北向和天向分量
    Vector3d wie_n(0,
        w_e * cos(latitude),
        w_e * sin(latitude));

    // 载体系下的测量矢量
    Vector3d g_b = Vector3d::Zero(); //avg_acc;      // 加速度计测量
    Vector3d wie_b = Vector3d::Zero(); //= avg_gyro;   // 陀螺仪测量

    for (int i = 0; i < 3; i++)
    {
        g_b(i) = acc[i];
        wie_b(i) = gyro[i];
    }

    /* cout << "理论重力矢量: " << g_n.transpose() << endl;
     cout << "理论地球自转矢量: " << wie_n.transpose() << endl;*/

     // 导航系下的正交归一化基向量
    Vector3d v_g, v_omega, v_gomega;

    // v_g = g_n / |g_n|
    v_g = vectorNormalize(g_n);

    // v_omega = (g_n × wie_n) / |g_n × wie_n|
    Vector3d temp1 = g_n.cross(wie_n);
    v_omega = vectorNormalize(temp1);

    // v_gomega = (g_n × wie_n × g_n) / |g_n × wie_n × g_n|
    Vector3d temp2 = temp1.cross(g_n);
    v_gomega = vectorNormalize(temp2);

    //载体系下的正交归一化基向量
    Vector3d w_g, w_omega, w_gomega;

    // w_g = g_b / |g_b|
    w_g = vectorNormalize(g_b);

    // w_omega = (g_b × wie_b) / |g_b × wie_b|
    Vector3d temp3 = g_b.cross(wie_b);
    w_omega = vectorNormalize(temp3);

    // w_gomega = (g_b × wie_b × g_b) / |g_b × wie_b × g_b|
    Vector3d temp4 = temp3.cross(g_b);
    w_gomega = vectorNormalize(temp4);

    // 构造导航系和载体系的基向量矩阵
    Matrix3d V, W;
    V << v_g, v_omega, v_gomega;
    W << w_g, w_omega, w_gomega;

    // 计算姿态矩阵 C_bn = V * W^T
    Matrix3d C_bn = V * W.transpose();

    //return C_bn;
     // 北东地坐标系的欧拉角提取公式
    roll = atan2(C_bn(1, 2), C_bn(2, 2));  // 横滚角
    pitch = -asin(C_bn(0, 2));               // 俯仰角
    yaw = atan2(C_bn(0, 1), C_bn(0, 0));    // 航向角

    // 转换为角度
    roll *= RAD_TO_DEG;
    pitch *= RAD_TO_DEG;
    yaw *= RAD_TO_DEG;

    // 确保航向角在0-360度范围内
    if (yaw < 0) yaw += 360.0;
}
//// 从姿态矩阵提取欧拉角（北东地坐标系）
//void extractEulerAnglesNED(const Matrix3d& C_bn, double& roll, double& pitch, double& yaw) {
//    // 北东地坐标系的欧拉角提取公式
//    roll = atan2(-C_bn(2, 1), C_bn(2, 2));  // 横滚角
//    pitch = -asin(C_bn(2, 0));               // 俯仰角
//    yaw = atan2(C_bn(1, 0), C_bn(0, 0));    // 航向角
//
//    // 转换为角度
//    roll *= RAD_TO_DEG;
//    pitch *= RAD_TO_DEG;
//    yaw *= RAD_TO_DEG;
//
//    // 确保航向角在0-360度范围内
//    if (yaw < 0) yaw += 360.0;
//}

//IMUDATA parseIMULine(const string& line) {
//    IMUDATA data;
//    // 1. 找到分号位置，提取数据部分
//    size_t semicolon_pos = line.find(';');
//    if (semicolon_pos == string::npos) {
//        return data; // 无效行
//    }
//
//    // 2. 找到星号位置（校验和开始）
//    size_t star_pos = line.find('*', semicolon_pos);
//    if (star_pos == string::npos) {
//        return data; // 无效行
//    }
//
//    // 3. 提取分号后到星号前的数据部分
//    string data_part = line.substr(semicolon_pos + 1, star_pos - semicolon_pos - 1);
//
//    // 4. 按逗号分割数据部分
//    stringstream ss(data_part);
//    string token;
//    vector<string> tokens;
//
//    while (getline(ss, token, ',')) {
//        tokens.push_back(token);
//    }
//
//    // 5. 检查是否有足够字段
//    //数据直接进行转换
//    //转换时，先乘以acc_scale。因为输出的为速度增量，
//    // 还需乘以CPT的采样率值100（如果是其它设备，乘以对应的采样率值即可），
//    // 得到比力值（m/s/s)。另外，y轴给出的负方向，需反号		
//
//    if (tokens.size() >= 11) {
//        data.week = stoi(tokens[2]);
//        data.time = stod(tokens[3]);
//        data.status_hex = tokens[4];
//        data.accel[2] = stoi(tokens[5]) * acc_scale * SampleRate;
//        data.accel[1] = -stoi(tokens[6]) * acc_scale * SampleRate;
//        data.accel[0] = stoi(tokens[7]) * acc_scale * SampleRate;
//        data.gyro[2] = stoi(tokens[8]) * gyo_scale * SampleRate;
//        data.gyro[1] = -stoi(tokens[9]) * gyo_scale * SampleRate;
//        data.gyro[0] = stoi(tokens[10]) * gyo_scale * SampleRate;
//    }
//
//    return data;
//}

//int main() {
//    // 读取IMU数据文件
//    ifstream file("group4.ASC");
//    if (!file.is_open()) {
//        cout << "无法打开文件" << endl;
//        return -1;
//    }
//
//    // 设置当地纬度
//    double latitude = 30.531651244;
//    Matrix3d C_bn_rt = Matrix3d::Zero();
//    double roll, pitch, yaw;
//    vector<IMUDATA> imu_data;
//    string line;
//
//    // 文件输出：每历元姿态角
//    string epoch_filename = "StaticAimResult_Epoch.txt";
//    ofstream epoch_output(epoch_filename);
//    if (!epoch_output) {
//        cerr << "无法创建文件: " << epoch_filename << endl;
//        return 0;
//    }
//    epoch_output << fixed << setprecision(8);
//
//    // 文件输出：每秒平均姿态角
//    string per_second_filename = "StaticAimResult_PerSecond.txt";
//    ofstream per_second_output(per_second_filename);
//    if (!per_second_output) {
//        cerr << "无法创建文件: " << per_second_filename << endl;
//        return 0;
//    }
//    per_second_output << fixed << setprecision(8);
//
//
//    // 按秒分组存储数据
//    map<int, vector<IMUDATA>> second_data_map;
//
//    while (getline(file, line)) {
//        if (line.find("%RAWIMUSXA") != string::npos) {
//            IMUDATA data = parseIMULine(line);
//            imu_data.push_back(data);
//
//            // 计算当前历元姿态角并输出
//            staticCoarseAlignmentNED(data.accel, data.gyro, latitude, roll, pitch, yaw);
//            epoch_output << data.time << "\t" << roll << "\t" << pitch << "\t" << yaw << endl;
//
//            // 按秒分组数据
//            int second = static_cast<int>(data.time);
//            second_data_map[second].push_back(data);
//        }
//    }
//    epoch_output.close();
//    file.close();
//
//    cout << "每历元计算姿态角完成，结果保存在 " << epoch_filename << endl;
//    cout << "成功读取 " << imu_data.size() << " 条IMU数据" << endl;
//
//    if (imu_data.empty()) {
//        cout << "没有有效数据" << endl;
//        return -1;
//    }
//
//    // 计算每秒平均姿态角
//    cout << "\n计算每秒平均姿态角..." << endl;
//    for (auto it = second_data_map.begin(); it != second_data_map.end(); ++it) {
//        int second = it->first;
//        vector<IMUDATA>& data_list = it->second;
//
//        if (data_list.size() < 100) {
//            cout << "警告: 第 " << second << " 秒只有 " << data_list.size() << " 个数据点" << endl;
//        }
//
//        // 计算该秒内加速度和角速度的平均值
//        double avg_acc[3] = { 0.0 }, avg_gyro[3] = { 0.0 };
//        for (const auto& data : data_list) {
//            for (int i = 0; i < 3; i++) {
//                avg_acc[i] += data.accel[i];
//                avg_gyro[i] += data.gyro[i];
//            }
//        }
//
//        for (int i = 0; i < 3; i++) {
//            avg_acc[i] /= data_list.size();
//            avg_gyro[i] /= data_list.size();
//        }
//
//        // 使用平均值计算姿态角
//        staticCoarseAlignmentNED(avg_acc, avg_gyro, latitude, roll, pitch, yaw);
//
//        // 输出每秒平均姿态角
//        per_second_output << second << "\t" << roll << "\t" << pitch << "\t" << yaw << endl;
//    }
//    per_second_output.close();
//    cout << "每秒平均姿态角计算完成，结果保存在 " << per_second_filename << endl;
//
//    // 计算整个数据集的平均值
//    double acc[3] = { 0.0 }, gyro[3] = { 0.0 };
//    CalAverage(imu_data, acc, gyro);
//    Vector3d avg_acc = Vector3d::Zero();
//    Vector3d avg_gyro = Vector3d::Zero();
//    for (int i = 0; i < 3; i++) {
//        avg_acc(i) = acc[i];
//        avg_gyro(i) = gyro[i];
//    }
//
//    cout << "\n加速度计平均值: " << avg_acc.transpose() << " m/s^2" << endl;
//    cout << "陀螺仪平均值: " << avg_gyro.transpose() << " rad/s" << endl;
//
//    // 进行静态解析粗对准(整个数据集平均值)
//    staticCoarseAlignmentNED(acc, gyro, latitude, roll, pitch, yaw);
//
//    // 输出结果
//    cout << "\n=== 整个数据集平均姿态角结果 ===" << endl;
//
//    cout << "\n欧拉角:" << endl;
//    cout << "横滚角 (Roll): " << roll << " 度" << endl;
//    cout << "俯仰角 (Pitch): " << pitch << " 度" << endl;
//    cout << "航向角 (Yaw): " << yaw << " 度" << endl;
//
//
//    return 0;
//}