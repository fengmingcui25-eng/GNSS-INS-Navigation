#include "INS.h"
#include <fstream>
//#include <sstream>
//#include <iomanip>
//#include<Eigen\Dense>





void CalAverage(const string& filename,double AverData[])
{
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << endl;
        return ;
    }
    string line;
    int count = 0;

    //初始化accel
    AverData[0] = AverData[1] = AverData[2] = AverData[3] = AverData[4] = AverData[5] = 0.0;
    while (getline(file, line)) {
        // 只处理以 %RAWIMUSA 开头的行
        if (line.find("%RAWIMUSA") != 0) {
            continue; // 跳过其他类型的数据
        }

        //解码
        IMUDATA imu = parseIMULine(line);
        AverData[0] += imu.accel[0];
        AverData[1] += imu.accel[1];
        AverData[2] += imu.accel[2];
        AverData[3] += imu.gyro[0];
        AverData[4] += imu.gyro[1];
        AverData[5] += imu.gyro[2];
        count++;
    }
    AverData[0] = AverData[0] / count;
    AverData[1] = AverData[1] / count;
    AverData[2] = AverData[2] / count;
    AverData[3] = AverData[3] / count;
    AverData[4] = AverData[4] / count;
    AverData[5] = AverData[5] / count;
}

void Calibration(BiasData* data)
{
    MatrixXd L = MatrixXd::Zero(3, 6);
    MatrixXd M = MatrixXd::Zero(3, 4);
    MatrixXd A = MatrixXd::Zero(4, 6);
    MatrixXd Mat_a = MatrixXd::Zero(3, 3);//误差补偿矩阵
    for (int i = 0; i < 6; i++)
    {
        L(0, i) = data->accel[3 * i];
        L(1, i) = data->accel[3 * i + 1];
        L(2, i) = data->accel[3 * i + 2];
        A(3, i) = 1;
    }
    A(0, 0) = A(1, 2) = A(2, 4) = GRAVITY;
    A(0, 1) = A(1, 3) = A(2, 5) = -GRAVITY;
    M = L * A.transpose() * (A * A.transpose()).inverse();

    for (int i = 0; i < 3; i++)
    {
        data->S_a[i] = M(i, i) - 1;
        data->b_a[i] = M(i, 3);
    }
    data->N_a[0] = M(1, 0); data->N_a[1] = M(0, 1);//xy,yx
    data->N_a[2] = M(2, 0); data->N_a[3] = M(0, 2);//xz,zx
    data->N_a[4] = M(2, 1); data->N_a[5] = M(1, 2);//yz,zy
    Mat_a = M.block(0, 0, 3, 3);
    int idx = 0;
    for(int i=0;i<3;i++)
        for (int j = 0; j < 3; j++)
        {
            data->Mat_a[idx] = Mat_a(i, j);//按行存储
            idx++;
        }
    data->b_g[0] = (data->gyro[0] + data->gyro[3]);//x
    data->b_g[1] = (data->gyro[7] + data->gyro[10]);//y
    data->b_g[2] = (data->gyro[14] + data->gyro[17]);//z
    data->S_g[0] = (data->gyro[0] - data->gyro[3]) / 2.0 / w_e / sin(Latitude / 180.0 * PI) - 1;
    data->S_g[1] = (data->gyro[7] - data->gyro[10]) / 2.0 / w_e / sin(Latitude / 180.0 * PI) - 1;
    data->S_g[2] = (data->gyro[14] - data->gyro[17]) / 2.0 / w_e / sin(Latitude / 180.0 * PI) - 1;
}
void CorrectError(IMUDATA* imu, BiasData* data)
{
    MatrixXd Mat_a = MatrixXd::Zero(3, 3);
    MatrixXd Mat_g = MatrixXd::Zero(3, 3);
    MatrixXd b_a = MatrixXd::Zero(3, 1);
    MatrixXd b_g = MatrixXd::Zero(3, 1);
    MatrixXd m_a = MatrixXd::Zero(3, 1);
    MatrixXd m_g = MatrixXd::Zero(3, 1);
    MatrixXd M_c_a = MatrixXd::Zero(3, 1);
    MatrixXd M_c_g = MatrixXd::Zero(3, 1);
    int idx = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            Mat_a(i, j) = data->Mat_a[idx];//按行存储
            idx++;
        }
    for (int i = 0; i < 3; i++)b_a(i) = data->b_a[i];
    m_a(0) = imu->accel[0]; m_a(1) = imu->accel[1]; m_a(2) = imu->accel[2];
    M_c_a = Mat_a.inverse() * (m_a - b_a);
    for (int i = 0; i < 3; i++)
    {
        data->accel_C[i] = M_c_a(i);
        imu->accel[i] = M_c_a(i);
    }
    
    for (int i = 0; i < 3; i++)Mat_g(i, i) = 1 + data->S_g[i];
    for (int i = 0; i < 3; i++)b_g(i) = data->b_g[i];
    m_g(0) = imu->gyro[0]; m_g(1) = imu->gyro[1]; m_g(2) = imu->gyro[2];
    M_c_g = Mat_g.inverse() * (m_g - b_g);
    for (int i = 0; i < 3; i++)
    {
        data->gyro_C[i] = M_c_g(i);
        imu->gyro[i] = M_c_g(i);
    }
}

/*星网宇达误差参数*/
void InitBias(BiasData* data)
{
    //加速度计误差参数
    data->b_a[0] = -0.009769454884589;
    data->b_a[1] = 0.0016622833567754153;
    data->b_a[2] = -0.01668140630928208;
    data->Mat_a[0] = 1 + 0.000238033676738;// 1 + sx;
    data->Mat_a[1] = 0.0039518320167125018;// ryx;
    data->Mat_a[2] = -0.007999129102578061;// rzx;
    data->Mat_a[3] = -0.00414191562;// rxy;
    data->Mat_a[4] = 1 + 0.0003581465048088095;// 1 + sy;
    data->Mat_a[5] = 0.000463114933640740;// rzy;
    data->Mat_a[6] = 0.008386382768;// rxz;
    data->Mat_a[7] = -0.000403957440478629;// ryz;
    data->Mat_a[8] = 1 + 0.000428946210247288;// 1 + sz;
    //陀螺误差参数
    data->b_g[0] = -5.891158223e-6;
    data->b_g[1] = 2.51317605437e-6;
    data->b_g[2] = -4.95370664719e-6;
    data->S_g[0] = -0.13552049256438814;
    data->S_g[1] = -0.099413199347388570;
    data->S_g[2] = -0.12411068696051542;
}

//int main() {
//    IMURes data;
//    vector<string> files = { "x_up.ASC", "x_down.ASC", "y_up.ASC",
//                        "y_down.ASC", "z_up.ASC", "z_down.ASC" };//静置五分钟数据文件
//    vector<string> filesR = { "x_tl_1.ASC", "x_tl_2.ASC", "y_tl_1.ASC",
//                        "y_tl_2.ASC", "z_tl_1.ASC", "z_tl_2.ASC" };//绕轴旋转数据文件
//    double AverData[6];
//    int i = 0;
//    for (const auto& filename : files) {
//        
//        CalAverage(filename, AverData);
//       
//        data.accel[3 * i] = AverData[0];
//        data.accel[3 * i + 1] = AverData[1];
//        data.accel[3 * i + 2] = AverData[2];
//        data.gyro[3 * i] = AverData[3];
//        data.gyro[3 * i + 1] = AverData[4];
//        data.gyro[3 * i + 2] = AverData[5];
//
//        i++;
//
//    }
//    //计算得到误差参数
//    Calibration(&data);
//    
//
//    //改正输出
//    for (const auto& filename : filesR)
//    {
//        ifstream file(filename);
//        if (!file.is_open()) {
//            cerr << "无法打开文件: " << filename << endl;
//            return 0;
//        }
//        string line;
//        string WriteFile_Acc = "AccelResult1" + filename + ".txt";
//        string WriteFile_G = "GyroResult1" + filename + ".txt";
//        ofstream outputAcc(WriteFile_Acc);
//        ofstream outputGyro(WriteFile_G);
//        if (!outputAcc) {
//            cerr << "无法创建文件: " << WriteFile_Acc << endl;
//            continue;
//        }
//        if (!outputGyro) {
//            cerr << "无法创建文件: " << WriteFile_G << endl;
//            continue;
//        }
//        outputAcc << fixed << setprecision(8);
//        outputGyro << fixed << setprecision(8);
//        while (getline(file, line)) {
//            // 只处理以 %RAWIMUSA 开头的行
//            if (line.find("%RAWIMUSA") != 0) {
//                continue; // 跳过其他类型的数据
//            }
//
//            //解码
//            IMUData imu = parseIMULine(line);
//            CorrectError(&imu, &data);
//            //原始观测值+改正值
//            outputAcc << imu.accel_x << "\t" << imu.accel_y << "\t" << imu.accel_z << "\t" << data.accel_C[0] << "\t" << data.accel_C[1] << "\t" << data.accel_C[2] << endl;
//            outputGyro << imu.gyro_x << "\t" << imu.gyro_y << "\t" << imu.gyro_z << "\t" << data.gyro_C[0] << "\t" << data.gyro_C[1] << "\t" << data.gyro_C[2] << endl;
//        }
//        cout << filename+"数据处理完成" << endl;
//        outputAcc.close();
//        outputGyro.close();
//    }
//    cout << "结果已保存" << endl;
//
//
//    
//    return 0;
//}