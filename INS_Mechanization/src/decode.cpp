#include<string.h>
#include<cmath>
#include"INS.h"


double R8(unsigned char* p)
{
    double r;
    memcpy(&r, p, 8);
    return r;
}

int Decode_IMUDATA(unsigned char* data, IMUDATA* obs)
{
    double accel[3] = { 0.0 }, gyro[3] = { 0.0 };
    memset(obs, 0, sizeof(IMUDATA));
    obs->time = R8(data);
    obs->gyro[0] = R8(data + 8);
    obs->gyro[1] = R8(data + 16);
    obs->gyro[2] = R8(data + 24);
    obs->accel[0] = R8(data + 32);
    obs->accel[1] = R8(data + 40);
    obs->accel[2] = R8(data + 48);
    return 0;
}

 /*obs->accel[0] = accel[0];
    obs->accel[1] = accel[1];
    obs->accel[2] = accel[2];
    obs->gyro[0] = gyro[0];
    obs->gyro[1] = gyro[1];
    obs->gyro[2] = gyro[2];*/

int Decode_RefRes(unsigned char* data, RefRes* obs)
{
    memset(obs, 0, sizeof(RefRes));
    obs->time = R8(data);
    obs->PosBLH[0] = R8(data + 8);//deg
    obs->PosBLH[1] = R8(data + 16);//deg
    obs->PosBLH[2] = R8(data + 24);
    obs->vel[0] = R8(data + 32);
    obs->vel[1] = R8(data + 40);
    obs->vel[2] = R8(data + 48);
    obs->attitude[0] = R8(data + 56);//deg
    obs->attitude[1] = R8(data + 64);//deg
    obs->attitude[2] = R8(data + 72);//deg
    double BLH[3];
    for (int i = 0; i < 3; i++)BLH[i] = obs->PosBLH[i];
    for (int i = 0; i < 2; i++)BLH[i] = BLH[i] * DEG_TO_RAD;
    BLHToXYZ(BLH, obs->PosXYZ);
    return 0;
}

IMUDATA parseIMULine(const string& line) {
    IMUDATA data;
    double accel[3] = { 0.0 }, gyro[3] = {0.0};
    // 1. 找到分号位置，提取数据部分
    size_t semicolon_pos = line.find(';');
    if (semicolon_pos == string::npos) {
        return data; // 无效行
    }

    // 2. 找到星号位置（校验和开始）
    size_t star_pos = line.find('*', semicolon_pos);
    if (star_pos == string::npos) {
        return data; // 无效行
    }

    // 3. 提取分号后到星号前的数据部分
    string data_part = line.substr(semicolon_pos + 1, star_pos - semicolon_pos - 1);

    // 4. 按逗号分割数据部分
    stringstream ss(data_part);
    string token;
    vector<string> tokens;

    while (getline(ss, token, ',')) {
        tokens.push_back(token);
    }

    // 5. 检查是否有足够字段
    //数据直接进行转换
    //转换时，先乘以acc_scale。因为输出的为速度增量，
    // 还需乘以CPT的采样率值100（如果是其它设备，乘以对应的采样率值即可），
    // 得到比力值（m/s/s)。另外，y轴给出的负方向，需反号		

    if (tokens.size() >= 9) {
        data.week = stoi(tokens[0]);
        data.time = stod(tokens[1]);
        data.status_hex = tokens[2];
        /*(m/s2)&&&&&&&&(rad/s)*/
        /*data.accel[2] = stoi(tokens[3]) * acc_scale * SampleRate;
        data.accel[1] = -stoi(tokens[4]) * acc_scale * SampleRate;
        data.accel[0] = stoi(tokens[5]) * acc_scale * SampleRate;
        data.gyro[2] = stoi(tokens[6]) * gyo_scale * SampleRate;
        data.gyro[1] = -stoi(tokens[7]) * gyo_scale * SampleRate;
        data.gyro[0] = stoi(tokens[8]) * gyo_scale * SampleRate;*/

        /*(m/s)&&&&&&&&(rad)*/
        accel[2] = stoi(tokens[3]) * acc_scale;// *SampleRate;
        accel[1] = -stoi(tokens[4]) * acc_scale;// *SampleRate;
        accel[0] = stoi(tokens[5]) * acc_scale;// *SampleRate;
        gyro[2] = stoi(tokens[6]) * gyo_scale;// *SampleRate;
        gyro[1] = -stoi(tokens[7]) * gyo_scale;// *SampleRate;
        gyro[0] = stoi(tokens[8]) * gyo_scale;// *SampleRate;
    }

    data.accel[0] = accel[1];
    data.accel[1] = accel[0];
    data.accel[2] = -accel[2];
    data.gyro[0] = gyro[1];
    data.gyro[1] = gyro[0];
    data.gyro[2] = -gyro[2];

    return data;
}