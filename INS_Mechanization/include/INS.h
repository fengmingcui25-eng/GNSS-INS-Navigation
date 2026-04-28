#pragma once

#include<iostream>
#include<string>
#include"const_def.h"
#include<Eigen\Dense>
#include <vector>
#include <map>
using namespace std;
using namespace Eigen;

//坐标转换
struct ElliPara//椭球参数
{
	double a;
	double b;
	double e;
	ElliPara()
	{
		a = 6378137.000;
		b = 6356752.314;
		e = 0.081819790992;
	}
};
void BLHToXYZ(const double BLH[3], double XYZ[3]);

void XYZToBLH(const ElliPara* para, const double XYZ[3], double BLH[3]);

void BLHToNEUMat(const double Blh[], Eigen::MatrixXd& Mat);

void CompSatElAz(const double Blh[], const double Xr[], const double Xs[], double* Elev, double* Azim);

void CompEnudPos(const double Blh[], const double X0[], const double Xr[], Eigen::MatrixXd& dNeu);


//解码示例数据
/*加速度计和陀螺数据*/
struct IMUDATA
{
	int week;
	double time;
	string status_hex;
	double gyro[3];//rad
	double accel[3];
	IMUDATA()
	{
		week = 0;
		time = 0.0;
		for (int i = 0; i < 3; i++)gyro[i] = accel[i] = 0.0;
	}
};
/*参考结果数据*/
struct RefRes
{
	double time;
	double PosXYZ[3];
	double PosBLH[3];
	double vel[3];
	double attitude[3];//单位：°
	/*double blh[3];*/
	RefRes()
	{
		time = 0.0;
		for (int i = 0; i < 3; i++)
			PosXYZ[i] = PosBLH[i] = vel[i] = attitude[i] = 0.0;
	}
};
/*标定参数*/
struct BiasData {
	double accel[18];
	double gyro[18];
	double S_a[3], S_g[3];
	double N_a[6];//xy,yx,xz,zx,yz,zy
	double b_a[3], b_g[3];
	double Mat_a[9];//误差补偿矩阵
	double accel_C[3];//改正后输出
	double gyro_C[3];//改正后输出
	BiasData()
	{
		for (int i = 0; i < 3; i++)
		{
			S_a[i] = 0.0; N_a[i] = 0.0; b_a[i] = 0.0; accel_C[i] = 0.0;
			S_g[i] = 0.0; b_g[i] = 0.0; gyro[i] = 0.0;
		}
		for (int i = 0; i < 9; i++)Mat_a[i] = 0.0;
		for (int i = 0; i < 18; i++)accel[i] = 0.0;
	}
};

double R8(unsigned char* p);
int Decode_IMUDATA(unsigned char* data, IMUDATA* obs);
int Decode_RefRes(unsigned char* data, RefRes* obs);
IMUDATA parseIMULine(const string& line);

void CalAverage(const string& filename, double AverData[]);
void Calibration(BiasData* data);
void CorrectError(IMUDATA* imu, BiasData* data);
void InitBias(BiasData* data);

Vector3d vectorNormalize(const Vector3d& vec);
void staticCoarseAlignmentNED(double acc[], double gyro[], double latitude_deg, double& roll, double& pitch, double& yaw);
void staticCoarseAlignmentENU(double acc[], double gyro[], double latitude_deg, double& roll, double& pitch, double& yaw);
//惯导状态
struct InsState
{
	Quaterniond q; // 姿态四元数
	Matrix3d C;    // 方向余弦矩阵
	Vector3d att;  // 姿态
	Vector3d vel;  //速度
	Vector3d pos;  //位置（BLH）,单位：rad，rad，m
	double static_time = 0.0;  //静止状态累计时间
};

Matrix3d Q2DCM(const Eigen::Quaterniond& q);
Quaterniond Euler2Q(const Eigen::Vector3d& a_);
Vector3d DCM2Euler(const Eigen::Matrix3d& C);
Matrix3d skewSymmetric(const Vector3d& v);
double GetGravity(double lat, double h);
static Vector3d GetAgc(const Vector3d& vel, const Vector3d& pos);
void AttitudeUpdate(InsState* state, double gyro_last[], double gyro[]);
void VelUpdate(InsState* state, double gyro_last[], double gyro[], double acc_last[], double acc[]);
void PosUpdate(InsState* state, double vel_last[]);
void ZuptCheck(double acc[], double gyro[], InsState* state);