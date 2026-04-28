#pragma once
#include <Eigen\Dense>
#include <Eigen/QR>
#include<cmath>
#include"sockets.h"
#include"const_def.h"

using namespace Eigen;

enum GNSSSys { GPS = 1, BDS, UNKS };
/// <summary>
/// //////////////////////时间转换
/// </summary>
struct COMMONTIME
{
	short Year;
	unsigned short Month;
	unsigned short Day;
	unsigned short Hour;
	unsigned short Minute;
	double Second;
};

struct MJDTIME
{
	int Days;
	double FracDay;
	MJDTIME()
	{
		Days = 0;
		FracDay = 0.0;
	}
};

struct GPSTIME
{
	unsigned short Week;
	double SecOfWeek;
	GPSTIME()
	{
		Week = 0;
		SecOfWeek = 0.0;
	}
};


/// <summary>
/// /////////////坐标转换
/// </summary>
union XYZ
{
	struct
	{
		double x;
		double y;
		double z;
	};
	double xyz[3];
};

union GEOCOOR
{
	struct
	{
		double longitude;
		double latitude;
		double height;
	};
	double BLH[3];//0-B;1-L;2-H
};

union NEU
{
	struct
	{
		double dN;//单位m
		double dE;
		double dU;
	};
	double Neu[3];//0-N;1-E;2-U
};

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



////////////////////////文件读取解码
struct SATOBSDATA
{
	short Prn;
	GNSSSys System;
	double P[2], L[2], D[2];//0->L1;1->L2;
	double cn0[2], LockTime[2];
	int/*unsigned char*/ half[2];
	bool Valid;
	SATOBSDATA()
	{
		Prn = 0;
		System = GPS;
		for (int i = 0; i < 2; i++)
			P[i] = L[i] = D[i] = LockTime[i] = 0.0;
		Valid = false;
	}
};
struct MWGF
{
	short Prn;//卫星号
	GNSSSys Sys;
	double MW;
	double GF;
	double PIF;
	int n; //平滑计数
	double LockTime[2];
	MWGF()
	{
		Prn = 0;
		n = 1;
		Sys = UNKS;
		MW = GF = PIF = 0.0;
		for (int i = 0; i < 2; i++)
			 LockTime[i] = 0.0;
	}
};

struct SATMIDRES
{
	double SatPos[3], SatVel[3];
	double SatClkOft, SatClkSft;
	double Elevation, Azimuth;
	double TropCorr;
	double Tgd1, Tgd2;
	bool Valid;   //false=没有星历或过期,true=计算成功
	SATMIDRES()
	{
		SatPos[0] = SatPos[1] = SatPos[2] = 0.0;
		SatVel[0] = SatVel[1] = SatVel[2] = 0.0;
		Elevation = pi / 2.0;
		Azimuth = 0.0;
		SatClkOft = SatClkSft = 0.0;
		Tgd1 = Tgd2 = TropCorr = 0.0;
		Valid = false;
	}
};

struct EPOCHOBSDATA
{
	GPSTIME Time;
	short SatNum;
	SATOBSDATA SatObs[MAXCHANNUM];
	MWGF ComObs[MAXCHANNUM];
	SATMIDRES SATPVT[MAXCHANNUM];
	double     Pos[3], BLH[3];      // 保存基站或NovAtel接收机定位结果
	EPOCHOBSDATA()
	{
		SatNum = 0;
		for (int i = 0; i < 3; i++)Pos[i] = BLH[i] = 0.0;
	}
};

struct GPSEPHREC
{
	short PRN;
	GNSSSys System;
	GPSTIME TOC, TOE;
	double   ClkBias, ClkDrift, ClkDriftRate;
	double   IODE, IODC;
	double   SqrtA, M0, e, OMEGA/*/升交点经度*/, i0, omega/*/近地幅角*/;
	double   Crs, Cuc, Cus, Cic, Cis, Crc;
	double   DetlaN, OMEGADot, iDot;
	int SVHealth;
	double   TGD1, TGD2;//群延差改正
};

/* 每个历元的定位结果结构体定义 */
struct POSRES
{
	GPSTIME Time;
	double Pos[3], Vel[3];
	double ClkGPS, ClkBDS;
	double PDOP, SigmaPos, SigmaVel;
	int SatNum_tracked;
	int SatNum_used;
	int BDSNum, GPSNum;
};



/// <summary>
/// ///////////////卫星位置计算
/// </summary>




/// <summary>
/// ////////////对流层延迟改正、粗差探测
/// </summary>
const int H0 = 0;//m
const double T0 = 15.0 + 273.16;//K
const double p0 = 1013.25;//mbar
const double RH0 = 0.5;




/*  每颗卫星的单差观测数据定义  */
struct SDSATOBS
{
	short    Prn;
	GNSSSys  System;
	short    Valid;
	double   dP[2], dL[2];   // m
	short    nBas, nRov;   // 存储单差观测值对应的基准和流动站的数值索引号

	SDSATOBS()
	{
		Prn = nBas = nRov = 0;
		System = UNKS;
		dP[0] = dL[0] = dP[1] = dL[1] = 0.0;
		Valid = -1;
	}
};

/*  每个历元的单差观测数据定义  */
struct SDEPOCHOBS
{
	GPSTIME    Time;
	short      SatNum, SatNumGPS, SatNumBDS;
	SDSATOBS   SdSatObs[MAXCHANNUM];
	MWGF       SdCObs[MAXCHANNUM];

	SDEPOCHOBS()
	{
		SatNum = SatNumGPS = SatNumBDS = 0;
	}
};

/*  双差相关的数据定义  */
struct DDCOBS
{
	int RefPrn[2], RefPos[2];         // 参考星卫星号与存储位置，0=GPS; 1=BDS;
	bool RefValid;                    //检测参考星是否变化;true为变化；false为未变化
	int Sats, DDSatNum[2];            // 待估的双差模糊度数量，0=GPS; 1=BDS
	double FixedAmb[MAXCHANNUM * 4];  // 包括双频最优解[0,AmbNum]和次优解[AmbNum,2*AmbNum]
	double ResAmb[2], Ratio;          // LAMBDA浮点解中的模糊度残差
	float  FixRMS[2];                 // 固定解定位中rms误差
	double dPos[3];                   // 基线向量
	double Position[3];               // 解算坐标
	bool bFixed;                      // true为固定，false为未固定

	DDCOBS()
	{
		int i;
		for (i = 0; i < 2; i++) {
			DDSatNum[i] = 0;    // 各卫星系统的双差数量
			RefPos[i] = RefPrn[i] = -1;
		}
		Sats = 0;              // 双差卫星总数
		Position[0] = Position[1] = Position[2] = 0.0;
		dPos[0] = dPos[1] = dPos[2] = 0.0;
		ResAmb[0] = ResAmb[1] = FixRMS[0] = FixRMS[1] = Ratio = 0.0;
		RefValid = false;
		bFixed = false;
		for (i = 0; i < MAXCHANNUM * 2; i++)
		{
			FixedAmb[2 * i + 0] = FixedAmb[2 * i + 1] = 0.0;
		}
	}
};

/* 每个历元单点定位和测速的结果及其精度指标 */
struct PPRESULT
{
	GPSTIME Time;
	double Position[3];
	double Velocity[3];
	double RcvClkOft[2];               /* 0 为GPS钟差; 1=BDS钟差 */
	double RcvClkSft;
	double PDOP, SigmaPos, SigmaVel;  // 精度指标
	short  GPSSatNum, BDSSatNum;      /* 单点定位使用的GPS卫星数 */
	short  AllSatNum;                /* 观测历元的所有卫星数   */
	bool   IsSuccess;                /* 单点定位是否成功, 1为成功, 0为失败 */

	PPRESULT()
	{
		for (int i = 0; i < 3; i++)		Position[i] = Velocity[i] = 0.0;
		RcvClkOft[0] = RcvClkOft[1] = RcvClkSft = 0.0;
		PDOP = SigmaPos = SigmaVel = 999.9;
		GPSSatNum = BDSSatNum = AllSatNum = 0;
		IsSuccess = false;
	}
};

struct FloatRes
{
	double pos[3];
	double Amb[MAXCHANNUM * 2];
	double Qxx[9];
	double Qnn[6400];//(2*MAXCHANNUM)^2
	double Qxn[3*MAXCHANNUM*2];
	int isFix;
	double Ratio;
};

/*  RTK定位的数据定义  */
struct RAWDAT {
	EPOCHOBSDATA BasEpk;
	EPOCHOBSDATA RovEpk;
	SDEPOCHOBS SdObs;
	DDCOBS DDObs;
	GPSEPHREC GpsEph[MAXGPSNUM], BdsEph[MAXBDSNUM];
	POSRES bestPos_rover, bestPos_base;
	FloatRes Res;
};

struct RTKEKF
{
	GPSTIME Time;
	double X[3 + MAXCHANNUM * 2];
	double P [(3 + MAXCHANNUM * 2)* (3 + MAXCHANNUM * 2)];//位置，权阵
	int Index[MAXCHANNUM], nSats, nPos[MAXCHANNUM];
	int FixAmb[MAXCHANNUM];          // 时间更新后上个历元已经固定并传递的模糊度， 1=已固定，-1=未固定或有周跳
	DDCOBS DDObs, CurDDObs;           // 上一个历元和当前历元的双差观测值信息
	SDEPOCHOBS SDObs;                 // 上一个历元的单差观测值
	double X0[3 + MAXCHANNUM * 2], P0[(3 + MAXCHANNUM * 2) * (3 + MAXCHANNUM * 2)];  // 状态备份
	bool IsInit;                      // 滤波是否初始化

	RTKEKF() {
		IsInit = false;
		nSats = 0;
		for (int i = 0; i < MAXCHANNUM; i++) nPos[i] = Index[i] = FixAmb[i] = -1;
		for (int i = 0; i < 3 + MAXCHANNUM * 2; i++) {
			X[i] = X0[i] = 0.0;
			for (int j = 0; j < 3 + MAXCHANNUM * 2; j++) P[i * (3 + MAXCHANNUM * 2) + j] = P0[i * (3 + MAXCHANNUM * 2) + j] = 0.0;
		}
	}
	
};


//时间转换
void CommonTimeToMjdTime(const COMMONTIME* ct, MJDTIME* mjd);

void MjdTimeToCommonTime(const MJDTIME* mjd, COMMONTIME* ct);

void MjdTimeToGPSTime(const MJDTIME* mjd, GPSTIME* Gpst);

void GPSTimeToMjdTime(const GPSTIME* Gpst, MJDTIME* mjd);

void CommonTimeToGPSTime(const COMMONTIME* ct, GPSTIME* Gpst);

void GPSTimeToCommonTime(const GPSTIME* Gpst, COMMONTIME* ct);
//坐标转换
void BLHToXYZ(const ElliPara* para, const double BLH[3], double XYZ[3]);

void XYZToBLH(const ElliPara* para, const double XYZ[3], double BLH[3]);

void BLHToNEUMat(const double Blh[], Eigen::MatrixXd& Mat);

void CompSatElAz(const double Blh[], const double Xr[], const double Xs[], double* Elev, double* Azim);

void CompEnudPos(const double Blh[], const double X0[], const double Xr[], Eigen::MatrixXd& dNeu);
//文件解码
double R8(unsigned char* p);
float R4(unsigned char* p);
int I4(unsigned char* p);
unsigned int UI4(unsigned char* p);
short I2(unsigned char* p);
unsigned short UI2(unsigned char* p);
void DecodeRange(unsigned char* data, int len, EPOCHOBSDATA* obs);
void DecodeGpsEphem(unsigned char* data, int len, GPSEPHREC geph[]);
void DecodeBdsEphem(unsigned char* data, int len, GPSEPHREC beph[]);
int Decode_psrpos(unsigned char* data, int len, POSRES* pos);
int DecodeNovOem7Dat(unsigned char buff[], int& len, EPOCHOBSDATA* obs, GPSEPHREC geph[], GPSEPHREC beph[], POSRES* pos);
//
//钟差计算
bool CompSatClkOff(const int Prn, const GNSSSys Sys, const GPSTIME* t, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, SATMIDRES* Mid);
//GPS卫星位置计算
int CompGPSSatPVT(const int Prn, GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid);
//BDS卫星位置计算
int CompBDSSatPVT(const int Prn, GPSTIME* t, const GPSEPHREC* Eph, SATMIDRES* Mid);
//开普勒计算
double solveKeplerEquation(double M_t, double e);
//地球自转改正
int CorrectEarthRotation(double Xr[], SATMIDRES* Mid);

//对流层延迟，高度角输入单位为角度
double Hopfield(const double H, const double Elev);
//粗差探测
void DetectOutlier(EPOCHOBSDATA* Obs);
void DetectCycleSlip(SDEPOCHOBS* Obs);//void DetectCycleSlip(SDEPOCHOBS* Obs, const EPOCHOBSDATA* BasEpk, const EPOCHOBSDATA* RovEpk);
///////////////单点定位、单点测速

//信号发射时刻卫星位置计算
void NormalizeGPSTime(GPSTIME* t);
void ComputeSatPVTAtSignalTrans(const EPOCHOBSDATA* Epk, GPSEPHREC* Eph, GPSEPHREC* BDSEph, double UserPos[], SATMIDRES Mid[]);

bool SPP(EPOCHOBSDATA* Epoch, RAWDAT* Raw, PPRESULT* Result);
void SPV(EPOCHOBSDATA* Epoch, PPRESULT* Res);

int TimeSyn(FILE* FObs_base, FILE* FObs_rover, RAWDAT* rawdata);
int TimeSynTCP(SOCKET sock_base, SOCKET sock_rover, RAWDAT* rawdata);
void SingleDiff(EPOCHOBSDATA* Obs_rover, EPOCHOBSDATA* Obs_base, SDEPOCHOBS* SDRes);
void DetRefSat(const EPOCHOBSDATA* epkA, const EPOCHOBSDATA* epkB, SDEPOCHOBS* SDObs, DDCOBS* DDObs);
bool RTKFloat(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov);
bool RTKFloatKalman(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov, RTKEKF* rtkekf);

//LAMBDA
void CopyArray(int size, double* dest, double* src);
void MatrixMultiply(int m1, int n1, int m2, int n2, const double M1[], const double M2[], double M3[]);
int MatrixInv(int n, double a[], double b[]);
int LD(int n, const double* Q, double* L, double* D);
void gauss(int n, double* L, double* Z, int i, int j);
void perm(int n, double* L, double* D, int j, double del, double* Z);
void reduction(int n, double* L, double* D, double* Z);
int search(int n, int m, const double* L, const double* D, const double* zs, double* zn, double* s);
int lambda(int n, int m, const double* a, const double* Q, double* Fix, double* s);