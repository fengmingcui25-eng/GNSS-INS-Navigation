#include"INS.h"
using namespace std;



ElliPara para;


//BLH转XYZ，输入BLH为弧度
void BLHToXYZ(const double BLH[3], double XYZ[3])
{
	double B, L, H;
	B = BLH[0];
	L = BLH[1];
	H = BLH[2];
	double N = para.a / sqrt(1 - pow(para.e * sin(B), 2));
	XYZ[0] = (N + H) * cos(B) * cos(L);
	XYZ[1] = (N + H) * cos(B) * sin(L);
	XYZ[2] = (N * (1 - pow(para.e, 2)) + H) * sin(B);
}

//XYZ转BLH，输出结果为弧度
void XYZToBLH(const ElliPara* para, const double XYZ[3], double BLH[3])
{
	double X, Y, Z;
	X = XYZ[0];
	Y = XYZ[1];
	Z = XYZ[2];

	if (X == 0 && Y == 0 && Z == 0) 
	{
		BLH[0] = 0;  // 纬度无效
		BLH[1] = 0;  // 经度无效
		BLH[2] = -para->a;  // 大地高 = -长半轴
	}
	else
	{
		BLH[1] = atan2(Y, X);
		double deltaZ, N;
		BLH[0] = atan(Z / sqrt(X * X + Y * Y));//B的初值
		//deltaZ = N * pow(para->e, 2) * sin(BLH[0]);
		deltaZ = Z * pow(para->e, 2);//给deltaZ赋初值
		double B;
		//B = BLH[0];
		int m = 1;
		do
		{
			B = atan((Z + deltaZ) / sqrt(X * X + Y * Y));
			N = para->a / sqrt(1 - pow(para->e * sin(B), 2));
			deltaZ = N * pow(para->e, 2) * sin(B);
			BLH[0] = B;
			m++;
		} while ((fabs(BLH[0] - B) > 1e-12) || m < 7);

		BLH[2] = sqrt(X * X + Y * Y + pow((Z + deltaZ), 2)) - N;
	}
}

//NEU转换矩阵，输入BLH为弧度
void BLHToNEUMat(const double Blh[], Eigen::MatrixXd& Mat)//BLHToNEUMat(Blh, R);
{
	//double B = Blh[0] * pi / 180.0;
	//double L = Blh[1] * pi / 180.0;
	double B = Blh[0];
	double L = Blh[1];
	Mat << -sin(L), cos(L), 0,
		-sin(B) * cos(L), -sin(B) * sin(L), cos(B),
		cos(B)* cos(L), cos(B)* sin(L), sin(B);
}

//定位误差计算
void CompEnudPos(const double Blh[], const double X0[], const double Xr[], Eigen::MatrixXd& dNeu)
{
	Eigen::MatrixXd delta_X(3,1);
	delta_X << Xr[0] - X0[0],
		Xr[1] - X0[1],
		Xr[2] - X0[2];
	Eigen::MatrixXd R(3, 3);
	BLHToNEUMat(Blh, R);
	dNeu = R * delta_X;
}

//卫星高度角方位角计算，输出结果是弧度
void CompSatElAz(const double Blh[], const double Xr[], const double Xs[], double* Elev, double* Azim)//CompSatElAz(Blh, Xr, Xs, &Elev, &Azim);
{
	Eigen::MatrixXd dNeu(3, 1);
	CompEnudPos(Blh, Xr, Xs, dNeu);
	*Elev = atan(dNeu(2) / sqrt(pow(dNeu(0), 2) + pow(dNeu(1), 2)));
	*Azim = atan2(dNeu(1), dNeu(0));
}
