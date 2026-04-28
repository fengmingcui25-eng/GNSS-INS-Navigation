#include"RTK.h"
void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if (rowToRemove < numRows) {
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
			matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
	}

	matrix.conservativeResize(numRows, numCols);
}
//删除矩阵某一列
void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols() - 1;

	if (colToRemove < numCols) {
		matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
			matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
	}

	matrix.conservativeResize(numRows, numCols);
}
void NormalizeGPSTime(GPSTIME* t)
{
	while (t->SecOfWeek < 0) {
		t->SecOfWeek += 604800.0; // 一周秒数
		t->Week--;
	}
	while (t->SecOfWeek >= 604800.0) {
		t->SecOfWeek -= 604800.0;
		t->Week++;
	}
}

void ComputeSatPVTAtSignalTrans(const EPOCHOBSDATA* Epk, GPSEPHREC* Eph, GPSEPHREC* BDSEph, double UserPos[], SATMIDRES Mid[])
{
	GPSTIME t_trGPS;
	//设置初值𝛿𝑡𝑗 = 0s
	double delta_tj = 0.0, v_delta_tj = 0.0;
	GPSEPHREC* eph = { 0 };
	for (int i = 0; i < Epk->SatNum; i++)
	{
		//计算𝑡𝑡𝑟[𝐺𝑃𝑆]
		t_trGPS.SecOfWeek = Epk->Time.SecOfWeek - Epk->ComObs[i].PIF / C_Light;//delta_tj=0.0
		t_trGPS.Week = Epk->Time.Week;
		NormalizeGPSTime(&t_trGPS);
		//计算卫星钟差𝛿𝑡
		//星历匹配（系统、prn号）
		if (Epk->SatObs[i].System == GPS)
		{
			for (int j = 0; j < MAXGPSNUM; j++)
			{
				eph = Eph + j;
				// eph = BDSEph + j;
				if (Epk->SatObs[i].Prn == eph->PRN)
					break;
				else
					continue;
			}
		}
		else if (Epk->SatObs[i].System == BDS)
		{
			for (int j = 0; j < MAXBDSNUM; j++)
			{
				eph = BDSEph + j;
				if (Epk->SatObs[i].Prn == eph->PRN)
					break;
				else
					continue;
			}
		}
		if (!CompSatClkOff(Epk->SatObs[i].Prn, Epk->SatObs[i].System, &t_trGPS, eph, eph, Mid + i))
		{
			Mid[i].Valid = false;
			continue;
		}
		//重新计算𝑡𝑡𝑟[𝐺𝑃𝑆]
		t_trGPS.SecOfWeek = Epk->Time.SecOfWeek - Epk->ComObs[i].PIF / C_Light - Mid[i].SatClkOft;
		t_trGPS.Week = Epk->Time.Week;
		NormalizeGPSTime(&t_trGPS);
		//计算卫星位置、速度、钟差和钟速（含相对论效应）
		delta_tj = Mid[i].SatClkOft;
		v_delta_tj = Mid[i].SatClkSft;

		if (Epk->SatObs[i].System == GPS)
		{
			CompGPSSatPVT(Epk->SatObs[i].Prn, &t_trGPS, eph, Mid + i);
		}
		else if (Epk->SatObs[i].System == BDS)
		{
			CompBDSSatPVT(Epk->SatObs[i].Prn, &t_trGPS, eph, Mid + i);
		}
		//CompGPSSatPVT函数中只计算了相对论效应改正，
		// CompSatClkOff这个函数计算了不含相对论改正的部分，二者相加是最后的卫星钟差
		Mid[i].SatClkOft = Mid[i].SatClkOft + delta_tj;
		Mid[i].SatClkSft = Mid[i].SatClkSft + v_delta_tj;

		//计算地球旋转角度𝜛𝑒 ∗ Δ𝑡1，计算地球自转改正
		CorrectEarthRotation(UserPos, &Mid[i]);

		//计算对流层延迟
		//double Xr[3] = { UserPos[0], UserPos[1], UserPos[2] };//测站坐标
		//ElliPara para;
		//para.a = 6378137.000;
		//para.b = 6356752.314;
		//para.e = 0.081819790992;
		//double BLH[3];
		//double Xs[3] = { Mid[i]->SatPos[0], Mid[i]->SatPos[1], Mid[i]->SatPos[2] };//卫星坐标
		//XYZToBLH(&para, Xr, BLH);

		//double Elev, Azim;
		//CompSatElAz(BLH, Xr, Xs, &Elev, &Azim);
		//
		//Mid[i]->TropCorr = Hopfield(BLH[2], Elev * 180.0 / pi);

		//计算完成后，将Valid标记为true
		Mid[i].Valid = true;
	}
	//卫星位置与钟差结果输出
}

//bool SPP(EPOCHOBSDATA* Epoch, GPSEPHREC* GPSEph, GPSEPHREC* BDSEph, POSRES* Res)
bool SPP(EPOCHOBSDATA* Epoch, RAWDAT* Raw, PPRESULT* Result)
{
	double deltaX[MAXCHANNUM], deltaY[MAXCHANNUM], deltaZ[MAXCHANNUM], S[MAXCHANNUM], w[MAXCHANNUM];
	//① 设定初始位置𝑋𝑅
	MatrixXd UserPos(5, 1);
	UserPos << -2267335.9727, 5008647.9932, 3222372.6377, 0.0, 0.0;// 0.0, 0.0, 0.0, 0.0, 0.0;//-2267335.9727 Y:5008647.9932 Z:3222372.6377 GPS Clk: -2.561 
	//double Xr[3];
	SATMIDRES Mid[MAXCHANNUM] = { };


	//double Elev, Azim;//高度角、方位角

	double ResPos[3];
	int iter = 0;
	int n;//矩阵维度
	bool converged = false;//判断收敛
	do
	{
		//统计有效卫星
	//② 计算信号发射时刻的卫星位置和钟差 𝑋计算地球自转改正和对流层延迟𝑇𝑗 等
		n = 0;
		int GPSNum = 0, BDSNum = 0;
		ResPos[0] = UserPos(0); ResPos[1] = UserPos(1); ResPos[2] = UserPos(2);
		ComputeSatPVTAtSignalTrans(Epoch,Raw->GpsEph , Raw->BdsEph, ResPos, Mid);

		double Xr[3] = { UserPos(0), UserPos(1), UserPos(2) };//测站坐标
		ElliPara para;
		para.a = 6378137.000;
		para.b = 6356752.314;
		para.e = 0.081819790992;
		double BLH[3];

		for (int i = 0; i < Epoch->SatNum; i++)
		{

			double Xs[3] = { Mid[i].SatPos[0], Mid[i].SatPos[1], Mid[i].SatPos[2] };//卫星坐标
			XYZToBLH(&para, Xr, BLH);
			double Elev, Azim;
			CompSatElAz(BLH, Xr, Xs, &Elev, &Azim);
			Mid[i].Elevation = Elev;
			Mid[i].TropCorr = Hopfield(BLH[2], Elev * 180.0 / pi);
		}

		for (int i = 0; i < Epoch->SatNum; i++)
		{
			if (Mid[i].Valid && Epoch->SatObs[i].Valid)
			{
				n++;
				if (Epoch->SatObs[i].System == GPS) GPSNum++;
				else if (Epoch->SatObs[i].System == BDS) BDSNum++;
			}
		}
		//n即为可用卫星


		MatrixXd B = MatrixXd::Zero(n, 5);
		MatrixXd W = MatrixXd::Zero(n, 1);
		MatrixXd dx;//改正数
		//MatrixXd N = MatrixXd::Zero(5, 5);
		//MatrixXd Y = MatrixXd::Zero(5, 1);

		int validIdx = 0; // 有效卫星计数器，每次迭代需要重置索引
		for (int i = 0; i < Epoch->SatNum; i++)
		{
			//		③ 对所有卫星的观测数据进行线性化
			//		(1) 卫星位置计算失败、观测数据不完整或有粗差，不参与定位计算
			//      已用n处理,但是由于结构体中未删除标记为false的数据，所以仍存在在循环中，用SatNum进行循环，对于错误数据，continue跳过
			if (Mid[i].Valid == false || Epoch->SatObs[i].Valid == false)
			{
				continue;
			}
			//		(2) 以初始位置为参考，对观测方程线性化，计算观测系数矩阵𝐵𝑛∗5和残差向量𝑤 𝑎 + 𝑏 ∗1，统计
		//		参与定位的各系统卫星数和所有卫星数
			S[i] = sqrt(pow(Mid[i].SatPos[0] - UserPos(0), 2) + pow(Mid[i].SatPos[1] - UserPos(1), 2) + pow(Mid[i].SatPos[2] - UserPos(2), 2));
			deltaX[i] = -Mid[i].SatPos[0] + UserPos(0);
			deltaY[i] = -Mid[i].SatPos[1] + UserPos(1);
			deltaZ[i] = -Mid[i].SatPos[2] + UserPos(2);
			if (Epoch->SatObs[i].System == GPS)
			{
				w[i] = Epoch->ComObs[i].PIF - S[i] - UserPos(3) - Mid[i].TropCorr + C_Light * Mid[i].SatClkOft;//数组的实际维度为n
			}
			else if (Epoch->SatObs[i].System == BDS)
			{
				if (GPSNum != 0)
					w[i] = Epoch->ComObs[i].PIF - S[i] - UserPos(4) - Mid[i].TropCorr + C_Light * Mid[i].SatClkOft - C_Light * pow(FG1_BDS, 2) / (pow(FG1_BDS, 2) - pow(FG3_BDS, 2)) * Mid[i].Tgd1;//数组的实际维度为n
				else if (GPSNum == 0)
					w[i] = Epoch->ComObs[i].PIF - S[i] - UserPos(3) - Mid[i].TropCorr + C_Light * Mid[i].SatClkOft - C_Light * pow(FG1_BDS, 2) / (pow(FG1_BDS, 2) - pow(FG3_BDS, 2)) * Mid[i].Tgd1;//数组的实际维度为n

			}

			B(validIdx, 0) = deltaX[i] / S[i];
			B(validIdx, 1) = deltaY[i] / S[i];
			B(validIdx, 2) = deltaZ[i] / S[i];
			if (Epoch->SatObs[i].System == GPS) B(validIdx, 3) = 1;
			else if (Epoch->SatObs[i].System == BDS) B(validIdx, 4) = 1;
			W(validIdx, 0) = w[i];
			validIdx++;
		}
		//		④ 卫星总数是否大于未知参数数量，如果卫星数不足，直接返回定位失败。
		if (n < 5) return false;
		//		⑤ 如果GPS或BDS卫星数量为0，重构法方程矩阵N和Y
		//N = B.transpose() * B;
		//Y = B.transpose() * W;

		//if (GPSNum == 0) {
		//	RemoveColumn(B, 4);
		//	RemoveRow(N, 4);
		//	RemoveColumn(N, 4);
		//	RemoveRow(Y, 4);
		//}
		//if (BDSNum == 0) {
		//	RemoveColumn(B, 5);
		//	RemoveRow(N, 5);
		//	RemoveColumn(N, 5);
		//	RemoveRow(Y, 5);
		//}

		if (GPSNum == 0)
		{
			MatrixXd B_new(n, 4);

			// 位置参数不变
			B_new.col(0) = B.col(0); // dX
			B_new.col(1) = B.col(1); // dY
			B_new.col(2) = B.col(2); // dZ

			// BDS钟差参数移到第3列（原第4列）
			B_new.col(3) = B.col(4); // BDS接收机钟差

			// 替换原矩阵
			B = B_new;
		}
		else if (BDSNum == 0)
		{
			MatrixXd B_new(n, 4);
			B_new.col(0) = B.col(0); // dX
			B_new.col(1) = B.col(1); // dY
			B_new.col(2) = B.col(2); // dZ
			B_new.col(3) = B.col(3); // GPS接收机钟差（保持不变）
			B = B_new;
		}

		//		⑥ 最小二乘求解 𝑥ො = 𝐁
		//		𝑇𝐏𝐁 −1𝑩
		//		𝑻𝐏𝐰
		dx = (B.transpose() * B).inverse() * B.transpose() * W;//N.inverse() * B.transpose()*W;
		//		⑦ 检查𝑥ො是否收敛，如果没有收敛，将初始位置设定为
		// 
		//		⑧ 迭代计算②③④⑤⑥⑦⑧，直到收敛
		UserPos(0) += dx(0);
		UserPos(1) += dx(1);
		UserPos(2) += dx(2);

		double dclkG = 0.0, dclkB = 0.0;

		if (GPSNum > 0 && BDSNum > 0) {
			UserPos(3) += dx(3); dclkG = dx(3);
			UserPos(4) += dx(4); dclkB = dx(4);
		}
		else if (BDSNum == 0) {
			UserPos(3) += dx(3); dclkG = dx(3);
		}
		else if (GPSNum == 0) {
			UserPos(4) += dx(3); dclkB = dx(3);
		}

		double deltaS;//dx平方和开根号
		deltaS = sqrt(pow(dx(0), 2) + pow(dx(1), 2) + pow(dx(2), 2));


		if (deltaS < 1E-3)//(dx(0) < 1E-3&&dx(1)<1E-3&&dx(2)<1E-3&&dx(3)<1E-3)
		{
			
			Result->Position[0] = UserPos(0);
			Result->Position[1] = UserPos(1);
			Result->Position[2] = UserPos(2);
			Result->RcvClkOft[0] = UserPos(3); 
			Result->RcvClkOft[1] = UserPos(4);
			Result->AllSatNum = Epoch->SatNum;
			Result->BDSSatNum = BDSNum;
			Result->GPSSatNum = GPSNum;

			//		⑨ 定位精度评价，计算PDOP和标准差。
			MatrixXd Qxx;
			Qxx = (B.transpose() * B).inverse();
			Result->PDOP = sqrt(Qxx(0, 0) + Qxx(1, 1) + Qxx(2, 2));
			MatrixXd V = MatrixXd::Zero(n, 1);
			V = B * dx - W;
			MatrixXd xigma2(1, 1);
			int t = 0;
			if (GPSNum == 0 || BDSNum == 0) t = 4;
			else t = 5;
			xigma2 = V.transpose() * V / (GPSNum + BDSNum - t);
			Result->SigmaPos = sqrt(xigma2(0));
			Result->Time = Epoch->Time;
			Result->IsSuccess = true;
			converged = true;
			break;
		}

		iter++;
	} while (iter < 15);

	for (int i = 0; i < MAXCHANNUM; i++)
	{
		if (Mid[i].Valid) {
			Epoch->SATPVT[i] = Mid[i];
		}
	}
	
	return converged;
}

void SPV(EPOCHOBSDATA* Epoch, PPRESULT* Res)
{
	MatrixXd RecPos(3, 1);
	RecPos << Res->Position[0], Res->Position[1], Res->Position[2];

	// 统计有效卫星
	int nValid = 0;
	for (int i = 0; i < Epoch->SatNum; i++)
	{
		if (Epoch->SatObs[i].Valid && Epoch->SATPVT[i].Valid)
		{
			nValid++;
		}
	}

	if (nValid < 4)
	{
		return;
	}
	MatrixXd B = MatrixXd::Zero(nValid, 4); // vx, vy, vz, 接收机钟差变化率
	VectorXd W = VectorXd::Zero(nValid);
	int validIdx = 0;
	for (int i = 0; i < Epoch->SatNum; i++) {
		if (!Epoch->SatObs[i].Valid || !Epoch->SATPVT[i].Valid) continue;

		double dx = Res->Position[0] - Epoch->SATPVT[i].SatPos[0];
		double dy = Res->Position[1] - Epoch->SATPVT[i].SatPos[1];
		double dz = Res->Position[2] - Epoch->SATPVT[i].SatPos[2];
		double rho = sqrt(dx * dx + dy * dy + dz * dz);

		double l = dx / rho;
		double m = dy / rho;
		double n = dz / rho;

		double rho_dot = -(dx * Epoch->SATPVT[i].SatVel[0] + dy * Epoch->SATPVT[i].SatVel[1] + dz * Epoch->SATPVT[i].SatVel[2]) / rho;

		double w_j = Epoch->SatObs[i].D[0] - rho_dot + C_Light * Epoch->SATPVT[i].SatClkSft;


		B(validIdx, 0) = l;
		B(validIdx, 1) = m;
		B(validIdx, 2) = n;
		B(validIdx, 3) = 1.0;  // 接收机钟差变化率系数
		W(validIdx) = w_j;
		validIdx++;
	}
	VectorXd X = (B.transpose() * B).inverse() * B.transpose() * W;

	// 保存结果到Res结构体
	Res->Velocity[0] = X(0);  // vx
	Res->Velocity[1] = X(1);  // vy
	Res->Velocity[2] = X(2);  // vz
	Res->RcvClkSft = X(3); // 接收机钟差变化率
	MatrixXd V = MatrixXd::Zero(nValid, 1);
	V = B * X - W;
	MatrixXd xigma2(1, 1);
	int t = 4;
	xigma2 = V.transpose() * V / (Epoch->SatNum - t);
	Res->SigmaVel = sqrt(xigma2(0));

}
