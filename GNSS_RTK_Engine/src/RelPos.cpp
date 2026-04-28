#include"RTK.h"
#include<iostream>
using namespace std;


//最小二乘求解
bool RTKFloat(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov)
{
	if (Raw->SdObs.Time.SecOfWeek == 470765)// 问题历元上一个历元   
		int wwwww = 1;

	if (Raw->SdObs.SatNum < 6)return false;//卫星少
	if (Raw->SdObs.SatNumBDS < 2 || Raw->SdObs.SatNumGPS < 2)return false;//双系统

	//1. 设置基站和流动站位置初值
	double X0_base[3], X0_rov[3];//基站和流动站初始位置
	MatrixXd B = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, 3 + (Raw->SdObs.SatNum - 2) * 2);
	MatrixXd W = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, 1);
	MatrixXd X_float = MatrixXd::Zero(3 , 1);
	MatrixXd X_fix = MatrixXd::Zero(3, 1);
	MatrixXd deltaX = MatrixXd::Zero(3, 1);
	MatrixXd X = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 1);
	MatrixXd P = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, (Raw->SdObs.SatNum - 2) * 4);
	MatrixXd dx;//改正数
	MatrixXd Q = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);//协因数阵
	double w[100]; 
	MatrixXd N = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 2, 1);
	MatrixXd QNN = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 2, (Raw->SdObs.SatNum - 2) * 2);
	MatrixXd QXN = MatrixXd::Zero(3, (Raw->SdObs.SatNum - 2) * 2);
	for (int i = 0; i < 3; i++)
		X0_base[i] = Raw->BasEpk.Pos[i];
	if (Rov->IsSuccess) {
		for (int i = 0; i < 3; i++)X0_rov[i] = Rov->Position[i];
	}
	else return false;//流动站单点定位未成功
	for (int i = 0; i < 3; i++)X(i) = X0_base[i];///////
	
	//参考星选取不成功
	if (Raw->DDObs.RefPos[0] == -1 || Raw->DDObs.RefPos[1] == -1|| Raw->DDObs.RefPrn[0] == -1 || Raw->DDObs.RefPrn[1] == -1)return false;

	/*为模糊度参数赋初值，初值不参与后续迭代*/
	int idx = 0;//独立索引，跳过参考星
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		//	① 参考星不用计算，存在半周不参与计算
		if (Raw->SdObs.SdSatObs[i].System == GPS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0])continue;
		if (Raw->SdObs.SdSatObs[i].System == BDS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1])continue;//参考星不参与计算
		if (!Raw->SdObs.SdSatObs[i].Valid)continue;//周跳不参与计算
		if (Raw->SdObs.SdSatObs[i].System == GPS)
		{
			X(2 * idx + 3) = ((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[0]) -
				(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[0])) / WL1_GPS;
			X(2 * idx + 4) = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[1]) -
				(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[1])) / WL2_GPS;
			N(2 * idx) = X(2 * idx + 3); N(2 * idx + 1) = X(2 * idx + 4);
		}
		else if (Raw->SdObs.SdSatObs[i].System == BDS)
		{
			X(2 * idx + 3) = ((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[0]) -
				(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[0])) / WL1_BDS;
			X(2 * idx + 4) = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[1]) -
				(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[1])) / WL3_BDS;
			N(2 * idx) = X(2 * idx + 3); N(2 * idx + 1) = X(2 * idx + 4);
		}
		idx++;
	}

	//	2. 计算GPS和BDS双差卫星数（单双频混用：各系统每个频率的双差卫
	//	星数）
	Raw->DDObs.Sats = Raw->SdObs.SatNum-2;//单差总数-1：单频,单观测值
	Raw->DDObs.DDSatNum[0] = Raw->SdObs.SatNumGPS - 1;//GPS,单频
	Raw->DDObs.DDSatNum[1] = Raw->SdObs.SatNumBDS - 1;//BDS,单频

	
	double S_bas[MAXCHANNUM], S_rov[MAXCHANNUM], DDS[MAXCHANNUM], DDl[MAXCHANNUM], DDm[MAXCHANNUM], DDn[MAXCHANNUM];
	
	//	3. 计算基站坐标到所有卫星的几何距离
	double S_bas_ref[2];//0-->GPS,1-->BDS,到参考星距离
	S_bas_ref[0] = sqrt(pow(X0_base[0] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nBas].SatPos[0], 2) +
		pow(X0_base[1] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nBas].SatPos[1], 2) +
		pow(X0_base[2] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nBas].SatPos[2], 2));
	S_bas_ref[1] = sqrt(pow(X0_base[0] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nBas].SatPos[0], 2) +
		pow(X0_base[1] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nBas].SatPos[1], 2) +
		pow(X0_base[2] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nBas].SatPos[2], 2));
	
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		if (!Raw->SdObs.SdSatObs[i].Valid)continue;
		S_bas[i] = sqrt(pow(X0_base[0] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[i].nBas].SatPos[0], 2) +
			pow(X0_base[1] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[i].nBas].SatPos[1], 2) +
			pow(X0_base[2] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[i].nBas].SatPos[2], 2));
	}

	int iter = 0;
	bool converged = false;//判断收敛
	do
	{//	4. 计算流动站到参考星的几何距离
		double S_rov_ref[2];//0-->GPS,1-->BDS,到参考星距离
		S_rov_ref[0] = 
			sqrt(pow(X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[0], 2) +
			pow(X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[1], 2) +
			pow(X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[2], 2));
		S_rov_ref[1] = 
			sqrt(pow(X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[0], 2) +
			pow(X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[1], 2) +
			pow(X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[2], 2));
		for (int i = 0; i < Raw->SdObs.SatNum; i++)
		{
			if (!Raw->SdObs.SdSatObs[i].Valid)continue;
			S_rov[i] =
				sqrt(pow(X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[0], 2) +
				pow(X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[1], 2) +
				pow(X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[2], 2));
		}
		//	5. 对单差观测值进行循环

		int validIdx = 0; //
		for (int i = 0; i < Raw->SdObs.SatNum; i++)
		{
			//	① 参考星不用计算，存在半周不参与计算
			if (Raw->SdObs.SdSatObs[i].System == GPS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0])continue;
			if (Raw->SdObs.SdSatObs[i].System == BDS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1])continue;//参考星不参与计算
			//if (!Raw->SdObs.SdSatObs[i].Valid)continue;//存在半周不参与计算;有周跳不参与计算
			if (!Raw->SdObs.SdSatObs[i].Valid)continue;
			//	② 线性化双差观测方程得到B矩阵和W向量
			if (Raw->SdObs.SdSatObs[i].System == GPS)
			{
				DDS[i] = S_rov[i] - S_rov_ref[0] - S_bas[i] + S_bas_ref[0];
				DDl[i] = (X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[0]) / S_rov[i] -
					(X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[0]) / S_rov_ref[0];
				DDm[i] = (X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[1]) / S_rov[i] -
					(X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[1]) / S_rov_ref[0];
				DDn[i] = (X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[2]) / S_rov[i] -
					(X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[2]) / S_rov_ref[0];
				W(4 * validIdx) = Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[0] - DDS[i];
				W(4 * validIdx + 1) = Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[1] - DDS[i];
				
				W(4 * validIdx + 2) = Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[0] - DDS[i] - WL1_GPS * X(2 * validIdx + 3);
				W(4 * validIdx + 3) = Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[1] - DDS[i] - WL2_GPS * X(2 * validIdx + 4);
				B(4 * validIdx + 2, 2 * validIdx + 3) = WL1_GPS;//1-->2,3;1-->3,4;2-->6,5;7,6
				B(4 * validIdx + 3, 2 * validIdx + 4) = WL2_GPS;
			}
			else if (Raw->SdObs.SdSatObs[i].System == BDS)
			{
				DDS[i] = S_rov[i] - S_rov_ref[1] - S_bas[i] + S_bas_ref[1];
				DDl[i] = (X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[0]) / S_rov[i] -
					(X0_rov[0] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[0]) / S_rov_ref[1];
				DDm[i] = (X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[1]) / S_rov[i] -
					(X0_rov[1] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[1]) / S_rov_ref[1];
				DDn[i] = (X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[2]) / S_rov[i] -
					(X0_rov[2] - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[2]) / S_rov_ref[1];
				W(4 * validIdx) = Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[0] - DDS[i];
				W(4 * validIdx + 1) = Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[1] - DDS[i];
				
				W(4 * validIdx + 2) = Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[0] - DDS[i] - WL1_BDS * X(2 * validIdx + 3);//
				W(4 * validIdx + 3) = Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[1] - DDS[i] - WL3_BDS * X(2 * validIdx + 4);
				B(4 * validIdx + 2, 2 * validIdx + 3) = WL1_BDS;//1-->2,3;1-->3,4;2-->6,5;7,6
				B(4 * validIdx + 3, 2 * validIdx + 4) = WL3_BDS;
			}
			//B[(Raw->SdObs.SatNum-2)*4,3+(Raw->SdObs.SatNum-2)*2]
			//W[(Raw->SdObs.SatNum-2)*4,1]
			B(4 * validIdx, 0) = B(4 * validIdx + 1, 0) = B(4 * validIdx + 2, 0) = B(4 * validIdx + 3, 0) = DDl[i];
			B(4 * validIdx, 1) = B(4 * validIdx + 1, 1) = B(4 * validIdx + 2, 1) = B(4 * validIdx + 3, 1) = DDm[i];
			B(4 * validIdx, 2) = B(4 * validIdx + 1, 2) = B(4 * validIdx + 2, 2) = B(4 * validIdx + 3, 2) = DDn[i];

			w[4 * validIdx] = W(4 * validIdx); w[4 * validIdx + 1] = W(4 * validIdx + 1);
			w[4 * validIdx + 2] = W(4 * validIdx + 2); w[4 * validIdx + 3] = W(4 * validIdx + 3);
			validIdx++;//validIdx最大Raw->SdObs.SatNum-3
		}
		/*if (iter == 0)
		{
			cout << "W" << endl;
			cout << W << endl << endl;
		}*/

		//	③ 计算权矩阵
		// 设置观测值权重
		double w_p = 2 * 0.3 * 0.3;   // 伪距权重因子
		double w_c = 2 * 0.01 * 0.01; // 载波相位权重因子
		// GPS系统权矩阵块
		if (Raw->DDObs.DDSatNum[0] > 0&& iter<1) {
			// 计算GPS系统的权重参数
			double diag[2], offdiag[2];
			diag[0] = Raw->DDObs.DDSatNum[0] / (w_p * (Raw->DDObs.DDSatNum[0] + 1));      // 伪距主对角线元素
			offdiag[0] = -1.0 / (w_p * (Raw->DDObs.DDSatNum[0] + 1)); // 伪距非对角线元素
			diag[1] = Raw->DDObs.DDSatNum[0] / (w_c * (Raw->DDObs.DDSatNum[0] + 1));       // 载波主对角线元素
			offdiag[1] = -1.0 / (w_c * (Raw->DDObs.DDSatNum[0] + 1));  // 载波非对角线元素

			for (int i = 0; i < Raw->DDObs.DDSatNum[0]; i++) {
				for (int j = 0; j < Raw->DDObs.DDSatNum[0]; j++) {
					// 计算权重值
					double value_p = (i == j) ? diag[0] : offdiag[0];
					double value_l = (i == j) ? diag[1] : offdiag[1];
					// 每个卫星的4个观测值
					int row = 4 * i;
					int col = 4 * j;
					P(row, col) = value_p;
					P(row+1, col+1) = value_p;
					P(row+2, col+2) = value_l;
					P(row+3, col+3) = value_l;
				}
			}
		}

		

		// BDS系统权矩阵块
		if (Raw->DDObs.DDSatNum[1] > 0 && iter < 1) {
			// 计算BDS系统的权重参数
			double diag[2], offdiag[2];
			diag[0] = Raw->DDObs.DDSatNum[1] / (w_p * (Raw->DDObs.DDSatNum[1] + 1));      // 伪距主对角线元素
			offdiag[0] = -1.0 / (w_p * (Raw->DDObs.DDSatNum[1] + 1)); // 伪距非对角线元素
			diag[1] = Raw->DDObs.DDSatNum[1] / (w_c * (Raw->DDObs.DDSatNum[1] + 1));       // 载波主对角线元素
			offdiag[1] = -1.0 / (w_c * (Raw->DDObs.DDSatNum[1] + 1));  // 载波非对角线元素


			int bds_idx = 4 * Raw->DDObs.DDSatNum[0];

			for (int i = 0; i < Raw->DDObs.DDSatNum[1]; i++) {
				for (int j = 0; j < Raw->DDObs.DDSatNum[1]; j++) {
					// 计算权重值
					double value_p = (i == j) ? diag[0] : offdiag[0];
					double value_l = (i == j) ? diag[1] : offdiag[1];

					int row = bds_idx + 4 * i;
					int col = bds_idx + 4 * j;
					P(row, col) = value_p;
					P(row + 1, col + 1) = value_p;
					P(row + 2, col + 2) = value_l;
					P(row + 3, col + 3) = value_l;
				}
			}
		}

		/*if (iter == 0)cout << W << endl;
		cout << endl;*/
		/*if (iter == 0)cout << P << endl;
		cout << endl;
		if (iter == 0)cout << B.transpose() * P * B << endl;
		cout << endl;*/
		//	6. 双差观测方程数量大于未知数，则求解
		if ((Raw->SdObs.SatNum - 2) * 4 < 3 + (Raw->SdObs.SatNum - 2) * 2)break;
		//	7. 最小二乘解算
		dx = (B.transpose() * P * B).inverse() * B.transpose() * P * W;

		/*JacobiSVD<MatrixXd> svd(B.transpose()* P* B, ComputeThinU | ComputeThinV);
		dx = svd.solve(B.transpose() * P * W);*/
		//	8. 更新流动站的位置和双差模糊度参数
		X = X + dx;
		for (int rovidx = 0; rovidx < 3; rovidx++)X0_rov[rovidx] = X(rovidx);
		//for (int rovidx = 0; rovidx < 3; rovidx++)X0(rovidx) = X(rovidx);
		for (int Nidx = 0; Nidx < (Raw->SdObs.SatNum - 2) * 2; Nidx++)N(Nidx) = X(Nidx + 3);
		//	9. 流动站位置增量大于阈值或迭代次数小于阈值，返回4迭代计算，否则浮点解计算完成。
		double deltaS;//dx平方和开根号
		deltaS = sqrt(pow(dx(0), 2) + pow(dx(1), 2) + pow(dx(2), 2));
		double amb_norm = 0;// dx.tail((Raw->SdObs.SatNum - 2) * 2).norm(); // 模糊度参数变化
		double sum_sq = 0.0;
		int amb_count = (Raw->SdObs.SatNum - 2) * 2;
		for (int i = 0; i < amb_count; i++) {
			double d = dx(3 + i);  // 模糊度参数的增量
			sum_sq += d * d;
		}
		amb_norm = sqrt(sum_sq);
		if ((deltaS < 1E-5&& amb_norm<1E-3 )||iter==29)
		{
            //	10. 保存双差浮点解模糊度及其协因数矩阵，用于LAMBDA模糊度固定
			Raw->Res.pos[0] = X(0);
			Raw->Res.pos[1] = X(1);
			Raw->Res.pos[2] = X(2);//浮点解计算位置
			for (int a = 0; a < 3; a++)X_float(a) = X(a);
			for (int a = 0; a < (Raw->SdObs.SatNum - 2) * 2; a++)
			{
				Raw->Res.Amb[a] = X(3 + a);//模糊度，未固定
				N(a)= X(3 + a);
			}
			

			//	11. 精度评价
			//Q = (B.transpose() * P * B).inverse();
			MatrixXd H = B.transpose() * P * B;
			// 4. 创建单位矩阵（与H同维）
			MatrixXd I = MatrixXd::Identity(H.rows(), H.cols());
			//H = H + 1E-6 * I;
			// 2. 使用列主元QR分解（更稳定）
			ColPivHouseholderQR<MatrixXd> qr(H);

			// 3. 检查矩阵是否可逆
			if (qr.rank() < H.cols()) {
				cerr << "警告：矩阵秩亏！秩 = " << qr.rank()
					<< "/" << H.cols() << endl;
			}

			// 5. 求解 Q = H（通过解方程 H*Q = I）
			MatrixXd Q = qr.solve(I);

			/*cout << B.transpose() * P * B;
			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;*/
			int Qxxidx = 0;
			for (int j = 0; j < 3; j++) //按列保存Qxx，（0,0），（1,0）,（2,0），（0,1）......
			{  // 按列循环
				for (int i = 0; i < 3; i++) {  // 按行循环
					Raw->Res.Qxx[Qxxidx++] = Q(i, j);
				}
			}
			int QnnIdx = 0;
			/*for (int i = 0;i< (Raw->SdObs.SatNum - 2) * 2; i++)*///按列保存Qnn
				for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
					for (int i = 0; i < (Raw->SdObs.SatNum - 2) * 2; i++)
				{
					Raw->Res.Qnn[QnnIdx] = Q(j + 3, i + 3);
					QnnIdx++;
				}
			for (int i = 0; i < (Raw->SdObs.SatNum - 2) * 2; i++)//
				for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
				{
					QNN(j, i) = Q(j + 3, i + 3);
				}
			
			int QxnIdx=0;		
			/*for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)*///按列保存Qxn
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
				{
					Raw->Res.Qxn[QxnIdx] = Q(i, j + 3);
					QxnIdx++;
				}
			for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)//按列保存Qxn
				for (int i = 0; i < 3; i++)
				{
					QXN(i, j) = Q(i, j + 3);
				}
				converged = true;
			break;
		}


		iter++;
	}while (iter < 30);
	
	
	MatrixXd FixedN = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 2, 1);
	 Raw->Res.isFix = lambda((Raw->SdObs.SatNum - 2) * 2, 2, Raw->Res.Amb, Raw->Res.Qnn, Raw->DDObs.FixedAmb, Raw->DDObs.ResAmb);
	 if (Raw->Res.isFix != 0) {
		 converged = false;//lambda函数运行失败
		 Raw->DDObs.bFixed = false;
		 return converged;
	 }
	 
	Raw->DDObs.Ratio = Raw->Res.Ratio = Raw->DDObs.ResAmb[1] / Raw->DDObs.ResAmb[0];
	if (Raw->DDObs.Ratio <= 3.0)
	{
		/*for (int a = 0; a < (Raw->SdObs.SatNum - 2) * 2; a++)
			FixedN(a) = Raw->DDObs.FixedAmb[a + (Raw->SdObs.SatNum - 2) * 2];*/
		for (int a = 0; a < 3; a++)Raw->DDObs.Position[a] = X_float(a);
		for (int a = 0; a < 3; a++)Raw->DDObs.dPos[a] = X_float(a) - X0_base[a];
		Raw->DDObs.bFixed = false;
		return converged;
	}
	else {
		for (int a = 0; a < (Raw->SdObs.SatNum - 2) * 2; a++)
			FixedN(a) = Raw->DDObs.FixedAmb[a];
		Raw->DDObs.bFixed = true;
		deltaX = QXN * QNN.inverse() * (FixedN - N);
		X_fix = X_float + QXN * QNN.inverse() * (FixedN - N);
		for (int a = 0; a < 3; a++)Raw->DDObs.dPos[a] = X_fix(a) - X0_base[a];
		for (int a = 0; a < 3; a++)Raw->DDObs.Position[a] = X_fix(a);
	}
	return converged;
}






//滤波初始化
void KfInit(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov, RTKEKF* rtkekf)
{
	//初始化流动站位置
	for (int i = 0; i < 3; i++)rtkekf->X[i] = Rov->Position[i];
	//初始化模糊度
	int idx = 0;//独立索引，跳过参考星
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		//	① 参考星不用计算，存在半周不参与计算
		if (Raw->SdObs.SdSatObs[i].System == GPS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0])continue;
		if (Raw->SdObs.SdSatObs[i].System == BDS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1])continue;//参考星不参与计算
		if (!Raw->SdObs.SdSatObs[i].Valid)continue;//周跳不参与计算
		if (Raw->SdObs.SdSatObs[i].System == GPS)
		{
			rtkekf->X[2 * idx + 3] = ((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[0]) -
				(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[0])) / WL1_GPS;
			rtkekf->X[2 * idx + 4] = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[1]) -
				(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[1])) / WL2_GPS;
		}
		else if (Raw->SdObs.SdSatObs[i].System == BDS)
		{
			rtkekf->X[2 * idx + 3] = ((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[0]) -
				(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[0])) / WL1_BDS;
			rtkekf->X[2 * idx + 4] = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[1]) -
				(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[1])) / WL3_BDS;
		}
		idx++;
	}
	
	int CurIdx = 0; // 当前历元有效卫星计数

	// 为所有非参考星卫星建立索引（包括第一历元）
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		// 跳过参考星
		if (Raw->SdObs.SdSatObs[i].System == GPS &&
			Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0]) continue;
		if (Raw->SdObs.SdSatObs[i].System == BDS &&
			Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1]) continue;

		// 保存当前卫星索引信息
		rtkekf->Index[CurIdx] = i;  // 观测数组索引
		rtkekf->nPos[CurIdx] = 3 + 2 * CurIdx;  // 状态向量位置

		CurIdx++;
	}
	MatrixXd P = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);
	for (int i = 0; i < 3; i++)P(i, i) = 64;//设置X精度为8m，默认第一历元
	for (int i = 3; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)P(i, i) = 25;//设置N精度为5m

	int ppidx = 0;
	for (int i = 0; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)
		for (int j = 0; j < 3 + (Raw->SdObs.SatNum - 2) * 2; j++)//?
		{
			rtkekf->P[ppidx] = P(j, i);
			ppidx++;
		}

	rtkekf->Time = Raw->SdObs.Time;
	rtkekf->nSats = idx;

	rtkekf->DDObs = Raw->DDObs;
	rtkekf->SDObs = Raw->SdObs;
	rtkekf->IsInit = true;

}

void TimePredict(RAWDAT* Raw, RTKEKF* rtkekf)
{
	//**********************1.时间预测**********************

	MatrixXd X = MatrixXd::Zero(3 + (rtkekf->SDObs.SatNum - 2) * 2, 1);/*前一时刻X*/
	MatrixXd X_predict = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 1);/*时间预测X*/
	MatrixXd x_new_slip = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 1);/*//用于设置有周跳和新升起的卫星的时间预测值*/
	MatrixXd phi = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (rtkekf->SDObs.SatNum - 2) * 2);/*状态转移矩阵*/
	MatrixXd P_last = MatrixXd::Zero(3 + (rtkekf->SDObs.SatNum - 2) * 2, 3 + (rtkekf->SDObs.SatNum - 2) * 2);/*//状态协方差P*/
	MatrixXd Q = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);//协因数阵
	MatrixXd G = MatrixXd::Identity(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);//过程噪声驱动矩阵 Γ
	MatrixXd P_predict = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);/*时间预测P*/



	for (int i = 0; i < 3 + (rtkekf->SDObs.SatNum - 2) * 2; i++)X(i) = rtkekf->X[i];
	int ppidx = 0;
	for (int i = 0; i < 3 + (rtkekf->SDObs.SatNum - 2) * 2; i++)
		for (int j = 0; j < 3 + (rtkekf->SDObs.SatNum - 2) * 2; j++)//?
		{
			P_last(j, i) = rtkekf->P[ppidx];
			ppidx++;
		}

	//默认流动站静止
	for (int i = 0; i < 3; i++)phi(i, i) = 1;

	// 保存上一历元的索引（仅在非第一历元需要）
	int prev_Index[MAXCHANNUM] = {0};
	int prev_nPos[MAXCHANNUM] = {0};
	int prev_nSats = rtkekf->nSats; // 上一历元有效卫星数
	if (rtkekf->SDObs.SatNum > 2) { // 非第一历元时保存
		for (int i = 0; i < MAXCHANNUM; i++) {
			prev_Index[i] = rtkekf->Index[i];
			prev_nPos[i] = rtkekf->nPos[i];
		}
	}

	// 初始化索引数组（无论是否第一历元）
	for (int i = 0; i < MAXCHANNUM; i++) {
		rtkekf->Index[i] = -1;
		rtkekf->nPos[i] = -1;
	}

	int CurIdx = 0; // 当前历元有效卫星计数

	// 为所有非参考星卫星建立索引（包括第一历元）
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{

		// 跳过参考星
		if (Raw->SdObs.SdSatObs[i].System == GPS &&
			Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0]) continue;
		if (Raw->SdObs.SdSatObs[i].System == BDS &&
			Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1]) continue;

		// 保存当前卫星索引信息
		rtkekf->Index[CurIdx] = i;  // 观测数组索引
		rtkekf->nPos[CurIdx] = 3 + 2 * CurIdx;  // 状态向量位置

		CurIdx++;
	}



	bool isMatch;
	CurIdx = 0;
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		isMatch = false;
		double w1, w2;
		short ref;
		if (Raw->SdObs.SdSatObs[i].System == GPS) {
			w1 = WL1_GPS; w2 = WL2_GPS; ref = Raw->DDObs.RefPos[0];
		}
		else if (Raw->SdObs.SdSatObs[i].System == BDS) {
			w1 = WL1_BDS; w2 = WL3_BDS; ref = Raw->DDObs.RefPos[1];
		}
		// 跳过参考星
		if (Raw->SdObs.SdSatObs[i].System == GPS &&
			Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0]) continue;
		if (Raw->SdObs.SdSatObs[i].System == BDS &&
			Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1]) continue;
		//有周跳
		if (!Raw->SdObs.SdSatObs[i].Valid)
		{
			//重新进行初始化，phi矩阵对应位置默认0
			x_new_slip(rtkekf->nPos[CurIdx]) =
				((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[ref].dL[0]) -
					(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[ref].dP[0])) / w1;
			x_new_slip(rtkekf->nPos[CurIdx] + 1) = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[ref].dL[1]) -
				(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[ref].dP[1])) / w2;
			Q(rtkekf->nPos[CurIdx], rtkekf->nPos[CurIdx]) = 1E5;
			Q(rtkekf->nPos[CurIdx] + 1, rtkekf->nPos[CurIdx] + 1) = 1E5;
			continue;
		}
		// 使用保存的上一历元索引进行匹配
		for (int j = 0; j < prev_nSats; j++)
		{
			if (prev_Index[j] == -1) continue; // 无效索引跳过

			SDSATOBS& prevSat = rtkekf->SDObs.SdSatObs[prev_Index[j]]; // 正确访问上一历元卫星
			if (Raw->SdObs.SdSatObs[i].System == prevSat.System &&
				Raw->SdObs.SdSatObs[i].Prn == prevSat.Prn) {
				//匹配到了，但是参考星发生了变化
				if (Raw->SdObs.SdSatObs[i].System == GPS && Raw->SdObs.SdSatObs[i].Prn == rtkekf->DDObs.RefPrn[0])
				{
					//重新进行初始化，phi矩阵对应位置默认0
					x_new_slip(rtkekf->nPos[CurIdx]) =
						((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[ref].dL[0]) -
							(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[ref].dP[0])) / w1;
					x_new_slip(rtkekf->nPos[CurIdx] + 1) = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[ref].dL[1]) -
						(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[ref].dP[1])) / w2;
					Q(rtkekf->nPos[CurIdx], rtkekf->nPos[CurIdx]) = 1E5;
					Q(rtkekf->nPos[CurIdx] + 1, rtkekf->nPos[CurIdx] + 1) = 1E5;
					break;
				}
				if (Raw->SdObs.SdSatObs[i].System == BDS && Raw->SdObs.SdSatObs[i].Prn == rtkekf->DDObs.RefPrn[1])
				{
					//重新进行初始化，phi矩阵对应位置默认0
					x_new_slip(rtkekf->nPos[CurIdx]) =
						((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[ref].dL[0]) -
							(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[ref].dP[0])) / w1;
					x_new_slip(rtkekf->nPos[CurIdx] + 1) = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[ref].dL[1]) -
						(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[ref].dP[1])) / w2;
					Q(rtkekf->nPos[CurIdx], rtkekf->nPos[CurIdx]) = 1E5;
					Q(rtkekf->nPos[CurIdx] + 1, rtkekf->nPos[CurIdx] + 1) = 1E5;
					break;
				}
				phi(rtkekf->nPos[CurIdx], prev_nPos[j]) = 1;
				phi(rtkekf->nPos[CurIdx] + 1, prev_nPos[j] + 1) = 1;
				// 设置过程噪声
				Q(rtkekf->nPos[CurIdx], rtkekf->nPos[CurIdx]) = 1E-10;
				Q(rtkekf->nPos[CurIdx] + 1, rtkekf->nPos[CurIdx] + 1) = 1E-10;

				isMatch = true;
				break;

			}
			//没有在上一个历元找到这颗卫星，可能是新升起的卫星，或者是参考星变化，重新初始化
			if (j == prev_nSats - 1 && !isMatch)
			{
				//重新进行初始化，phi矩阵对应位置默认0
				x_new_slip(rtkekf->nPos[CurIdx]) =
					((Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[ref].dL[0]) -
						(Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[ref].dP[0])) / w1;
				x_new_slip(rtkekf->nPos[CurIdx] + 1) = ((Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[ref].dL[1]) -
					(Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[ref].dP[1])) / w2;
				Q(rtkekf->nPos[CurIdx], rtkekf->nPos[CurIdx]) = 1E5;
				Q(rtkekf->nPos[CurIdx] + 1, rtkekf->nPos[CurIdx] + 1) = 1E5;
				break;
			}
		}
		CurIdx++;
	}

	for (int i = 0; i < 3; i++) Q(i, i) = 36;//设置SPP精度为6m
	for (int i = 0; i < 3; i++) G(i, i) = 0.5;


	X_predict = phi * X + x_new_slip;
	P_predict = phi * P_last * phi.transpose() + G * Q * G.transpose();

	for (int a = 0; a < 3 + (Raw->SdObs.SatNum - 2) * 2; a++)rtkekf->X0[a] = X_predict(a);
	int aP_predict = 0;
	for(int a=0;a< 3 + (Raw->SdObs.SatNum - 2) * 2;a++)
		for (int b = 0; b < 3 + (Raw->SdObs.SatNum - 2) * 2; b++)
		{
			rtkekf->P0[aP_predict] = P_predict(a, b);
			aP_predict++;
		}
	
}

void ObsUpdate(RAWDAT* Raw, RTKEKF* rtkekf)
{
	//**********************2.测量更新**********************
	MatrixXd X_predict = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 1);/*时间预测X*/
	MatrixXd P_predict = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);/*时间预测P*/
	MatrixXd X_update = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 1);/*测量更新X*/
	MatrixXd P_update = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);/*测量更新P*/
	MatrixXd I = MatrixXd::Identity(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);
	//B矩阵即为Hk
	MatrixXd B = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, 3 + (Raw->SdObs.SatNum - 2) * 2);
	//P矩阵为Rk的逆矩阵
	MatrixXd P = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, (Raw->SdObs.SatNum - 2) * 4);
	MatrixXd Rk = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, (Raw->SdObs.SatNum - 2) * 4);
	//W=Vk=Lk-Hk*Xk,k-1
	MatrixXd W = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 4, 1);
	

	for (int a = 0; a < 3 + (Raw->SdObs.SatNum - 2) * 2; a++)X_predict(a) = rtkekf->X0[a];
	int aP_predict = 0;
	for (int a = 0; a < 3 + (Raw->SdObs.SatNum - 2) * 2; a++)
		for (int b = 0; b < 3 + (Raw->SdObs.SatNum - 2) * 2; b++)
		{
			P_predict(a, b) = rtkekf->P0[aP_predict];
			aP_predict++;
		}

	//	 计算GPS和BDS双差卫星数（单双频混用：各系统每个频率的双差卫
	//	星数）
	//Raw->DDObs.Sats = Raw->SdObs.SatNum - 2;//单差总数-1：单频,单观测值
	Raw->DDObs.DDSatNum[0] = Raw->SdObs.SatNumGPS - 1;//GPS,单频
	Raw->DDObs.DDSatNum[1] = Raw->SdObs.SatNumBDS - 1;//BDS,单频
	//方向余弦
	double S_bas[MAXCHANNUM], S_rov[MAXCHANNUM], DDS[MAXCHANNUM], DDl[MAXCHANNUM], DDm[MAXCHANNUM], DDn[MAXCHANNUM];

	double X0_base[3];
	for (int i = 0; i < 3; i++)X0_base[i] = Raw->BasEpk.Pos[i];
	//	 计算基站坐标到所有卫星的几何距离
	double S_bas_ref[2];//0-->GPS,1-->BDS,到参考星距离
	S_bas_ref[0] = sqrt(pow(X0_base[0] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nBas].SatPos[0], 2) +
		pow(X0_base[1] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nBas].SatPos[1], 2) +
		pow(X0_base[2] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nBas].SatPos[2], 2));
	S_bas_ref[1] = sqrt(pow(X0_base[0] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nBas].SatPos[0], 2) +
		pow(X0_base[1] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nBas].SatPos[1], 2) +
		pow(X0_base[2] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nBas].SatPos[2], 2));

	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		if (!Raw->SdObs.SdSatObs[i].Valid)continue;
		S_bas[i] = sqrt(pow(X0_base[0] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[i].nBas].SatPos[0], 2) +
			pow(X0_base[1] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[i].nBas].SatPos[1], 2) +
			pow(X0_base[2] - Raw->BasEpk.SATPVT[Raw->SdObs.SdSatObs[i].nBas].SatPos[2], 2));
	}

	//	 计算流动站到参考星的几何距离
	double S_rov_ref[2];//0-->GPS,1-->BDS,到参考星距离
	S_rov_ref[0] =
		sqrt(pow(X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[0], 2) +
			pow(X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[1], 2) +
			pow(X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[2], 2));
	S_rov_ref[1] =
		sqrt(pow(X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[0], 2) +
			pow(X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[1], 2) +
			pow(X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[2], 2));
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		if (!Raw->SdObs.SdSatObs[i].Valid)continue;

		S_rov[i] =
			sqrt(pow(X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[0], 2) +
				pow(X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[1], 2) +
				pow(X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[2], 2));
	}

	int validIdx = 0; //
	for (int i = 0; i < Raw->SdObs.SatNum; i++)
	{
		//	① 参考星不用计算，存在半周不参与计算
		if (Raw->SdObs.SdSatObs[i].System == GPS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[0])continue;
		if (Raw->SdObs.SdSatObs[i].System == BDS && Raw->SdObs.SdSatObs[i].Prn == Raw->DDObs.RefPrn[1])continue;//参考星不参与计算
		if (!Raw->SdObs.SdSatObs[i].Valid)continue;//存在半周不参与计算;有周跳不参与计算
		//	② 线性化双差观测方程得到B矩阵和W向量
		if (Raw->SdObs.SdSatObs[i].System == GPS)
		{
			DDS[i] = S_rov[i] - S_rov_ref[0] - S_bas[i] + S_bas_ref[0];
			DDl[i] = (X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[0]) / S_rov[i] -
				(X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[0]) / S_rov_ref[0];
			DDm[i] = (X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[1]) / S_rov[i] -
				(X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[1]) / S_rov_ref[0];
			DDn[i] = (X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[2]) / S_rov[i] -
				(X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].nRov].SatPos[2]) / S_rov_ref[0];
			W(4 * validIdx) = Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[0] - DDS[i];
			W(4 * validIdx + 1) = Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dP[1] - DDS[i];
			W(4 * validIdx + 2) = Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[0] - DDS[i] - WL1_GPS * X_predict(2 * validIdx + 3);
			W(4 * validIdx + 3) = Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[0]].dL[1] - DDS[i] - WL2_GPS * X_predict(2 * validIdx + 4);
			B(4 * validIdx + 2, 2 * validIdx + 3) = WL1_GPS;//1-->2,3;1-->3,4;2-->6,5;7,6
			B(4 * validIdx + 3, 2 * validIdx + 4) = WL2_GPS;
		}
		else if (Raw->SdObs.SdSatObs[i].System == BDS)
		{
			DDS[i] = S_rov[i] - S_rov_ref[1] - S_bas[i] + S_bas_ref[1];
			DDl[i] = (X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[0]) / S_rov[i] -
				(X_predict(0) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[0]) / S_rov_ref[1];
			DDm[i] = (X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[1]) / S_rov[i] -
				(X_predict(1) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[1]) / S_rov_ref[1];
			DDn[i] = (X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[i].nRov].SatPos[2]) / S_rov[i] -
				(X_predict(2) - Raw->RovEpk.SATPVT[Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].nRov].SatPos[2]) / S_rov_ref[1];
			W(4 * validIdx) = Raw->SdObs.SdSatObs[i].dP[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[0] - DDS[i];
			W(4 * validIdx + 1) = Raw->SdObs.SdSatObs[i].dP[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dP[1] - DDS[i];
			W(4 * validIdx + 2) = Raw->SdObs.SdSatObs[i].dL[0] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[0] - DDS[i] - WL1_BDS * X_predict(2 * validIdx + 3);
			W(4 * validIdx + 3) = Raw->SdObs.SdSatObs[i].dL[1] - Raw->SdObs.SdSatObs[Raw->DDObs.RefPos[1]].dL[1] - DDS[i] - WL3_BDS * X_predict(2 * validIdx + 4);
			B(4 * validIdx + 2, 2 * validIdx + 3) = WL1_BDS;//1-->2,3;1-->3,4;2-->6,5;7,6
			B(4 * validIdx + 3, 2 * validIdx + 4) = WL3_BDS;
		}

		B(4 * validIdx, 0) = B(4 * validIdx + 1, 0) = B(4 * validIdx + 2, 0) = B(4 * validIdx + 3, 0) = DDl[i];
		B(4 * validIdx, 1) = B(4 * validIdx + 1, 1) = B(4 * validIdx + 2, 1) = B(4 * validIdx + 3, 1) = DDm[i];
		B(4 * validIdx, 2) = B(4 * validIdx + 1, 2) = B(4 * validIdx + 2, 2) = B(4 * validIdx + 3, 2) = DDn[i];


		validIdx++;//validIdx最大Raw->SdObs.SatNum-3
	}
	//	③ 计算权矩阵
	// 设置观测值权重
	double w_p = 2 * 0.3 * 0.3;   // 伪距权重因子
	double w_c = 2 * 0.01 * 0.01; // 载波相位权重因子

	// GPS系统权矩阵块
	if (Raw->DDObs.DDSatNum[0] > 0)
	{
		// 计算GPS系统的权重参数
		double diag[2], offdiag[2];
		diag[0] = Raw->DDObs.DDSatNum[0] / (w_p * (Raw->DDObs.DDSatNum[0] + 1));      // 伪距主对角线元素
		offdiag[0] = -1.0 / (w_p * (Raw->DDObs.DDSatNum[0] + 1)); // 伪距非对角线元素
		diag[1] = Raw->DDObs.DDSatNum[0] / (w_c * (Raw->DDObs.DDSatNum[0] + 1));       // 载波主对角线元素
		offdiag[1] = -1.0 / (w_c * (Raw->DDObs.DDSatNum[0] + 1));  // 载波非对角线元素

		for (int i = 0; i < Raw->DDObs.DDSatNum[0]; i++) {
			for (int j = 0; j < Raw->DDObs.DDSatNum[0]; j++) {
				// 计算权重值
				double value_p = (i == j) ? diag[0] : offdiag[0];
				double value_l = (i == j) ? diag[1] : offdiag[1];
				// 每个卫星的4个观测值
				int row = 4 * i;
				int col = 4 * j;
				P(row, col) = value_p;
				P(row + 1, col + 1) = value_p;
				P(row + 2, col + 2) = value_l;
				P(row + 3, col + 3) = value_l;
			}
		}
	}



	// BDS系统权矩阵块
	if (Raw->DDObs.DDSatNum[1] > 0)
	{
		// 计算BDS系统的权重参数
		double diag[2], offdiag[2];
		diag[0] = Raw->DDObs.DDSatNum[1] / (w_p * (Raw->DDObs.DDSatNum[1] + 1));      // 伪距主对角线元素
		offdiag[0] = -1.0 / (w_p * (Raw->DDObs.DDSatNum[1] + 1)); // 伪距非对角线元素
		diag[1] = Raw->DDObs.DDSatNum[1] / (w_c * (Raw->DDObs.DDSatNum[1] + 1));       // 载波主对角线元素
		offdiag[1] = -1.0 / (w_c * (Raw->DDObs.DDSatNum[1] + 1));  // 载波非对角线元素


		int bds_idx = 4 * Raw->DDObs.DDSatNum[0];

		for (int i = 0; i < Raw->DDObs.DDSatNum[1]; i++) {
			for (int j = 0; j < Raw->DDObs.DDSatNum[1]; j++) {
				// 计算权重值
				double value_p = (i == j) ? diag[0] : offdiag[0];
				double value_l = (i == j) ? diag[1] : offdiag[1];

				int row = bds_idx + 4 * i;
				int col = bds_idx + 4 * j;
				P(row, col) = value_p;
				P(row + 1, col + 1) = value_p;
				P(row + 2, col + 2) = value_l;
				P(row + 3, col + 3) = value_l;
			}
		}
	}

	Rk = P.inverse();


	//Kalman增益
	MatrixXd Kk;
	Kk = P_predict * B.transpose() * (B * P_predict * B.transpose() + Rk).inverse();
	X_update = X_predict + Kk * W;
	P_update = (I - Kk * B) * P_predict;
	
	//保存结果
	 // 保存当前历元观测数据到RTKEKF
	rtkekf->DDObs = Raw->DDObs;  // 保存当前双差观测值
	rtkekf->SDObs = Raw->SdObs;     // 保存当前单差观测值
	for (int i = 0; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)
	{
		rtkekf->X[i] = X_update(i);
	}

	

	//rtkekf->P = P_update;
	int ppidx = 0;
	for (int i = 0; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)
		for (int j = 0; j < 3 + (Raw->SdObs.SatNum - 2) * 2; j++)//?
		{
			rtkekf->P[ppidx] = P_update(j, i);
			ppidx++;
		}

	rtkekf->nSats = Raw->SdObs.SatNum - 2;

	int Qxxidx = 0;
	/*for (int j = 0; j < 3; j++)*/ //按列保存Qxx，（0,0），（1,0）,（2,0），（0,1）......
	  // 按列循环
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
		{  // 按行循环
			Raw->Res.Qxx[Qxxidx++] = P_update(i, j);
		}

	int QnnIdx = 0;
	/*for (int i = 0; i < (Raw->SdObs.SatNum - 2) * 2; i++)*///按列保存Qnn
	for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
		for (int i = 0; i < (Raw->SdObs.SatNum - 2) * 2; i++)
		{
			Raw->Res.Qnn[QnnIdx] = P_update(j + 3, i + 3);
			QnnIdx++;
		}

	int QxnIdx = 0;
	/*for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)*///按列保存Qxn
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
		{
			Raw->Res.Qxn[QxnIdx] = P_update(i, j + 3);
			QxnIdx++;
		}
}


//Kalman滤波求解
bool RTKFloatKalman(RAWDAT* Raw, PPRESULT* Base, PPRESULT* Rov, RTKEKF* rtkekf)
{
	//卫星数太少不定位；单系统不定位
	if (Raw->SdObs.SatNum < 6)return false;
	if (Raw->SdObs.SatNumBDS < 2 || Raw->SdObs.SatNumGPS < 2)return false;
	//参考星选取不成功
	if (Raw->DDObs.RefPos[0] == -1 || Raw->DDObs.RefPos[1] == -1 || Raw->DDObs.RefPrn[0] == -1 || Raw->DDObs.RefPrn[1] == -1)return false;

	//*************************1.初始化******************************
	//是否初始化
	//参考星变化应该重新初始化
	//上个历元如果没有固定，重新初始化
	if (!rtkekf->IsInit||Raw->DDObs.RefValid|| !Raw->DDObs.bFixed)
	{
		KfInit(Raw, Base, Rov, rtkekf);
		//return false;
	}



	//*************************2.时间预测******************************
	TimePredict(Raw, rtkekf);
	//*************************3.测量更新******************************
	ObsUpdate(Raw, rtkekf);

	//*************************4.模糊度固定******************************

	MatrixXd X_update = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 1);/*测量更新X*/
	MatrixXd P_update = MatrixXd::Zero(3 + (Raw->SdObs.SatNum - 2) * 2, 3 + (Raw->SdObs.SatNum - 2) * 2);/*测量更新P*/
	MatrixXd N = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 2, 1);
	MatrixXd QNN = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 2, (Raw->SdObs.SatNum - 2) * 2);
	MatrixXd QXN = MatrixXd::Zero(3, (Raw->SdObs.SatNum - 2) * 2);
	MatrixXd X_fixed = MatrixXd::Zero(3, 1);
	MatrixXd deltaX = MatrixXd::Zero(3, 1);

	for (int i = 0; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)
	{
		X_update(i) = rtkekf->X[i];
	}



	//rtkekf->P = P_update;
	int ppidx = 0;
	for (int i = 0; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)
		for (int j = 0; j < 3 + (Raw->SdObs.SatNum - 2) * 2; j++)//?
		{
			P_update(j, i) = rtkekf->P[ppidx];
			ppidx++;
		}

	for (int i = 0; i < 3; i++)X_fixed(i) = X_update(i);//暂存，待更新
	

	for (int a = 0; a < (Raw->SdObs.SatNum - 2) * 2; a++)
	{
		Raw->Res.Amb[a] = X_update(3 + a);//模糊度，未固定
		N(a) = X_update(3 + a);
	}

	for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
		for (int i = 0; i < (Raw->SdObs.SatNum - 2) * 2; i++)
		{
			QNN(j, i) = P_update(j + 3, i + 3);
		}

	/*for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)*///按列保存Qxn
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < (Raw->SdObs.SatNum - 2) * 2; j++)
		{
			QXN(i, j) = P_update(i, j + 3);
		}

	MatrixXd FixedN = MatrixXd::Zero((Raw->SdObs.SatNum - 2) * 2, 1);
	Raw->Res.isFix = lambda((Raw->SdObs.SatNum - 2) * 2, 2, Raw->Res.Amb, Raw->Res.Qnn, Raw->DDObs.FixedAmb, Raw->DDObs.ResAmb);
	if (Raw->Res.isFix != 0) {
		Raw->DDObs.bFixed = false;
		return false;//没固定
	}
	if (Raw->DDObs.ResAmb[0] < 1E-10)Raw->DDObs.ResAmb[0] = 1E-6;
	Raw->DDObs.Ratio = Raw->Res.Ratio = Raw->DDObs.ResAmb[1] / Raw->DDObs.ResAmb[0];
	if (Raw->DDObs.Ratio <= 3.0)
	{
		/*for (int a = 0; a < (Raw->SdObs.SatNum - 2) * 2; a++)
			FixedN(a) = Raw->DDObs.FixedAmb[a + (Raw->SdObs.SatNum - 2) * 2];*/
		for (int a = 0; a < 3; a++)Raw->DDObs.Position[a] = X_update(a);
		for (int a = 0; a < 3; a++)Raw->DDObs.dPos[a] = X_update(a) - Raw->BasEpk.Pos[a];
		Raw->DDObs.bFixed = false;
		return true;//直接返回浮点解
	}
	else {
		for (int a = 0; a < (Raw->SdObs.SatNum - 2) * 2; a++)FixedN(a) = Raw->DDObs.FixedAmb[a];
		deltaX = QXN * QNN.inverse() * (FixedN - N);
		X_fixed = X_fixed + QXN * QNN.inverse() * (FixedN - N);
		for (int a = 0; a < 3; a++)Raw->DDObs.dPos[a] = X_fixed(a) - Raw->BasEpk.Pos[a];
		for (int a = 0; a < 3; a++)Raw->DDObs.Position[a] = X_fixed(a);
		//保存固定解
		for (int i = 0; i < 3 + (Raw->SdObs.SatNum - 2) * 2; i++)
		{
			if (i < 3)
				rtkekf->X[i] = X_fixed(i);
			else
				rtkekf->X[i] = FixedN(i - 3);
		}
		Raw->DDObs.bFixed = true;
	}
	
	
	return true;
}