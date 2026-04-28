#include"RTK.h"
#include<iostream>
using namespace std;

void DetRefSat(const EPOCHOBSDATA* epkA, const EPOCHOBSDATA* epkB, SDEPOCHOBS* SDObs, DDCOBS* DDObs)
{
    // 存储上一历元的参考星PRN号
    int prev_ref_prn[2] = { DDObs->RefPrn[0], DDObs->RefPrn[1] };
    bool ref_found[2] = { false, false }; // 标记是否已为各系统找到参考星

    // 第一部分：尝试继承上一历元的参考星（通过PRN匹配）
    for (int i = 0; i < SDObs->SatNum; ++i) {
        // 确定卫星系统
        int sys_idx = (SDObs->SdSatObs[i].System == GPS) ? 0 : 1;

        // 如果该卫星不是上一历元的参考星，跳过
        if (SDObs->SdSatObs[i].Prn != prev_ref_prn[sys_idx]) continue;

        // 获取流动站和基准站原始观测值的索引
        int nRov = SDObs->SdSatObs[i].nRov;
        int nBas = SDObs->SdSatObs[i].nBas;

        // 检查所有条件是否满足
        if (SDObs->SdSatObs[i].Valid &&                                       // 周跳探测
            epkB->SATPVT[nRov].Elevation > 30.0 / 180.0 * pi &&                  // 流动站高度角 > 30度
            epkA->SATPVT[nBas].Elevation > 30.0 / 180.0 * pi)                     // 基准站高度角 > 30度
        {
            // 条件全部满足，继承该参考星，参考星未发生改变
            DDObs->RefPos[sys_idx] = i;  // 更新为当前历元的索引
            DDObs->RefPrn[sys_idx] = SDObs->SdSatObs[i].Prn;
            ref_found[sys_idx] = true;
            DDObs->RefValid = false;
        }
    }

    // 如果两个系统都成功继承了参考星，直接返回
    if (ref_found[0] && ref_found[1]) {
        return;
    }


    int i, j, n;
    double Sum[2] = { 0.0 }, MaxSum[2] = { 0.0 };
    int RefPrn[2] = { -1 }, RefIndex[2] = { -1 };

    /*基准星选取*/
    for (int i = 0; i < SDObs->SatNum; i++)
    {
        if (!SDObs->SdSatObs[i].Valid || !epkA->SATPVT[SDObs->SdSatObs[i].nBas].Valid || !epkB->SATPVT[SDObs->SdSatObs[i].nRov].Valid) continue;
        if (epkA->SatObs[SDObs->SdSatObs[i].nBas].LockTime[0] < 6 || epkA->SatObs[SDObs->SdSatObs[i].nBas].LockTime[1] < 6 || 
            epkB->SatObs[SDObs->SdSatObs[i].nRov].LockTime[0] < 6 || epkB->SatObs[SDObs->SdSatObs[i].nRov].LockTime[1] < 6) continue;

        n = SDObs->SdSatObs[i].System == GPS ? 0 : 1;
        Sum[n] = epkB->SATPVT[SDObs->SdSatObs[i].nRov].Elevation*180.0/pi + epkA->SATPVT[SDObs->SdSatObs[i].nBas].Elevation * 180.0 / pi + epkA->SatObs[SDObs->SdSatObs[i].nBas].cn0[0] +
            epkA->SatObs[SDObs->SdSatObs[i].nBas].cn0[1]+ epkB->SatObs[SDObs->SdSatObs[i].nRov].cn0[0]+ epkB->SatObs[SDObs->SdSatObs[i].nRov].cn0[1];

        if (Sum[n] > MaxSum[n]) {
            MaxSum[n] = Sum[n];
            RefPrn[n] = SDObs->SdSatObs[i].Prn;
            RefIndex[n] = i;
        }
    }

    for (j = 0; j < 2; j++) {
        if (MaxSum[j] < 200.0) DDObs->RefPos[j] = -1;
        else
        {
            DDObs->RefPos[j] = RefIndex[j];
            DDObs->RefPrn[j] = RefPrn[j];
           
        }
    }
    //参考星发生变化
    DDObs->RefValid = true;
    //cout << "参考星变化" << "\t" << "原参考星" << prev_ref_prn[0] << "  " << prev_ref_prn[1] << "新参考星" << DDObs->RefPrn[0] << "  " << DDObs->RefPrn[1] << endl;
}

