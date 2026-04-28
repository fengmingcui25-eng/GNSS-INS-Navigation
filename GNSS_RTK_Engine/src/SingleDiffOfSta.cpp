#include"RTK.h"

void SingleDiff(EPOCHOBSDATA* Obs_rover, EPOCHOBSDATA* Obs_base, SDEPOCHOBS* SDRes)
{
    
	//根据流动站卫星数进行循环，对于每颗流动站的卫星
    int satnum = 0, satnumG = 0, satnumB = 0;
    for (int i = 0; i < Obs_rover->SatNum; i++)
    {
        //判断历元间相位数据的连续性，未发生中断时，Locktime(k+1)>=Locktime(k),数据正常，反之有周跳
       /* if (Obs_base->SatObs[i].LockTime[0] < OBas.SatObs[i].LockTime[0] || Obs_base->SatObs[i].LockTime[1] < OBas.SatObs[i].LockTime[1]
            || Obs_rover->SatObs[i].LockTime[0] < ORov.SatObs[i].LockTime[0] || Obs_rover->SatObs[i].LockTime[1] < ORov.SatObs[i].LockTime[1])
        {
            SDRes->SdSatObs[i].Valid = false;
            continue;
        }*/
        //Parity：判断相位数据是否存在半周，Parity=0，可能存在半周
        if (Obs_rover->SatObs[i].half[0] == 0|| Obs_rover->SatObs[i].half[1] == 0)
        {
            SDRes->SdSatObs[i].Valid = false;
            continue;
        }
       
        //1. 检查卫星号和系统号是否正常，卫星观测值是否完整，如果不正常，continue；

        if (fabs(Obs_rover->SatObs[i].P[0]) < 1E-3 || fabs(Obs_rover->SatObs[i].P[1]) < 1E-3 ||
            fabs(Obs_rover->SatObs[i].L[0]) < 1E-3 || fabs(Obs_rover->SatObs[i].L[1]) < 1E-3)
        {
            Obs_rover->SatObs[i].Valid = false;
            continue;
        }
        //检查高度角，如果高度角太低，不采用
        if (Obs_rover->SATPVT[i].Elevation * 180.0 / pi < 20)continue;
        //流动站观测到的卫星Locktime小于6秒，意味着刚升起，不参与模糊度固定
        if (Obs_rover->SatObs[i].LockTime[0] < 6 || Obs_rover->SatObs[i].LockTime[1] < 6)continue;
        //卫星位置计算是否成功
        if (!Obs_rover->SATPVT[i].Valid)continue;
        //2. 在基站观测值中查找与该流动站相同的卫星，如果未找到，continue；
        //SDRes->SatNum = SDRes->SatNumGPS = SDRes->SatNumBDS = 0;
        if (!Obs_rover->SatObs[i].Valid)continue;
        for (int j = 0; j < Obs_base->SatNum; j++)
        {
            if (Obs_base->SatObs[j].half[0] == 0 || Obs_base->SatObs[j].half[1] == 0)
            {
                SDRes->SdSatObs[j].Valid = false;
                continue;
            }
            if (fabs(Obs_base->SatObs[j].P[0]) < 1E-3 || fabs(Obs_base->SatObs[j].P[1]) < 1E-3 ||
                fabs(Obs_base->SatObs[j].L[0]) < 1E-3 || fabs(Obs_base->SatObs[j].L[1]) < 1E-3)
            {
                Obs_base->SatObs[j].Valid = false;
                continue;
            }
            if (!Obs_base->SATPVT[j].Valid)continue;
            if (Obs_base->SatObs[j].System == Obs_rover->SatObs[i].System
                && Obs_base->SatObs[j].Prn == Obs_rover->SatObs[i].Prn)
            {
                //3. 对同类型和同频率的观测值求差并保存，保存索引号、卫星号和系统号等相关信息，累加单差观测值卫星数量
                if (Obs_base->SatObs[j].Valid)//单差无周跳
                {
                    SDRes->SdSatObs[satnum].dP[0] = Obs_rover->SatObs[i].P[0] - Obs_base->SatObs[j].P[0];
                    SDRes->SdSatObs[satnum].dP[1] = Obs_rover->SatObs[i].P[1] - Obs_base->SatObs[j].P[1];
                    SDRes->SdSatObs[satnum].dL[0] = Obs_rover->SatObs[i].L[0] - Obs_base->SatObs[j].L[0];
                    SDRes->SdSatObs[satnum].dL[1] = Obs_rover->SatObs[i].L[1] - Obs_base->SatObs[j].L[1];
                    SDRes->SdSatObs[satnum].nBas = j;
                    SDRes->SdSatObs[satnum].nRov = i;
                    SDRes->SdSatObs[satnum].System = Obs_rover->SatObs[i].System;
                    SDRes->SdSatObs[satnum].Prn = Obs_rover->SatObs[i].Prn;
                    //SDRes->SdSatObs[satnum].Valid = true;
                    satnum++;//只表示单频
                    if (Obs_rover->SatObs[i].System == GPS) satnumG++;
                    else if (Obs_rover->SatObs[i].System == BDS) satnumB++;
                    
                    
                    break;
                }
            }
            else continue;
        }
       
    }
    SDRes->SatNum = satnum; SDRes->SatNumGPS = satnumG; SDRes->SatNumBDS = satnumB;
    //4. 循环结束后，对单差观测值的观测时刻和卫星数量赋值。
    SDRes->Time.Week = Obs_rover->Time.Week;
    SDRes->Time.SecOfWeek = Obs_rover->Time.SecOfWeek;
    
}