#include<string>
#include"RTK.h"
void DetectOutlier(EPOCHOBSDATA* Obs)
{
    MWGF CurComObs[MAXCHANNUM];  // 当前历元计算结果

    for (int i = 0; i < Obs->SatNum; i++)
    {
        //判断历元间相位数据的连续性，未发生中断时，Locktime(k+1)>=Locktime(k),数据正常，反之有周跳
        Obs->ComObs[i].LockTime[0] = Obs->SatObs[i].LockTime[0];
        Obs->ComObs[i].LockTime[1] = Obs->SatObs[i].LockTime[1];
        if (Obs->SatObs[i].LockTime[0] < Obs->ComObs[i].LockTime[0] && Obs->SatObs[i].LockTime[1] < Obs->ComObs[i].LockTime[1])
        {
            Obs->SatObs[i].Valid = false;
            continue;
        }
        
        ////Parity：判断相位数据是否存在半周，Parity=0，可能存在半周
        //if(Obs->SatObs[i].half[0]==0&& Obs->SatObs[i].half[1]==0)
        //{
        //    Obs->SatObs[i].Valid = false;
        //    continue;
        //}
        
        // 1. 检查该卫星的双频伪距和相位数据是否有效和完整，若不全或为0，将Valid标记为false，
        //continue
        if (fabs(Obs->SatObs[i].P[0]) < 1E-3 || fabs(Obs->SatObs[i].P[1]) < 1E-3 ||
            fabs(Obs->SatObs[i].L[0]) < 1E-3 || fabs(Obs->SatObs[i].L[1]) < 1E-3)
        {
            Obs->SatObs[i].Valid = false;
            continue;
        }

        // 2. 计算当前历元该卫星的GF和MW组合值
        double f1, f2;
        switch (Obs->SatObs[i].System) 
        {
        case GPS:
            f1 = FG1_GPS;
            f2 = FG2_GPS;
            break;
        case BDS:
            f1 = FG1_BDS;
            f2 = FG3_BDS;
            break;
        default:
            Obs->SatObs[i].Valid = false;
            continue;
        }

        // 计算当前MW和GF（未平滑）
        CurComObs[i].Prn = Obs->SatObs[i].Prn;
        CurComObs[i].Sys = Obs->SatObs[i].System;
        CurComObs[i].MW = (f1 * Obs->SatObs[i].L[0] - f2 * Obs->SatObs[i].L[1]) / (f1 - f2)
            - (f1 * Obs->SatObs[i].P[0] + f2 * Obs->SatObs[i].P[1]) / (f1 + f2);
        CurComObs[i].GF = Obs->SatObs[i].L[0] - Obs->SatObs[i].L[1];

        // 3. 从上个历元的MWGF数据中查找该卫星的GF和MW组合值
        MWGF* prev = nullptr;
        for (int j = 0; j < MAXCHANNUM; j++)
        {
            if (Obs->ComObs[j].Prn == Obs->SatObs[i].Prn &&
                Obs->ComObs[j].Sys == Obs->SatObs[i].System)
            {
                prev = &Obs->ComObs[j];  // 找到匹配项
                break;
            }
        }

        // 4. 计算当前历元该卫星GF与上一历元对应GF的差值dGF
        //5. 计算当前历元该卫星MW与上一历元对应MW平滑值的差值dMW
        double dGF = 0.0, dMW = 0.0;
        if (prev != nullptr)  // 找到历史数据才计算差值
        {
            dGF = fabs(CurComObs[i].GF - prev->GF);
            dMW = fabs(CurComObs[i].MW - prev->MW);  // 注意比较的是平滑后的MW
        }

        // 6. 检查dGF和dMW是否超限，限差阈值建议为5cm和3m。若超限，标记为粗差，将Valid标记为
        //false ，若不超限，标记为可用将Valid标记为true，并计算该卫星的MW平滑值
        if (dGF > 0.05 || dMW > 3.0)
        {
            Obs->SatObs[i].Valid = false;  // 超限标记为粗差
        }
        else
        {
            Obs->SatObs[i].Valid = true;

            // 更新MW平滑值
            if (prev != nullptr)
            {
                // 连续跟踪：递推平滑
                CurComObs[i].n = prev->n + 1;
                CurComObs[i].MW = (prev->MW * prev->n + CurComObs[i].MW) / CurComObs[i].n;
            }
            else
            {
                CurComObs[i].n = 1;
            }
        }

        // 7. 对于可用的观测数据，计算伪距的IF组合观测值，用于SPP
        if (Obs->SatObs[i].Valid)
        {
            double f1_sq = f1 * f1;
            double f2_sq = f2 * f2;
            CurComObs[i].PIF = (f1_sq * Obs->SatObs[i].P[0] - f2_sq * Obs->SatObs[i].P[1])
                / (f1_sq - f2_sq);
        }
    }

    // 8. 所有卫星循环计算完成之后，将CurComObs内存拷贝到ComObs，即函数运行结束后， ComObs
    //保存了当前历元的GF和MW平滑值。
    memcpy(Obs->ComObs, CurComObs, sizeof(MWGF) * MAXCHANNUM);
}


void DetectCycleSlip(SDEPOCHOBS* Obs)
{
    MWGF CurComObs[MAXCHANNUM];  // 当前历元计算结果
    
    for (int i = 0; i < Obs->SatNum; i++)
    {
        //// 1. 检查该卫星的双频伪距和相位数据是否有效和完整，若不全或为0，将Valid标记为false，
        ////continue
        //if (fabs(Obs->SdSatObs[i].dP[0]) < 1E-3 || fabs(Obs->SdSatObs[i].dP[1]) < 1E-3 ||
        //    fabs(Obs->SdSatObs[i].dL[0]) < 1E-3 || fabs(Obs->SdSatObs[i].dL[1]) < 1E-3)
        //{
        //    Obs->SdSatObs[i].Valid = false;
        //    Obs->SatNum = Obs->SatNum - 1;//只记录有效卫星数
        //    continue;
        //}

        // 2. 计算当前历元该卫星的GF和MW组合值
        double f1, f2;
        switch (Obs->SdSatObs[i].System)
        {
        case GPS:
            f1 = FG1_GPS;
            f2 = FG2_GPS;
            break;
        case BDS:
            f1 = FG1_BDS;
            f2 = FG3_BDS;
            break;
        default:
            Obs->SdSatObs[i].Valid = false;
            Obs->SatNum = Obs->SatNum - 1;//只记录有效卫星数
            continue;
        }

        // 计算当前MW和GF（未平滑）
        CurComObs[i].Prn = Obs->SdSatObs[i].Prn;
        CurComObs[i].Sys = Obs->SdSatObs[i].System;
        CurComObs[i].MW = (f1 * Obs->SdSatObs[i].dL[0] - f2 * Obs->SdSatObs[i].dL[1]) / (f1 - f2)
            - (f1 * Obs->SdSatObs[i].dP[0] + f2 * Obs->SdSatObs[i].dP[1]) / (f1 + f2);
        CurComObs[i].GF = Obs->SdSatObs[i].dL[0] - Obs->SdSatObs[i].dL[1];

        // 3. 从上个历元的MWGF数据中查找该卫星的GF和MW组合值
        MWGF* prev = nullptr;
        for (int j = 0; j < MAXCHANNUM; j++)
        {
            if (Obs->SdCObs[j].Prn == Obs->SdSatObs[i].Prn &&
                Obs->SdCObs[j].Sys == Obs->SdSatObs[i].System)
            {
                prev = &Obs->SdCObs[j];  // 找到匹配项
                break;
            }
        }

        // 4. 计算当前历元该卫星GF与上一历元对应GF的差值dGF
        //5. 计算当前历元该卫星MW与上一历元对应MW平滑值的差值dMW
        double dGF = 0.0, dMW = 0.0;
        if (prev != nullptr)  // 找到历史数据才计算差值
        {
            dGF = fabs(CurComObs[i].GF - prev->GF);
            dMW = fabs(CurComObs[i].MW - prev->MW);  // 注意比较的是平滑后的MW
        }

        // 6. 检查dGF和dMW是否超限，限差阈值建议为5cm和3m。若超限，标记为粗差，将Valid标记为
        //false ，若不超限，标记为可用将Valid标记为true，并计算该卫星的MW平滑值
        if (dGF > 0.05 || dMW > 3.0)
        {
            Obs->SdSatObs[i].Valid = false;  // 超限标记为粗差
            Obs->SatNum = Obs->SatNum - 1;//只记录有效卫星数
        }
        else
        {
            Obs->SdSatObs[i].Valid = true;

            // 更新MW平滑值
            if (prev != nullptr)
            {
                // 连续跟踪：递推平滑
                CurComObs[i].n = prev->n + 1;
                CurComObs[i].MW = (prev->MW * prev->n + CurComObs[i].MW) / CurComObs[i].n;
            }
            else
            {
                CurComObs[i].n = 1;
            }
        }

        
    }

    // 8. 所有卫星循环计算完成之后，将CurComObs内存拷贝到ComObs，即函数运行结束后， ComObs
    //保存了当前历元的GF和MW平滑值。
    memcpy(Obs->SdCObs, CurComObs, sizeof(MWGF) * MAXCHANNUM);
}