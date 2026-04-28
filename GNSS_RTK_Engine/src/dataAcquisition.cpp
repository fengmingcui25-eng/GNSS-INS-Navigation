#include <stdio.h>
#include <stdint.h>
#include<string.h>
#include<cmath>
#include"RTK.h"
double R8(unsigned char* p)
{
    double r;
    memcpy(&r, p, 8);
    return r;
}
float R4(unsigned char* p)
{
    float r;
    memcpy(&r, p, 4);
    return r;
}
int I4(unsigned char* p)
{
    int r;
    memcpy(&r, p, 4);
    return r;
}
unsigned int UI4(unsigned char* p)
{
    unsigned int r;
    memcpy(&r, p, 4);
    return r;
}
short I2(unsigned char* p)
{
    short r;
    memcpy(&r, p, 2);
    return r;
}
unsigned short UI2(unsigned char* p)
{
    unsigned short r;
    memcpy(&r, p, 2);
    return r;
}

// CRC32校验函数
unsigned int crc32(const unsigned char* buff, int len)
{
    int i, j;
    unsigned int crc = 0;
    for (i = 0; i < len; i++)
    {
        crc ^= buff[i];
        for (j = 0; j < 8; j++)
        {
            if (crc & 1)
                crc = (crc >> 1) ^ POLYCRC32;
            else
                crc >>= 1;
        }
    }
    return crc;
}

void DecodeRange(unsigned char* data, int len, EPOCHOBSDATA* obs)
{
    int i, j, n, k, ObsNum, Freq, Prn;
    GNSSSys sys;
    double wl;
    int PhaseLockFlag, CodeLockedFlag, ParityFlag, SatSystem, SigType;
    unsigned int ChanStatus;
    unsigned char* p = data + 28;

    //1. 从消息头中解码得到观测时刻，该时刻为接收机钟表面时，用GPSTIME结构体表示。
    obs->Time.Week = UI2(data + 14);
    obs->Time.SecOfWeek = UI4(data + 16)*1E-3;
    //    2. 解码得到观测值数量，为所有卫星所有信号观测值的总数。
    ObsNum = UI4(p);  // 观测卫星数量
    memset(obs->SatObs, 0, MAXCHANNUM * sizeof(SATOBSDATA));
    //    3. 对所有信号观测值进行循环解码
    for (i = 0, p += 4; i < ObsNum; i++, p += 44) {
        //    ① 解码得到跟踪状态标记，从中取出Phase lock flag / Code locked
        //    flag / Parity known flag / Satellite system / signal type等数据
        ChanStatus = UI4(p + 40);
        ParityFlag = (ChanStatus >> 11) & 0x01;
        PhaseLockFlag = (ChanStatus >> 10) & 0x01;
        CodeLockedFlag = (ChanStatus >> 12) & 0x01;
        SatSystem = (ChanStatus >> 16) & 0x07;
        SigType = (ChanStatus >> 21) & 0x1F;
        //    ② 如果卫星系统不是GPS或BDS， continue至①

    //    ③ 如果GPS卫星的信号类型不是L1 C / A或者L2P（Y）， BDS卫星不是B1I
    //    或B3I， continue至①，并记录信号频率类型，第一频率s = 0，第二频
    //    率s = 1；
        if (SatSystem == 0) {
            sys = GPS;
            if (SigType == 0) {
                Freq = 0; wl = WL1_GPS;
            }
            else if (SigType == 9) { Freq = 1; wl = WL2_GPS; }
            else continue;
        }
        else if (SatSystem == 4) {
            sys = BDS;
            if (SigType == 0 || SigType == 4) {
                Freq = 0; wl = WL1_BDS;
            }
            else if (SigType == 2 || SigType == 6) { Freq = 1; wl = WL3_BDS; }
            else continue;
        }
        else continue;

        //    ④ 解码得到卫星号Prn以及卫星系统号，在当前观测值结构体中进行搜索，
        //    如果找到相同的卫星，将解码的观测值填充到该卫星对应的数组中；如
        //    果在当前已解码的卫星数据中没有发现，则填充到现有数据的末尾。
        Prn = UI2(p);
        for (j = 0; j < MAXCHANNUM; j++) {
            if (obs->SatObs[j].System == sys && obs->SatObs[j].Prn == Prn) {
                n = j; break;
            }
            if (obs->SatObs[j].Prn == 0) {
                k = n = j; break;
            }
        }
        obs->SatObs[n].Prn = Prn;
        obs->SatObs[n].System = sys;
        obs->SatObs[n].P[Freq] = CodeLockedFlag == 1 ? R8(p + 4) : 0.0;
        obs->SatObs[n].L[Freq] = -wl * (PhaseLockFlag == 1 ? R8(p + 16) : 0.0);
        obs->SatObs[n].D[Freq] = -wl * R4(p + 28);
        obs->SatObs[n].cn0[Freq] = R4(p + 32);
        obs->SatObs[n].LockTime[Freq] = R4(p + 36);
        obs->SatObs[n].half[Freq] = ParityFlag;
    }
    obs->SatNum = k+1;
}
//
void DecodeGpsEphem(unsigned char* data, int len, GPSEPHREC geph[])
{
    unsigned char* p = data + 28; // 跳过28字节消息头
    GPSEPHREC eph;
    memset(&eph, 0, sizeof(GPSEPHREC));
    eph.System = GPS;

    eph.PRN = UI4(p);       // PRN号 (4字节)
    eph.SVHealth = UI4(p + 12);//// 卫星健康状态 (0=健康)
    eph.IODE= UI4(p + 16);     // IODE 
    eph.TOE.Week = UI4(p + 24);// TOE的GPS周
    eph.TOC.Week = eph.TOE.Week;// TOC的GPS周
    eph.TOE.SecOfWeek = R8(p + 32); // TOE的周内秒
    eph.SqrtA = sqrt(R8(p + 40));// 轨道长半轴平方根 (√m)//////////?
    eph.DetlaN = R8(p + 48);// 平均运动修正量 (rad/s)
    eph.M0 = R8(p + 56);// 平近点角 (弧度)
    eph.e = R8(p + 64);// 偏心率
    eph.omega= R8(p + 72);// 近地点角距 (弧度)
    eph.Cuc = R8(p + 80);
    eph.Cus = R8(p + 88);
    eph.Crc = R8(p + 96);
    eph.Crs = R8(p + 104);
    eph.Cic = R8(p + 112);
    eph.Cis = R8(p + 120);
    eph.i0 = R8(p + 128);// 轨道倾角 (弧度)
    eph.iDot = R8(p + 136);// 轨道倾角变化率 (rad/s)
    eph.OMEGA = R8(p + 144);// 升交点赤经 (弧度)
    eph.OMEGADot = R8(p + 152);// 升交点赤经变化率 (rad/s)
    eph.IODC = R8(p + 160);     // IODC 
    eph.TOC.SecOfWeek = R8(p + 164);// TOC的周内秒
    eph.TGD1 = R8(p + 172);     // 群延差 TGD (秒)
    eph.ClkBias = R8(p + 180);  // 卫星钟差 (秒)
    eph.ClkDrift = R8(p + 188); // 卫星钟速 (秒/秒)
    eph.ClkDriftRate = R8(p + 196); // 卫星钟漂 (秒/秒²)
   

    // 存入数组（PRN范围为1-32）
    if (eph.PRN >= 1 && eph.PRN <= MAXGPSNUM) {
        geph[eph.PRN - 1] = eph;
    }
}
//
void DecodeBdsEphem(unsigned char* data, int len, GPSEPHREC beph[])
{
    unsigned char* p = data + 28; // 跳过28字节消息头
    GPSEPHREC eph;
    memset(&eph, 0, sizeof(GPSEPHREC));
    eph.System = BDS;

    eph.PRN = UI4(p);       // PRN号 (4字节)
    eph.TOE.Week = UI4(p + 4);
    eph.TOC.Week = eph.TOE.Week;
    eph.SVHealth = UI4(p + 16);//// 卫星健康状态 (0=健康)
    eph.TGD1 = R8(p + 20);
    eph.TGD2 = R8(p + 28);
    eph.TOC.SecOfWeek = UI4(p + 40);

    eph.ClkBias = R8(p + 44);  // 卫星钟差 (秒)
    eph.ClkDrift = R8(p + 52); // 卫星钟速 (秒/秒)
    eph.ClkDriftRate = R8(p + 60); // 卫星钟漂 (秒/秒²)

    eph.TOE.SecOfWeek = UI4(p + 72);
    eph.SqrtA = R8(p + 76);
    eph.e = R8(p + 84);// 偏心率
    eph.omega = R8(p + 92);// 近地点角距 (弧度)
    eph.DetlaN = R8(p + 100);// 平均运动修正量 (rad/s)
    eph.M0 = R8(p + 108);// 平近点角 (弧度)
    eph.OMEGA = R8(p + 116);// 升交点赤经 (弧度)
    eph.OMEGADot = R8(p + 124);// 升交点赤经变化率 (rad/s)
    eph.i0 = R8(p + 132);// 轨道倾角 (弧度)
    eph.iDot = R8(p + 140);// 轨道倾角变化率 (rad/s)
    eph.Cuc = R8(p + 148);
    eph.Cus = R8(p + 156);
    eph.Crc = R8(p + 164);
    eph.Crs = R8(p + 172);
    eph.Cic = R8(p + 180);
    eph.Cis = R8(p + 188);

    
   
    // 存入数组（BDS PRN范围为1-63）
    if (eph.PRN >= 1 && eph.PRN <= MAXBDSNUM) {
        beph[eph.PRN - 1] = eph;
    }
}

int Decode_psrpos(unsigned char* data, int len, POSRES* pos)
{
    pos->Time.Week = I2(data+14);         
    pos->Time.SecOfWeek = (double)(UI4(data + 16)*1E-3); 
    unsigned char* p = data + 28; // 跳过28字节消息头
    memset(pos, 0, sizeof(POSRES));

    pos->Pos[0] = R8(p + 8);        // 纬度 Latitude (deg)
    pos->Pos[1] = R8(p + 16);        // 经度 Longitude (deg)
    pos->Pos[2] = R8(p + 24);        // 高度 Height (m)

    return 0;
}

int Decode_bestpos(unsigned char* data, int len, EPOCHOBSDATA* obs)
{
    unsigned char* p = data + 28; // 跳过28字节消息头
    //memset(obs, 0, sizeof(POSRES));

    obs->BLH[0] = R8(p + 8) / 180.0 * pi;        // 纬度 Latitude (rad)
    obs->BLH[1] = R8(p + 16) / 180.0 * pi;        // 经度 Longitude (rad)
    obs->BLH[2] = R8(p + 24);        // 高度 Height (m)
    ElliPara para;
    para.a = 6378137.000;
    para.b = 6356752.314;
    para.e = 0.081819790992;
    BLHToXYZ( &para, obs->BLH, obs->Pos);
    return 0;
}


int DecodeNovOem7Dat(unsigned char buff[], int& len, EPOCHOBSDATA* obs, GPSEPHREC geph[], GPSEPHREC beph[], POSRES* pos)
{
    int i,  MsgLen, MsgId, Status;  // 记录待处理数据的起始位置

    i = 0;
    Status = 0;

    while (1)
    {
        //1. 设置循环变量i = 0，开始查找AA 44 12同步字符
        for (; i < len - 2; i++) {
            if (buff[i] == 0xAA && buff[i + 1] == 0x44 && buff[i + 2] == 0x12) break;
        }

        //    2. 找到同步字符后，获取消息头长度的字符28字节， 若字节数量不足即i + 28 > len，跳出循环（break）至第6步
        if (i + 28 > len) break;

        //    3. 从消息头中解码消息长度MsgLen和消息类型MsgID，获得整条消息buff[i, i + 28 + MsgLen + 4]，若字节数量不足即i + 28 + MsgLen + 4 > len，跳出循环（break） 至第6步
        MsgId = UI2(buff + i + 4);
        MsgLen = UI2(buff + i + 8);
        if (i + 28 + MsgLen + 4 > len) break;

        //    4. CRC检验，若不通过，跳过同步字符3个字节，即i = i + 3，返回到第1步
        if (crc32(buff + i, 28 + MsgLen) != UI4(buff + i + 28 + MsgLen))
        {
            i += 3;  // 跳过当前AA 44 12
            continue;
        }

        //    5. 根据消息ID，调用对应的解码函数，若解码得到观测值，跳出循环至第6步；否则，跳过整条消息，即i = i + 28 + MsgLen + 4，返回第1步
       switch (MsgId)
        {
        case 43:
            DecodeRange(buff + i, MsgLen, obs);
            Status = 1;
            break;
        case 7:
            DecodeGpsEphem(buff + i, MsgLen, geph);
            break;
        case 1696:
            DecodeBdsEphem(buff + i, MsgLen, beph);
            break;
        case 47:
            Decode_psrpos(buff + i, MsgLen, pos);
            break;
        case 42:
            Decode_bestpos(buff + i, MsgLen, obs);
            break;
        default:
            //printf("Unknown message ID: 0x%04X\n", MsgId);
            break;
        }
        i = i + 28 + MsgLen + 4;
        if (Status == 1 && FILEMODE == 1) break;
    }

    //    6. 循环结束后，将不足一条消息的剩余字节，拷贝至buff缓冲区的开始
    //    处，将len设置为剩余字节数量，并返回给主函数。

    memcpy(buff + 0, buff + i, len - i);
    len = len - i;
    return Status;
}

