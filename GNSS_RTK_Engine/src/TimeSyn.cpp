#include<iostream>
#include"RTK.h"



int TimeSyn(FILE* FObs_base, FILE* FObs_rover, RAWDAT* rawdata)
{
    double dt;
    int lenR_base,lenR_rover;        /*已接收数据长度，实际接收长度*/
    static int lenD_base = 0, lenD_rover = 0;  /*剩余长度*/
    static unsigned char buff_base[MAXRAWLEN], buff_rover[MAXRAWLEN];      /*最大缓冲区*/

    // 获取流动站的观测数据，若不成功，继续获取，若文件结束，返回;若成功，得到观测时刻。
    while (!feof(FObs_rover))
    {
        if ((lenR_rover = fread(buff_rover + lenD_rover, sizeof(unsigned char), MAXRAWLEN - lenD_rover, FObs_rover)) < MAXRAWLEN - lenD_rover) return -1;
        lenD_rover += lenR_rover;
        if (DecodeNovOem7Dat(buff_rover, lenD_rover, &rawdata->RovEpk, rawdata->GpsEph, rawdata->BdsEph, &rawdata->bestPos_rover) == 1) break;
    }

    //流动站观测时刻与当前基站观测时刻比较，在限差范围内，文件同步成功，返回1；
    dt = (rawdata->RovEpk.Time.Week - rawdata->BasEpk.Time.Week) * 604800 + rawdata->RovEpk.Time.SecOfWeek - rawdata->BasEpk.Time.SecOfWeek;
    if (fabs(dt) < 0.01) return 1;

    //若不在限差范围内，获取基站观测数据，两站的观测时间求差，在限差范围内，同步成功返回1；
    while (!feof(FObs_base))
    {
        if ((lenR_base = fread(buff_base + lenD_base, sizeof(unsigned char), MAXRAWLEN - lenD_base, FObs_base)) < MAXRAWLEN - lenD_base) return -1;
        lenD_base += lenR_base;
        if (DecodeNovOem7Dat(buff_base, lenD_base, &rawdata->BasEpk, rawdata->GpsEph, rawdata->BdsEph, &rawdata->bestPos_base) == 1)
        {
            dt = (rawdata->RovEpk.Time.Week - rawdata->BasEpk.Time.Week) * 604800 + rawdata->RovEpk.Time.SecOfWeek - rawdata->BasEpk.Time.SecOfWeek;
            if (fabs(dt) < 0.01) return 1;
            //如果不在限差范围内，若基站时间在后，返回0；若基站时间在前，循环获取基站数据，直到成功。
            else if(dt < 0.5) return 0;
        //else;
        }
    }
}

// （TCP模式）时间同步函数
/*返回值=1，时间同步成功
* 返回值=0，基站时间在流动站时间后面
* 返回值=-2，连接关闭
* 返回值=-3，网络连接错误*/
int TimeSynTCP(SOCKET sock_base, SOCKET sock_rover, RAWDAT* rawdata)
{
    double dt;
    int lenR_base, lenR_rover;
    static int lenD_base = 0, lenD_rover = 0;
    static unsigned char buff_base[MAXRAWLEN], buff_rover[MAXRAWLEN];

    // 设置非阻塞模式,如果没有数据可读，立即返回错误
    u_long mode = 1;
    ioctlsocket(sock_rover, FIONBIO, &mode);
    ioctlsocket(sock_base, FIONBIO, &mode);

    
    while (true) {
        lenR_rover = recv(sock_rover, (char*)(buff_rover + lenD_rover),MAXRAWLEN - lenD_rover, 0);

        if (lenR_rover > 0) {
            lenD_rover += lenR_rover;

            if (DecodeNovOem7Dat(buff_rover, lenD_rover, &rawdata->RovEpk, rawdata->GpsEph, rawdata->BdsEph, &rawdata->bestPos_rover) == 1) {
                // 成功解码观测值，跳出循环
                break;
            }
            // 如果是0（其他消息）继续接收更多数据
  
        }
        else if (lenR_rover == 0) {
            return -2;
        }
        //没有数据，但是SOCKET正常
        else if (WSAGetLastError() == WSAEWOULDBLOCK) {
            Sleep(10);
            continue;
        }
        else {
            return -3;
        }
    }

    // 计算时间差
    dt = (rawdata->RovEpk.Time.Week - rawdata->BasEpk.Time.Week) * 604800 +
        rawdata->RovEpk.Time.SecOfWeek - rawdata->BasEpk.Time.SecOfWeek;

    if (fabs(dt) < 10) return 1;

    while (true) {
        lenR_base = recv(sock_base, (char*)(buff_base + lenD_base),MAXRAWLEN - lenD_base, 0);

        if (lenR_base > 0) {
            lenD_base += lenR_base;

            if (DecodeNovOem7Dat(buff_base, lenD_base, &rawdata->BasEpk, rawdata->GpsEph, rawdata->BdsEph, &rawdata->bestPos_base) == 1) {
                // 成功解码基站观测值，进行时间比较
                dt = (rawdata->RovEpk.Time.Week - rawdata->BasEpk.Time.Week) * 604800 +
                    rawdata->RovEpk.Time.SecOfWeek - rawdata->BasEpk.Time.SecOfWeek;

                if (fabs(dt) < 10) return 1;
                if (dt > 10) continue;  // 时间差太大，继续接收更新的数据
                return 0;  // 基站在流动站后
            }
            
        }
        else if (lenR_base == 0) {
            return -2;
        }
        else if (WSAGetLastError() == WSAEWOULDBLOCK) {
            Sleep(10);
            continue;
        }
        else {
            return -3;
        }
    }
}