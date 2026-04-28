#include <iostream>
#include <fstream>
#include <iomanip>
#include "sockets.h"
#include "RTK.h"

using namespace std;
using namespace Eigen;

// 模式选择：文件模式为1，TCP模式为其他

int main()
{
    // 文件模式相关变量
    FILE* FObs_base = nullptr;
    FILE* FObs_rover = nullptr;

    // TCP模式相关变量
    SOCKET sock_base = INVALID_SOCKET;
    SOCKET sock_rover = INVALID_SOCKET;

    // 基站和流动站IP/端口
    const char* BASE_IP = BaseIP;// "47.114.134.129";
    const unsigned short BASE_PORT = BasePort;
    const char* ROVER_IP = RoverIP;
    const unsigned short ROVER_PORT = RoverPort1;

    RAWDAT rawdata;
    PPRESULT PosBas, PosRov;
    RTKEKF rtkekf;

    double XYZ_LQ[3] = { 0.0 }, XYZ_KF[3] = { 0.0 }, X0[3] = { 0.0 },BLH[3] = { 0.0 }, X0_LQ[3] = { 0.0 }, X0_KF[3] = { 0.0 };
    int nIdx_LQ = 1, nIdx_KF = 1; // 平滑索引
    ElliPara para;
    para.a = 6378137.000;
    para.b = 6356752.314;
    para.e = 0.081819790992;
    MatrixXd ENU = MatrixXd::Zero(3, 1);

    ofstream outputFileLQ(LQ_File);
    ofstream outputFileRF(KF_File);
    outputFileLQ << fixed << setprecision(6);
    outputFileRF << fixed << setprecision(6);
    cout << fixed << setprecision(6);

    // ================== 模式初始化 ==================
#if FILEMODE == 1
    // 文件模式初始化
    if ((FObs_base = fopen(ShortFileBaseRT, "rb")) == NULL) {
        printf("Cannot open Base obs file. \n");
        return 0;
    }
    if ((FObs_rover = fopen(ShortFileRoverRT, "rb")) == NULL) {
        printf("Cannot open Rover obs file. \n");
        fclose(FObs_base);
        return 0;
    }
    printf("File mode initialized successfully!\n");
#else
    // TCP模式初始化
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) {
        cerr << "WSAStartup failed!" << endl;
        return -1;
    }

    // 连接基站
    if (!OpenSocket(sock_base, BASE_IP, BASE_PORT)) {
        cerr << "Failed to connect to base station!" << endl;
        WSACleanup();
        return -1;
    }
    cout << "Connected to base station!" << endl;

    // 连接流动站
    if (!OpenSocket(sock_rover, ROVER_IP, ROVER_PORT)) {
        cerr << "Failed to connect to rover!" << endl;
        CloseSocket(sock_base);
        WSACleanup();
        return -1;
    }
    cout << "Connected to rover!" << endl;

    // 设置非阻塞模式
    u_long mode = 1;
    ioctlsocket(sock_base, FIONBIO, &mode);
    ioctlsocket(sock_rover, FIONBIO, &mode);

    // 初始化缓冲区
    memset(&rawdata, 0, sizeof(RAWDAT));
    printf("TCP mode initialized successfully!\n");
#endif

    // ================== 主处理循环 ==================
    while (true) {
#if FILEMODE == 1
        // 文件模式结束检查
        if (feof(FObs_base) || feof(FObs_rover)) {
            cout << "End of files reached." << endl;
            break;
        }
#endif

        // 时间同步
        int syncResult;
#if FILEMODE == 1
        syncResult = TimeSyn(FObs_base, FObs_rover, &rawdata);
#else
        syncResult = TimeSynTCP(sock_base, sock_rover, &rawdata);
#endif

        if (syncResult == 1) {
            // 时间同步成功，处理数据
            DetectOutlier(&rawdata.BasEpk);
            DetectOutlier(&rawdata.RovEpk);

            // 单点定位
            if (SPP(&rawdata.BasEpk, &rawdata, &PosBas) &&
                SPP(&rawdata.RovEpk, &rawdata, &PosRov)) {

                SPV(&rawdata.BasEpk, &PosBas);
                SPV(&rawdata.RovEpk, &PosRov);

                SingleDiff(&rawdata.RovEpk, &rawdata.BasEpk, &rawdata.SdObs);
                DetectCycleSlip(&rawdata.SdObs);
                DetRefSat(&rawdata.BasEpk, &rawdata.RovEpk, &rawdata.SdObs, &rawdata.DDObs);

                //if (RTKFloat(&rawdata, &PosBas, &PosRov)) {
                //    // 处理定位结果
                //    
                //    for (int i = 0; i < 3; i++) XYZ_LQ[i] = rawdata.DDObs.Position[i];

                //    if (rawdata.DDObs.bFixed) {
                //        X0_LQ[0] = (X0_LQ[0] * (nIdx_LQ - 1) + XYZ_LQ[0]) / nIdx_LQ;
                //        X0_LQ[1] = (X0_LQ[1] * (nIdx_LQ - 1) + XYZ_LQ[1]) / nIdx_LQ;
                //        X0_LQ[2] = (X0_LQ[2] * (nIdx_LQ - 1) + XYZ_LQ[2]) / nIdx_LQ;
                //        nIdx_LQ++;
                //    }

                //    XYZToBLH(&para, XYZ_LQ, BLH);
                //    CompEnudPos(rawdata.BasEpk.BLH, X0, rawdata.DDObs.dPos, ENU);

                //    // 输出结果
                //    cout << "LQ:";
                //    cout << rawdata.RovEpk.Time.Week << "\t"
                //        << rawdata.RovEpk.Time.SecOfWeek << "\t"
                //        << BLH[0] * RAD_To_DEG << "\t" << BLH[1] * RAD_To_DEG << "\t" << BLH[2] << "\t" << rawdata.DDObs.bFixed << "\t" << rawdata.DDObs.Ratio << "\t"
                //        << XYZ_LQ[0] << "\t" << XYZ_LQ[1] << "\t" << XYZ_LQ[2] << "\t"
                //        << rawdata.DDObs.dPos[0] << "\t" << rawdata.DDObs.dPos[1] << "\t" << rawdata.DDObs.dPos[2] << "\t"
                //        << ENU(0) << "\t" << ENU(1) << "\t"
                //        << ENU(2) << endl;

                //    outputFileLQ << rawdata.RovEpk.Time.Week << "\t"
                //        << rawdata.RovEpk.Time.SecOfWeek << "\t"
                //        << BLH[0] * RAD_To_DEG << "\t" << BLH[1] * RAD_To_DEG << "\t" << BLH[2] << "\t" << rawdata.DDObs.bFixed << "\t"
                //        << XYZ_LQ[0] << "\t" << XYZ_LQ[1] << "\t" << XYZ_LQ[2] << "\t"
                //        << rawdata.DDObs.dPos[0] << "\t" << rawdata.DDObs.dPos[1] << "\t" << rawdata.DDObs.dPos[2] << "\t"
                //        << ENU(0) << "\t" << ENU(1) << "\t"
                //        << ENU(2) << endl;
                //}
                if (RTKFloatKalman(&rawdata, &PosBas, &PosRov, &rtkekf))
                {
                    // 处理定位结果
               
                    for (int i = 0; i < 3; i++) XYZ_KF[i] = rawdata.DDObs.Position[i];

                    if (rawdata.DDObs.bFixed) {
                        X0_KF[0] = (X0_KF[0] * (nIdx_KF - 1) + XYZ_KF[0]) / nIdx_KF;
                        X0_KF[1] = (X0_KF[1] * (nIdx_KF - 1) + XYZ_KF[1]) / nIdx_KF;
                        X0_KF[2] = (X0_KF[2] * (nIdx_KF - 1) + XYZ_KF[2]) / nIdx_KF;
                        nIdx_KF++;
                    }

                    XYZToBLH(&para, XYZ_KF, BLH);
                    CompEnudPos(rawdata.BasEpk.BLH, X0, rawdata.DDObs.dPos, ENU);

                    // 输出结果
                    cout << "RF:";
                    cout << rawdata.RovEpk.Time.Week << "\t"
                        << rawdata.RovEpk.Time.SecOfWeek << "\t"
                        << BLH[0] * RAD_To_DEG << "\t" << BLH[1] * RAD_To_DEG << "\t" << BLH[2] << "\t" << rawdata.DDObs.bFixed << "\t" << rawdata.DDObs.Ratio << "\t"
                        << X0_KF[0] << "\t" << X0_KF[1] << "\t" << X0_KF[2] << "\t"
                        << rawdata.DDObs.dPos[0] << "\t" << rawdata.DDObs.dPos[1] << "\t" << rawdata.DDObs.dPos[2] << "\t"
                        << ENU(0) << "\t" << ENU(1) << "\t"
                        << ENU(2) << endl;

                    outputFileRF << rawdata.RovEpk.Time.Week << "\t"
                        << rawdata.RovEpk.Time.SecOfWeek << "\t"
                        << BLH[0] * RAD_To_DEG << "\t" << BLH[1] * RAD_To_DEG << "\t" << BLH[2] << "\t" << rawdata.DDObs.bFixed << "\t"
                        << X0_KF[0] << "\t" << X0_KF[1] << "\t" << X0_KF[2] << "\t"
                        << rawdata.DDObs.dPos[0] << "\t" << rawdata.DDObs.dPos[1] << "\t" << rawdata.DDObs.dPos[2] << "\t"
                        << ENU(0) << "\t" << ENU(1) << "\t"
                        << ENU(2) << endl;
                }
            }

            // 重置数据结构
            memset(&rawdata.BasEpk, 0, sizeof(EPOCHOBSDATA));
            memset(&rawdata.RovEpk, 0, sizeof(EPOCHOBSDATA));
            memset(&rawdata.SdObs, 0, sizeof(SDEPOCHOBS));
            rawdata.SdObs.SatNum = 0;
            rawdata.SdObs.SatNumGPS = 0;
            rawdata.SdObs.SatNumBDS = 0;
        }
        else if (syncResult < 0) {
            // 错误处理
#if FILEMODE == 1
            if (syncResult == -1) {
                cout << "File read error or end of file reached." << endl;
                break;
            }
#else
            switch (syncResult) {
            case -2:
                cout << "TCP connection closed." << endl;
                break;
            case -3:
                cout << "TCP receive error." << endl;
                break;
            default:
                cout << "Unknown time sync error: " << syncResult << endl;
            }
            break;
#endif
        }
#if FILEMODE != 1
        else if (syncResult == 0) {
            // TCP模式下时间同步失败但未出错，短暂等待后继续
            Sleep(10);
        }
#endif
    }

    // ================== 清理资源 ==================
#if FILEMODE == 1
    if (FObs_base) fclose(FObs_base);
    if (FObs_rover) fclose(FObs_rover);
#else
    CloseSocket(sock_base);
    CloseSocket(sock_rover);
    WSACleanup();
#endif

    outputFileLQ.close();
    outputFileRF.close();
    cout << "结果已保存" << endl;

    return 0;
}