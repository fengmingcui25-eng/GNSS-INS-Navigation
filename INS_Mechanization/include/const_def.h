#pragma once


#define PI  3.141592653589793238
#define RAD_TO_DEG  (180.0/PI)
#define DEG_TO_RAD  (PI/180.0)
#define MAXRAWLEN  40960/*最大缓冲区*/

#define acc_scale  1.5258789063E-06/*转换因子*//*星网宇达*/
#define gyo_scale  1.0850694444E-07/*转换因子*//*星网宇达*/
#define SampleRate  100//采样率/**************************************************
//#define SampleRate_EX 200//示例数据采样率
#define dt  (1.0/SampleRate)
//#define dt_EX  (1.0/SampleRate_EX)
#define Latitude  30.531651244//纬度，单位：度
#define w_e  7.2921151467E-5//地球自转角速度,单位：rad/s

#define GRAVITY  9.7803267715          // 重力加速度 (m/s^2)
//#define latitude_local  30.531651244/*当地纬度*/
#define WorkMode  1               /*工作状态：0-->示例数据；1-->采集数据*/

// 物理常数
#define R_WGS84	  6378137.0               // 长半轴
#define F_WGS84	  1.0/298.257223563		  // 扁率
#define E2_WGS84    0.0066943800229         // 第一偏心率平方
