#pragma once
#define R_WGS84  6378137.0          /* Radius Earth [m]; WGS-84  */
#define F_WGS84  1.0/298.257223563  /* Flattening; WGS-84   */

#define R_CGS2K  6378137.0          /* Radius Earth [m]; CGCS2000  */
#define F_CGS2K  1.0/298.257222101  /* Flattening; CGCS2000   */

#define GPST_BDT_SEC  14         /* GPS时与北斗时的差值[s] */
#define GPST_BDT_WEEK  1356         /* GPS时与北斗时的差值[week] */

#define MAXCHANNUM 36                 //最大通道数量
#define MAXGPSNUM  32                 //GPS卫星最大数量
#define MAXBDSNUM 63                  //BDS卫星最大数量
#define MAXRAWLEN 40960

#define FILEMODE 1               //文件模式为1，TCP模式为其他
#define POLYCRC32   0xEDB88320u /* CRC32 polynomial */
#define C_Light 299792458.0      
#define  FG1_GPS  1575.42E6             /* L1信号频率 */
#define  FG2_GPS  1227.60E6             /* L2信号频率 */
#define  WL1_GPS  (C_Light/FG1_GPS)
#define  WL2_GPS  (C_Light/FG2_GPS)

#define  FG1_BDS  1561.098E6               /* B1信号的基准频率 */
#define  FG3_BDS  1268.520E6               /* B3信号的基准频率 */
#define  WL1_BDS  (C_Light/FG1_BDS)
#define  WL3_BDS  (C_Light/FG3_BDS)       // 波长

#define u_GPS 3.986005E14
#define u_BDS 3.986004418E14
#define Omega_e_dot_GPS 7.2921151467E-5
#define u_BDS 3.986004418E14
#define Omega_e_dot_BDS 7.2921150E-5
#define omega_e  7.2921151467E-5
#define F -4.442807633E-10

const double pi = 3.1415926535898;
const double c_V = 2.99792458E8;
#define RAD_To_DEG    180.0/pi

#define ZeroFileBase     "oem719-202202021500-base.bin"         //零基线基站文件名
#define ZeroFileRover    "oem719-202202021500-rover.bin"        //零基线流动站文件名
#define ShortFileBase    "oem719-202503261140-base.bin"         //短基线基站文件名
#define ShortFileRover   "oem719-202503261140-rover.bin"        //短基线流动站文件名
#define ShortFileBaseRT  "oem719-202510311730-base.bin"         //20h实时数据
#define ShortFileRoverRT "oem719-202510311730-rover.bin"        //20h实时数据


#define BaseIP          "47.114.134.129"                       //基站IP
#define RoverIP         "8.148.22.229"                         //流动站IP
#define BasePort        7190                                 //基站端口
#define RoverPort1      5002                                 //流动站端口1
#define RoverPort2      7002                                 //流动站端口2

#define LQ_File         "testLQ.txt"                        //最小二乘方法写入文件名LqResults.txt
#define KF_File         "testKF.txt"                        //Kalman滤波方法写入文件名KfResults.txt