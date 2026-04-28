#include"RTK.h"
using namespace std;



//通用时转儒略日
void CommonTimeToMjdTime(const COMMONTIME* ct, MJDTIME* mjd)
{
	short y;
	unsigned short m;
	if (ct->Month <= 2)
	{
		y = ct->Year - 1;
		m = ct->Month + 12;
	}
	else if (ct->Month > 2)
	{
		y = ct->Year;
		m = ct->Month;
	}

	double JD = int(365.25 * y) + int(30.6001 * (m + 1)) + ct->Day + (ct->Hour + ct->Minute / 60 + ct->Second / 3600.0) / 24.0 + 1720981.5;
	double MJD = JD - 2400000.5;
	mjd->Days = int(MJD);
	mjd->FracDay = (ct->Hour + ct->Minute / 60.0 + ct->Second / 3600.0) / 24.0;
	
}

//儒略日转通用时
void MjdTimeToCommonTime(const MJDTIME* mjd, COMMONTIME* ct)
{
	int a = int(mjd->Days + 2400000.5 + mjd->FracDay + 0.5);
	int b = a + 1537;
	int c = int((b - 122.1) / 365.25);
	int d = int(365.25 * c);
	int e = int((b - d) / 30.6001);
	ct->Day = b - d - int(30.6001 * e);
	ct->Month = e - 1 - 12 * int(e / 14);
	ct->Year = c - 4715 - int((7 + ct->Month) / 10);
	ct->Hour = int(mjd->FracDay * 24.0);
	ct->Minute = int((mjd->FracDay * 24.0 - ct->Hour) * 60);
	ct->Second = fmod((mjd->FracDay) * 24 * 3600.0, 60.0);
}


//儒略日转Gps时
void MjdTimeToGPSTime(const MJDTIME* mjd, GPSTIME* Gpst)
{
	Gpst->Week = int((mjd->Days + mjd->FracDay - 44244) / 7);
	Gpst->SecOfWeek = (mjd->Days + mjd->FracDay - 44244 - Gpst->Week * 7) * 86400;
}

//GPS时转儒略日
void GPSTimeToMjdTime(const GPSTIME* Gpst, MJDTIME* mjd)
{
	double MJD = 44244.0 + Gpst->Week * 7 + Gpst->SecOfWeek / 86400;
	mjd->Days = int(MJD);
	mjd->FracDay = Gpst->SecOfWeek / 86400;
}

//通用时转GPS时
void CommonTimeToGPSTime(const COMMONTIME* ct, GPSTIME* Gpst)
{
	MJDTIME mjd;
	CommonTimeToMjdTime( ct,  &mjd);
	MjdTimeToGPSTime(&mjd, Gpst);
}

//GPS时转通用时
void GPSTimeToCommonTime(const GPSTIME* Gpst, COMMONTIME* ct)
{
	MJDTIME mjd;
	GPSTimeToMjdTime( Gpst,  &mjd);
	MjdTimeToCommonTime(&mjd, ct);
}

/*int main()
{
	//测试
	/*MJDTIME mjd;
	COMMONTIME ct;
	mjd.Days = 60775; 
	mjd.FracDay = 0.50960648148;
	MjdTimeToCommonTime(&mjd, &ct);
	cout << ct.Year << endl;
	cout << ct.Month << endl;
	cout << ct.Day << endl;
	cout << ct.Hour << endl;
	cout << ct.Minute << endl;
	cout << ct.Second << endl;*/
	/*GEOCOOR Blh;
	Blh.BLH[0] = 0.1;//30.527819302;
	Blh.BLH[1] = 0.01;// 114.356564227;
	Blh.BLH[2] = 80.804;
	Eigen::MatrixXd R(3, 3);
	BLHToNEUMat(Blh.BLH, R);
	cout << R << endl;
	return 0;
}*/