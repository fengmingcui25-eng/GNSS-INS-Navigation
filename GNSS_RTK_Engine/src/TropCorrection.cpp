#include"RTK.h"

double Hopfield(const double H, const double Elev)//高度角输入为角度
{
	if (H < -500.0 || H > 10000.0) 
	{
		return 0.0;  // 返回 0 延迟（等效于暂不改正）
	}
	double deltaTrop,Kd, Kw, 
		hd, hw, e, T, p, RH;
	T = T0 - 0.0065 * (H - H0);
	p = p0 * pow((1 - 0.0000226 * (H - H0)), 5.225);
	RH = RH0 * exp( -0.0006396 * (H - H0));
	e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);
	hd = 40136 + 148.72 * (T0 - 273.16);
	hw = 11000.0;
	Kd = 155.2 * 1E-7 * p / T * (hd - H);
	Kw = 155.2 * 1E-7 * 4810 / T / T * e * (hw - H);
	deltaTrop = Kd / sin(sqrt(Elev * Elev + 6.25) * pi / 180.0) + Kw / sin(sqrt(Elev * Elev + 2.25) * pi / 180.0);
    return deltaTrop;
	
}

