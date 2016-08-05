#ifndef hamop_h
#define hamop_h hamop_h
#include<complex>
#include<iostream>
#include"grid.h"

using namespace std;

class density;
class fluid;

class hamop
{
 public:
  hamop(grid g,
	double (*fpx)(double),
	double (*fpy)(double),
	double (*fpz)(double),
	double (*fpsx)(double, double, double, double),
	double (*fpsy)(double, double, double, double),
	double (*fpsz)(double, double, double, double),
	double (*fpixy)(double, double, double, double),
	double (*fpimx)(long, long, long, double, grid),
	double (*fpimy)(long, long, long, double, grid),
	double (*mit)(double)
       );
  double vecpot_x(double time);
  double vecpot_y(double time);
  double vecpot_z(double time);
  double scalarpot(double x, double y, double z, double time);
  double scalarpotx(double x, double y, double z, double time);
  double scalarpoty(double x, double y, double z, double time);
  double scalarpotz(double x, double y, double z, double time);
  double interactionpotxy(double x, double y, double z, double time);
  double imagpot(long xindex, long yindex, long zindex, double time, grid g);
  double imagpotx(long xindex, long yindex, long zindex, double time, grid g);
  double imagpoty(long xindex, long yindex, long zindex, double time, grid g);
  double mitosis(double time);


 private:
  double delta_x, delta_y, delta_z;
  double (*hamopvecpotx)(double);
  double (*hamopvecpoty)(double);
  double (*hamopvecpotz)(double);
  double (*hamopscalarpotx)(double, double, double, double);
  double (*hamopscalarpoty)(double, double, double, double);
  double (*hamopscalarpotz)(double, double, double, double);
  double (*hamopinteractionpotxy)(double, double, double, double);
  double (*hamopimagpotx)(long, long, long, double, grid);
  double (*hamopimagpoty)(long, long, long, double, grid);
  double (*hamopmitosis)(double);


};



#endif // hamop_h




