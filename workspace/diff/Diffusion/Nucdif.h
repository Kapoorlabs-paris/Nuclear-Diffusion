
#ifndef Nucdif_h
#define Nucdif_h Nucdif_h
#include<iostream>
#include<complex>
#include<cmath>
#include<fstream>
#include<string>


using namespace std;

class density;
class fluid;
class grid;
class hamop;

double vecpot_x(double time);
double vecpot_y(double time);
double vecpot_z(double time);
double scalarpotx(double x, double y, double z, double time);
double scalarpoty(double x, double y, double z, double time);
double scalarpotz(double x, double y, double z, double time);
double interactionpotxy(double x, double y, double z, double time);
double imagpotx(long xindex, long yindex, long zindex, double time, grid g);
double imagpoty(long xindex, long yindex, long zindex, double time, grid g);
double mitosis(double time);




#endif // qprop_h
