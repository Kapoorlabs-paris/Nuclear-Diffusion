#include"hamop.h"
#include"density.h"
#include"fluid.h"


hamop::hamop(grid g,
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

	    )
{
  delta_x=g.delt_x();
  delta_y=g.delt_y();
  delta_z=g.delt_z();

  hamopvecpotx=fpx;
  hamopvecpoty=fpy;
  hamopvecpotz=fpz;
  hamopscalarpotx=fpsx;
  hamopscalarpoty=fpsy;
  hamopscalarpotz=fpsz;
  hamopinteractionpotxy=fpixy;
  hamopimagpotx=fpimx;
  hamopimagpoty=fpimy;
  hamopmitosis = mit;

}

double hamop::vecpot_x(double time )
{
  return hamopvecpotx(time);

}

double hamop::vecpot_y(double time )
{
  return hamopvecpoty(time);

}

double hamop::vecpot_z(double time )
{
  return hamopvecpotz(time);

}

double hamop::mitosis(double time)
{
	return hamopmitosis(time);
}

double hamop::scalarpot(double x, double y, double z, double time )
{
  return hamopscalarpotx(x,y,z,time)+hamopscalarpoty(x,y,z,time)+hamopscalarpotz(x,y,z,time);
}

double hamop::scalarpotx(double x, double y, double z, double time )
{
  return hamopscalarpotx(x,y,z,time);
}
double hamop::scalarpoty(double x, double y, double z, double time )
{
  return hamopscalarpoty(x,y,z,time);
}
double hamop::scalarpotz(double x, double y, double z, double time )
{
  return hamopscalarpotz(x,y,z,time);
}

double hamop::interactionpotxy(double x, double y, double z, double time )
{
  return hamopinteractionpotxy(x,y,z,time);
}

double hamop::imagpot(long xindex, long yindex, long zindex, double time, grid g)
{
  return hamopimagpotx(xindex,yindex,zindex,time,g)+hamopimagpoty(xindex,yindex,zindex,time,g);
}

double hamop::imagpotx(long xindex, long yindex, long zindex, double time, grid g)
{
  return hamopimagpotx(xindex,yindex,zindex,time,g);
}

double hamop::imagpoty(long xindex, long yindex, long zindex, double time, grid g)
{
  return hamopimagpoty(xindex,yindex,zindex,time,g);
}








