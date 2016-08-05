/*
 * density.cpp
 *
 *  Created on: Apr 28, 2016
 *      Author: varunkapoor
 */

#include "density.h"
#include"hamop.h"
#include"fluid.h"
#include"grid.h"
#include<cmath>
#define THRESH 1e-6
#define OOS 1.0/6.0
#define TOT 2.0/3.0
#define FOT 5.0/3.0
#define SQRTT sqrt(2.0)

#define dcomplex complex<double>

void density::nullify() {
	for (long i = 0; i < dens_dim; i++)
		start[i] = 0.0 * start[i]; //complex(0.0,0.0);
}

// ------------------ for other overloaded version see below ------------------
void density::init(grid g, int inittype, double width, int k, double time,
		FILE* filename, int output_of_interest) {
	long xindex, yindex, zindex, index_i, index, i;
	double x, y, z;
	double offs = 0.0;
	double realpart, imagpart;
	long ctrl_id;

	if (inittype == 99) {
		for (i = 0; i < output_of_interest; i++) {
			cout << "Scanning output no. " << i << "... - ignoring ..." << endl;
			for (zindex = 0; zindex < g.ngps_z(); zindex++) {
				for (xindex = 0; xindex < g.ngps_x(); xindex++) {
					for (yindex = 0; yindex < g.ngps_y(); yindex++) {
						ctrl_id = fscanf(filename, "%lf %lf", &realpart,
								&imagpart);
					};
				};
			};
		};
		cout << "Now I'm storing!" << endl;
		for (zindex = 0; zindex < g.ngps_z(); zindex++) {
			for (xindex = 0; xindex < g.ngps_x(); xindex++) {
				for (yindex = 0; yindex < g.ngps_y(); yindex++) {
					index = g.index(xindex, yindex, zindex);
					ctrl_id = fscanf(filename, "%lf %lf", &realpart, &imagpart);
					//		  printf("%e %e\n",realpart,imagpart);
					start[index] = dcomplex(realpart,imagpart);
				};
			};
		};
	};
}

void density::init(grid g, int inittype, double width, int k, double time) {
	long xindex, yindex, zindex, index_i;
	double x, y, z;
	dcomplex st;
	double offs = 0.0;
	int center = 0;

	st = width + 0.5 * time / width;

	srand((long) (width));
	density tmp(g.ngps_x() * g.ngps_y() * g.ngps_z());
	tmp.nullify();
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		x = g.x(xindex);
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			y = g.y(yindex);
			for (zindex = 0; zindex < g.ngps_z(); zindex++) {
				z = g.z(zindex);
				index_i = g.index(xindex, yindex, zindex);
				switch (inittype) {
				case 0:

					start[index_i] = exp(-0.5 * (x) * (x) * width)
							* exp(-0.5 * (y) * (y) * width)
							+ exp(
									-0.5 * (x - g.ngps_x() / 10)
											* (x - g.ngps_x() / 10) * width)
									* exp(
											-0.5 * (y - g.ngps_y() / 10)
													* (y - g.ngps_y() / 10)
													* width)
							+ exp(
									-0.5 * (x + g.ngps_x() / 10)
											* (x + g.ngps_x() / 10) * width)
									* exp(
											-0.5 * (y + g.ngps_y() / 10)
													* (y + g.ngps_y() / 10)
													* width)
							+ exp(
									-0.5 * (x - g.ngps_x() / 10)
											* (x - g.ngps_x() / 10) * width)
									* exp(
											-0.5 * (y + g.ngps_y() / 10)
													* (y + g.ngps_y() / 10)
													* width)
							+ exp(
									-0.5 * (x + g.ngps_x() / 10)
											* (x + g.ngps_x() / 10) * width)
									* exp(
											-0.5 * (y - g.ngps_y() / 10)
													* (y - g.ngps_y() / 10)
													* width);

					break;
				case 1:
					start[index_i] = exp(-0.5 * x * x * width)
							* exp(-0.5 * y * y * width);
					break;
				case 2:
					start[index_i] = exp(-0.5 * (x) * (x) * width)
							* exp(-0.5 * (y) * (y) * width)
							+ exp(
									-0.5 * (x - g.ngps_x() / 10)
											* (x - g.ngps_x() / 10) * width)
									* exp(
											-0.5 * (y - g.ngps_y() / 10)
													* (y - g.ngps_y() / 10)
													* width)
							+ exp(
									-0.5 * (x + g.ngps_x() / 10)
											* (x + g.ngps_x() / 10) * width)
									* exp(
											-0.5 * (y + g.ngps_y() / 10)
													* (y + g.ngps_y() / 10)
													* width);
					break;
				case 3:
					start[index_i] = sqrt(x*x+y*y)*exp(-0.5 * x * x * width)
					* exp(-0.5 * y * y * width) ;
					break;
				case 4:
					start[index_i] = 1.0;
					break;

				case 5:
					start[index_i] = 0;
					break;
				case 6:
					start[index_i] = cos(3.1415926 / width * x)
							* cos(3.1415926 / width * y);
					break;
				case 7:
					start[index_i] = cos(3.1415926 / width * x)
							* sin(2.0 * 3.1415926 / width * y);
					break;
				case 8:
					start[index_i] = exp(
							-width
									* ((sqrt(x * x + y * y) - 2.6)
											* (sqrt(x * x + y * y) - 2.6)));
					//start[index_i]=exp(-width*((x-6.71)*(x-6.71)+y*y)) + exp(-width*((x+6.71)*(x+6.71)+y*y));
					break;
				case 9:
					start[index_i] = dcomplex(0.0,0.0);
				};
			};
		};
	};

}

density &density::operator=(const density &v) {
	if (this != &v) {
		delete[] start;
		dens_dim = v.dens_dim;
		start = new dcomplex[dens_dim];
		for (long i = 0; i < dens_dim; i++)
			start[i] = v.start[i];
	}
	return *this;
}

density &density::operator=(const fluid &v) {
	delete[] start;
	dens_dim = v.dens_size();
	start = new dcomplex[dens_dim];
	for (long i = 0; i < dens_dim; i++)
		start[i] = v[i];

	return *this;
}

density operator +(const density &v, const density &w) {

	density temp(v.dens_size());
	for (long i = 0; i < v.dens_size(); i++) {
		temp[i] = v[i] + w[i];
	};

	return temp;

}

density operator -(const density &v, const density &w) {

	density temp(v.dens_size());
	for (long i = 0; i < v.dens_size(); i++) {
		temp[i] = v[i] - w[i];
	};

	return temp;

}

density operator +(const density &v, const fluid &w) {
	density temp(v.dens_size());
	for (long i = 0; i < v.dens_size(); i++) {
		temp[i] = v[i] + w[i];
	};
	return temp;
}

density operator +(const fluid &w, const density &v) {
	density temp(v.dens_size());
	for (long i = 0; i < v.dens_size(); i++) {
		temp[i] = v[i] + w[i];
	};
	return temp;

}

density& density::operator *=(double z) {
	for (long i = 0; i < dens_dim; i++)
		start[i] = start[i] * z;
	return *this;
}

density& density::operator *=(dcomplex z)
{
	for(long i=0; i<dens_dim; i++)
	start[i]=start[i]*z;
	return *this;
}

density operator *(double z, const density &v) {
	density temp = v;
	return temp *= z;
}

density operator *(dcomplex z, const density &v)
{
	density temp=v;
	return temp *= z;
}

density operator *(const density &v, double z) {
	density temp = v;
	return temp *= z;
}

density operator *(const density &v, dcomplex z)
{
	density temp=v;
	return temp *= z;
}

dcomplex operator * (const density &v, const density &w )
{
	dcomplex result(0.0,0.0);
	for(long i=0; i<v.dens_size(); i++)
	{
		result+=conj(v[i])*w[i];
	};
	return result;
}

density operator /(const density &v, const density &w) {

	density temp(v.dens_size());
	for (long i = 0; i < v.dens_size(); i++) {
		temp[i] = v[i] / w[i];
	};

	return temp;

}

density operator /(const density &v, double z) {

	density temp(v.dens_size());
	for (long i = 0; i < v.dens_size(); i++) {
		temp[i] = v[i] / z;
	};

	return temp;

}

ostream& operator<<(ostream& os, const density& v) {
	for (long i = 0; i < v.dens_size(); i++) {
		os << real(v[i]) << " " << imag(v[i]) << endl;
	}
	return os;
}

istream& operator>>(istream& is, density& v) {
	double tmpre, tmpim;
	for (long i = 0; i < v.dens_size(); i++) {
		is >> tmpre >> tmpim;
		v[i] = dcomplex(tmpre,tmpim);
	}
	return is;
}

void density::regrid(grid g, grid g_small, const density &v) {
	long xindex, yindex, zindex, index, index_small, xshift, yshift;

	xshift = g.offs_x() - g_small.offs_x();
	yshift = g.offs_y() - g_small.offs_y();

	if ((xshift + g_small.ngps_x() <= g.ngps_x())
			&& (yshift + g_small.ngps_y() <= g.ngps_y())) {

		for (xindex = 0; xindex < g_small.ngps_x(); xindex++) {
			for (yindex = 0; yindex < g_small.ngps_y(); yindex++) {
				for (zindex = 0; zindex < g_small.ngps_z(); zindex++) {
					index = g.index(xindex + xshift, yindex + yshift, zindex);
					index_small = g_small.index(xindex, yindex, zindex);
					start[index] = v[index_small];
				};
			};
		};

		cout << "regridding in density sucessful" << endl;

	} else {
		cout << "Ooops! There's a problem in regridding data! " << endl;
	};

}
;

void density::propagator(double timestep, double time, double diffconst, grid g,
		hamop hamil, const density &staticpot_x, const density &staticpot_y,
		const density &staticpot_xy) {

	switch (g.dimens()) {

	case 2:

		do_cn_step_xy_muller(timestep, time, diffconst, g, hamil, staticpot_x,
				staticpot_y, staticpot_xy,0, 0);
		break;

	};

}

void density::do_cn_step_xy_muller(double timestep, double time,
		double diffconst, grid g, hamop hamil, const density &staticpot_x,
		const density &staticpot_y, const density &staticpot_xy,
		 long yindex, long zindex) {
	long xindex;
	long index, index_xp, index_xm, index_yp, index_ym;
	double x, y, z, halfvecpotvecpot, vecpot;
	density rhsone(g.ngps_x() * g.ngps_y() * g.ngps_z());
	density rhsone_x(g.ngps_x());
	density rhstwo_x(g.ngps_x());
	density rhsone_y(g.ngps_y());
	density rhstwo_y(g.ngps_y());
	double oneoverhsquare;
	density ax(g.ngps_x());
	density bx(g.ngps_x());
	density cx(g.ngps_x());
	density ay(g.ngps_y());
	density by(g.ngps_y());
	density cy(g.ngps_y());
	density timestepstaticpot_xy(g.ngps_x() * g.ngps_y() * g.ngps_z());
	density timestephalfoversixstaticpot_x(g.ngps_x());
	density fivetimestephalfoverthreestaticpot_x(g.ngps_x());
	density timestephalfoversixstaticpot_y(g.ngps_y());
	density fivetimestephalfoverthreestaticpot_y(g.ngps_y());
	density timestepnonlinearpot(g.ngps_x() * g.ngps_y());
	dcomplex aa,bb,cc,timestephalf,timestephalfoversix,fivetimestephalfoverthree;
	double aaa, bbb, ccc;
	double b_upperleft, b_lowerright;
	double vecpotwithprefactor;
	double lambda = sqrt(3.0) - 2.0;
	double llambda = -sqrt(3.0) + 2.0;

	z = g.z(zindex);
	timestephalf = timestep / 2;
	timestephalfoversix = timestephalf / 6.0;
	fivetimestephalfoverthree = 5.0 * timestephalf / 3.0;

	timestephalfoversixstaticpot_x = timestephalfoversix * staticpot_x;
	fivetimestephalfoverthreestaticpot_x = fivetimestephalfoverthree
			* staticpot_x;
	timestephalfoversixstaticpot_y = timestephalfoversix * staticpot_y;
	fivetimestephalfoverthreestaticpot_y = fivetimestephalfoverthree
			* staticpot_y;
	timestepstaticpot_xy = timestep * staticpot_xy;


	// ============================= S_x ========================================

	vecpotwithprefactor =
			(timestep / (8.0 * g.delt_x()) * hamil.vecpot_x(time));

	aaa = OOS - vecpotwithprefactor;
	ccc = OOS + vecpotwithprefactor;
	bbb = TOT;

	// Calculate the rhs vector S_-x *this

	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		xindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_xp = g.index(xindex + 1, yindex, zindex);

		rhsone[index] = aaa * start[index_xp]
				+ ((4.0 + lambda) / 6.0 - lambda * vecpotwithprefactor)
						* start[index];

		for (xindex = 1; xindex < g.ngps_x() - 1; xindex++) {
			index = g.index(xindex, yindex, zindex);
			index_xp = g.index(xindex + 1, yindex, zindex);
			index_xm = g.index(xindex - 1, yindex, zindex);

			rhsone[index] = aaa * start[index_xp] + bbb * start[index]
					+ ccc * start[index_xm];
		};

		xindex = g.ngps_x() - 1;
		index = g.index(xindex, yindex, zindex);
		index_xm = g.index(xindex - 1, yindex, zindex);

		rhsone[index] = ccc * start[index_xm]
				+ ((4.0 + lambda) / 6.0 - llambda * vecpotwithprefactor)
						* start[index];
	}


	// The matrix  S_+x
	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_x[xindex] = rhsone[index];
		};

		aaa = OOS + vecpotwithprefactor;
		ccc = (OOS - vecpotwithprefactor);
		bbb = TOT;

		b_upperleft = ((4.0 + lambda) / 6.0 + lambda * vecpotwithprefactor);
		b_lowerright = ((4.0 + lambda) / 6.0 + llambda * vecpotwithprefactor);

		rhstwo_x.solve_toep(aaa, bbb, b_upperleft, b_lowerright, ccc, rhsone_x,
				g.ngps_x());

		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_x[xindex];
		};
	};

	// ============================= S_y ========================================

	vecpotwithprefactor =
			(timestep / (8.0 * g.delt_y()) * hamil.vecpot_y(time));

	aaa = OOS - vecpotwithprefactor;
	ccc = OOS + vecpotwithprefactor;
	bbb = TOT;

	// Calculate the rhs vector S_-y *this

	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		yindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_yp = g.index(xindex, yindex + 1, zindex);

		rhsone[index] = aaa * start[index_yp]
				+ ((4.0 + lambda) / 6.0 - lambda * vecpotwithprefactor)
						* start[index];

		for (yindex = 1; yindex < g.ngps_y() - 1; yindex++) {
			index = g.index(xindex, yindex, zindex);
			index_yp = g.index(xindex, yindex + 1, zindex);
			index_ym = g.index(xindex, yindex - 1, zindex);

			rhsone[index] = aaa * start[index_yp] + bbb * start[index]
					+ ccc * start[index_ym];
		};

		yindex = g.ngps_y() - 1;
		index = g.index(xindex, yindex, zindex);
		index_ym = g.index(xindex, yindex - 1, zindex);

		rhsone[index] = ccc * start[index_ym]
				+ ((4.0 + lambda) / 6.0 - llambda * vecpotwithprefactor)
						* start[index];
	}

	// The matrix  S_+y
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_y[yindex] = rhsone[index];
		};

		aaa = OOS + vecpotwithprefactor;
		ccc = (OOS - vecpotwithprefactor);
		bbb = TOT;

		b_upperleft = ((4.0 + lambda) / 6.0 + lambda * vecpotwithprefactor);
		b_lowerright = ((4.0 + lambda) / 6.0 + llambda * vecpotwithprefactor);

		rhstwo_y.solve_toep(aaa, bbb, b_upperleft, b_lowerright, ccc, rhsone_y,
				g.ngps_y());

		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_y[yindex];
		};
	};

	// ==================================== W_x =================================

	oneoverhsquare = 1.0 / (g.delt_x() * g.delt_x());
	aa = -OOS - 0.5 * timestephalf * oneoverhsquare * diffconst;
	cc = aa;
	bb = -FOT + 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

	// ---------- W_-x
	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		xindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_xp = g.index(xindex + 1, yindex, zindex);
		rhsone[index] = (aa + 0.5 * timestephalfoversixstaticpot_x[xindex + 1])
				* start[index_xp]
				+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_x[xindex])
						* start[index];

		for (xindex = 1; xindex < g.ngps_x() - 1; xindex++) {
			index = g.index(xindex, yindex, zindex);
			index_xp = g.index(xindex + 1, yindex, zindex);
			index_xm = g.index(xindex - 1, yindex, zindex);
			rhsone[index] = (aa
					+ 0.5 * timestephalfoversixstaticpot_x[xindex + 1])
					* start[index_xp]
					+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_x[xindex])
							* start[index]
					+ (cc + 0.5 * timestephalfoversixstaticpot_x[xindex - 1])
							* start[index_xm];
		};

		xindex = g.ngps_x() - 1;
		index = g.index(xindex, yindex, zindex);
		index_xm = g.index(xindex - 1, yindex, zindex);
		rhsone[index] =
				(bb + 0.5 * fivetimestephalfoverthreestaticpot_x[xindex])
						* start[index]
						+ (cc + 0.5 * timestephalfoversixstaticpot_x[xindex - 1])
								* start[index_xm];

	};

	// --------------- W_+x
	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		y = g.y(yindex);
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_x[xindex] = rhsone[index];
		};
		aa = -OOS + 0.5 * timestephalf * oneoverhsquare * diffconst;
		cc = aa;
		bb = -FOT - 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

		xindex = 0;
		x = g.x(xindex);
		index = g.index(xindex, yindex, zindex);
		index_xp = g.index(xindex + 1, yindex, zindex);
		ax[xindex] = aa - 0.5 * timestephalfoversixstaticpot_x[xindex + 1];
		cx[xindex] = 1.0; // not used
		bx[xindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_x[xindex];

		for (xindex = 1; xindex < g.ngps_x() - 1; xindex++) {
			index = g.index(xindex, yindex, zindex);
			index_xp = g.index(xindex + 1, yindex, zindex);
			index_xm = g.index(xindex - 1, yindex, zindex);

			ax[xindex] = aa - 0.5 * timestephalfoversixstaticpot_x[xindex + 1];
			cx[xindex] = cc - 0.5 * timestephalfoversixstaticpot_x[xindex - 1];
			bx[xindex] = bb
					- 0.5 * fivetimestephalfoverthreestaticpot_x[xindex];

		};

		xindex = g.ngps_x() - 1;
		x = g.x(xindex);
		index = g.index(xindex, yindex, zindex);
		index_xm = g.index(xindex - 1, yindex, zindex);
		ax[xindex] = aa
				- 0.5 * timestephalfoversix
						* hamil.scalarpotx(x + g.delt_x(), y, z, time);
		cx[xindex] = cc - 0.5 * timestephalfoversixstaticpot_x[xindex - 1];
		bx[xindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_x[xindex];

		rhstwo_x.solve(ax, bx, cx, rhsone_x, g.ngps_x());

		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_x[xindex];
		};

	};

	// ==================================== W_y =================================

	oneoverhsquare = 1.0 / (g.delt_y() * g.delt_y());
	aa = -OOS - 0.5 * timestephalf * oneoverhsquare * diffconst;
	cc = aa;
	bb = -FOT + 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

	// ---------- W_-y
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		yindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_yp = g.index(xindex, yindex + 1, zindex);
		rhsone[index] = (aa + 0.5 * timestephalfoversixstaticpot_y[yindex + 1])
				* start[index_yp]
				+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_y[yindex])
						* start[index];

		for (yindex = 1; yindex < g.ngps_y() - 1; yindex++) {
			index = g.index(xindex, yindex, zindex);
			index_yp = g.index(xindex, yindex + 1, zindex);
			index_ym = g.index(xindex, yindex - 1, zindex);
			rhsone[index] = (aa
					+ 0.5 * timestephalfoversixstaticpot_y[yindex + 1])
					* start[index_yp]
					+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_y[yindex])
							* start[index]
					+ (cc + 0.5 * timestephalfoversixstaticpot_y[yindex - 1])
							* start[index_ym];
		};

		yindex = g.ngps_y() - 1;
		index = g.index(xindex, yindex, zindex);
		index_ym = g.index(xindex, yindex - 1, zindex);
		rhsone[index] =
				(bb + 0.5 * fivetimestephalfoverthreestaticpot_y[yindex])
						* start[index]
						+ (cc + 0.5 * timestephalfoversixstaticpot_y[yindex - 1])
								* start[index_ym];

	};

	// --------------- W_+y
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		x = g.x(xindex);
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_y[yindex] = rhsone[index];
		};
		aa = -OOS + 0.5 * timestephalf * oneoverhsquare * diffconst;
		cc = aa;
		bb = -FOT - 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

		yindex = 0;
		y = g.y(yindex);
		index = g.index(xindex, yindex, zindex);
		index_yp = g.index(xindex, yindex + 1, zindex);
		ay[yindex] = aa - 0.5 * timestephalfoversixstaticpot_y[yindex + 1];
		cy[yindex] = 1.0; // not used
		by[yindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_y[yindex];

		for (yindex = 1; yindex < g.ngps_y() - 1; yindex++) {
			index = g.index(xindex, yindex, zindex);
			index_yp = g.index(xindex, yindex + 1, zindex);
			index_ym = g.index(xindex, yindex - 1, zindex);

			ay[yindex] = aa - 0.5 * timestephalfoversixstaticpot_y[yindex + 1];
			cy[yindex] = cc - 0.5 * timestephalfoversixstaticpot_y[yindex - 1];
			by[yindex] = bb
					- 0.5 * fivetimestephalfoverthreestaticpot_y[yindex];

		};

		yindex = g.ngps_y() - 1;
		y = g.y(yindex);
		index = g.index(xindex, yindex, zindex);
		index_ym = g.index(xindex, yindex - 1, zindex);
		ay[yindex] = aa
				- 0.5 * timestephalfoversix
						* hamil.scalarpoty(x, y + g.delt_y(), z, time);
		cy[yindex] = cc - 0.5 * timestephalfoversixstaticpot_y[yindex - 1];
		by[yindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_y[yindex];

		rhstwo_y.solve(ay, by, cy, rhsone_y, g.ngps_y());

		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_y[yindex];
		};

	};

	// ============================ exp(V(xy)) =================================

	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = exp(-timestepstaticpot_xy[index]) * start[index];
		};
	};


	// ==================================== W_x =================================

	oneoverhsquare = 1.0 / (g.delt_x() * g.delt_x());
	aa = -OOS - 0.5 * timestephalf * oneoverhsquare * diffconst;
	cc = aa;
	bb = -FOT + 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

	// ---------- W_-x
	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		xindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_xp = g.index(xindex + 1, yindex, zindex);
		rhsone[index] = (aa + 0.5 * timestephalfoversixstaticpot_x[xindex + 1])
				* start[index_xp]
				+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_x[xindex])
						* start[index];

		for (xindex = 1; xindex < g.ngps_x() - 1; xindex++) {
			index = g.index(xindex, yindex, zindex);
			index_xp = g.index(xindex + 1, yindex, zindex);
			index_xm = g.index(xindex - 1, yindex, zindex);
			rhsone[index] = (aa
					+ 0.5 * timestephalfoversixstaticpot_x[xindex + 1])
					* start[index_xp]
					+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_x[xindex])
							* start[index]
					+ (cc + 0.5 * timestephalfoversixstaticpot_x[xindex - 1])
							* start[index_xm];
		};

		xindex = g.ngps_x() - 1;
		index = g.index(xindex, yindex, zindex);
		index_xm = g.index(xindex - 1, yindex, zindex);
		rhsone[index] =
				(bb + 0.5 * fivetimestephalfoverthreestaticpot_x[xindex])
						* start[index]
						+ (cc + 0.5 * timestephalfoversixstaticpot_x[xindex - 1])
								* start[index_xm];

	};

	// --------------- W_+x
	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		y = g.y(yindex);
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_x[xindex] = rhsone[index];
		};
		aa = -OOS + 0.5 * timestephalf * oneoverhsquare * diffconst;
		cc = aa;
		bb = -FOT - 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

		xindex = 0;
		x = g.x(xindex);
		index = g.index(xindex, yindex, zindex);
		index_xp = g.index(xindex + 1, yindex, zindex);
		ax[xindex] = aa - 0.5 * timestephalfoversixstaticpot_x[xindex + 1];
		cx[xindex] = 1.0; // not used
		bx[xindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_x[xindex];

		for (xindex = 1; xindex < g.ngps_x() - 1; xindex++) {
			index = g.index(xindex, yindex, zindex);
			index_xp = g.index(xindex + 1, yindex, zindex);
			index_xm = g.index(xindex - 1, yindex, zindex);

			ax[xindex] = aa - 0.5 * timestephalfoversixstaticpot_x[xindex + 1];
			cx[xindex] = cc - 0.5 * timestephalfoversixstaticpot_x[xindex - 1];
			bx[xindex] = bb
					- 0.5 * fivetimestephalfoverthreestaticpot_x[xindex];

		};

		xindex = g.ngps_x() - 1;
		x = g.x(xindex);
		index = g.index(xindex, yindex, zindex);
		index_xm = g.index(xindex - 1, yindex, zindex);
		ax[xindex] = aa
				- 0.5 * timestephalfoversix
						* hamil.scalarpotx(x + g.delt_x(), y, z, time);
		cx[xindex] = cc - 0.5 * timestephalfoversixstaticpot_x[xindex - 1];
		bx[xindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_x[xindex];

		rhstwo_x.solve(ax, bx, cx, rhsone_x, g.ngps_x());

		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_x[xindex];
		};

	};

	// ==================================== W_y =================================

	oneoverhsquare = 1.0 / (g.delt_y() * g.delt_y());
	aa = -OOS - 0.5 * timestephalf * oneoverhsquare * diffconst;
	cc = aa;
	bb = -FOT + 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

	// ---------- W_-y
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		yindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_yp = g.index(xindex, yindex + 1, zindex);
		rhsone[index] = (aa + 0.5 * timestephalfoversixstaticpot_y[yindex + 1])
				* start[index_yp]
				+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_y[yindex])
						* start[index];

		for (yindex = 1; yindex < g.ngps_y() - 1; yindex++) {
			index = g.index(xindex, yindex, zindex);
			index_yp = g.index(xindex, yindex + 1, zindex);
			index_ym = g.index(xindex, yindex - 1, zindex);
			rhsone[index] = (aa
					+ 0.5 * timestephalfoversixstaticpot_y[yindex + 1])
					* start[index_yp]
					+ (bb + 0.5 * fivetimestephalfoverthreestaticpot_y[yindex])
							* start[index]
					+ (cc + 0.5 * timestephalfoversixstaticpot_y[yindex - 1])
							* start[index_ym];
		};

		yindex = g.ngps_y() - 1;
		index = g.index(xindex, yindex, zindex);
		index_ym = g.index(xindex, yindex - 1, zindex);
		rhsone[index] =
				(bb + 0.5 * fivetimestephalfoverthreestaticpot_y[yindex])
						* start[index]
						+ (cc + 0.5 * timestephalfoversixstaticpot_y[yindex - 1])
								* start[index_ym];

	};

	// --------------- W_+y
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		x = g.x(xindex);
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_y[yindex] = rhsone[index];
		};
		aa = -OOS + 0.5 * timestephalf * oneoverhsquare * diffconst;
		cc = aa;
		bb = -FOT - 0.5 * timestephalf * (2.0 * oneoverhsquare) * diffconst;

		yindex = 0;
		y = g.y(yindex);
		index = g.index(xindex, yindex, zindex);
		index_yp = g.index(xindex, yindex + 1, zindex);
		ay[yindex] = aa - 0.5 * timestephalfoversixstaticpot_y[yindex + 1];
		cy[yindex] = 1.0; // not used
		by[yindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_y[yindex];

		for (yindex = 1; yindex < g.ngps_y() - 1; yindex++) {
			index = g.index(xindex, yindex, zindex);
			index_yp = g.index(xindex, yindex + 1, zindex);
			index_ym = g.index(xindex, yindex - 1, zindex);

			ay[yindex] = aa - 0.5 * timestephalfoversixstaticpot_y[yindex + 1];
			cy[yindex] = cc - 0.5 * timestephalfoversixstaticpot_y[yindex - 1];
			by[yindex] = bb
					- 0.5 * fivetimestephalfoverthreestaticpot_y[yindex];

		};

		yindex = g.ngps_y() - 1;
		y = g.y(yindex);
		index = g.index(xindex, yindex, zindex);
		index_ym = g.index(xindex, yindex - 1, zindex);
		ay[yindex] = aa
				- 0.5 * timestephalfoversix
						* hamil.scalarpot(x, y + g.delt_y(), z, time);
		cy[yindex] = cc - 0.5 * timestephalfoversixstaticpot_y[yindex - 1];
		by[yindex] = bb - 0.5 * fivetimestephalfoverthreestaticpot_y[yindex];

		rhstwo_y.solve(ay, by, cy, rhsone_y, g.ngps_y());

		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_y[yindex];
		};

	};

	// ============================= S_x ========================================

	vecpotwithprefactor = (
			timestep / (8.0 * g.delt_x()) * hamil.vecpot_x(time));

	aaa = OOS - vecpotwithprefactor;
	ccc = OOS + vecpotwithprefactor;
	bbb = TOT;

	// Calculate the rhs vector S_-x *this

	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		xindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_xp = g.index(xindex + 1, yindex, zindex);

		rhsone[index] = aaa * start[index_xp]
				+ ((4.0 + lambda) / 6.0 - lambda * vecpotwithprefactor)
						* start[index];

		for (xindex = 1; xindex < g.ngps_x() - 1; xindex++) {
			index = g.index(xindex, yindex, zindex);
			index_xp = g.index(xindex + 1, yindex, zindex);
			index_xm = g.index(xindex - 1, yindex, zindex);

			rhsone[index] = aaa * start[index_xp] + bbb * start[index]
					+ ccc * start[index_xm];
		};

		xindex = g.ngps_x() - 1;
		index = g.index(xindex, yindex, zindex);
		index_xm = g.index(xindex - 1, yindex, zindex);

		rhsone[index] = ccc * start[index_xm]
				+ ((4.0 + lambda) / 6.0 - llambda * vecpotwithprefactor)
						* start[index];
	}

	// The matrix  S_+x
	for (yindex = 0; yindex < g.ngps_y(); yindex++) {
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_x[xindex] = rhsone[index];
		};

		aaa = OOS + vecpotwithprefactor;
		ccc = (OOS - vecpotwithprefactor);
		bbb = TOT;

		b_upperleft = ((4.0 + lambda) / 6.0 + lambda * vecpotwithprefactor);
		b_lowerright = ((4.0 + lambda) / 6.0 + llambda * vecpotwithprefactor);

		rhstwo_x.solve_toep(aaa, bbb, b_upperleft, b_lowerright, ccc, rhsone_x,
				g.ngps_x());

		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_x[xindex];
		};
	};

	// ============================= S_y ========================================

	vecpotwithprefactor = (
			timestep / (8.0 * g.delt_y()) * hamil.vecpot_y(time));

	aaa = OOS - vecpotwithprefactor;
	ccc = OOS + vecpotwithprefactor;
	bbb = TOT;

	// Calculate the rhs vector S_-y *this

	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		yindex = 0;
		index = g.index(xindex, yindex, zindex);
		index_yp = g.index(xindex, yindex + 1, zindex);

		rhsone[index] = aaa * start[index_yp]
				+ ((4.0 + lambda) / 6.0 - lambda * vecpotwithprefactor)
						* start[index];

		for (yindex = 1; yindex < g.ngps_y() - 1; yindex++) {
			index = g.index(xindex, yindex, zindex);
			index_yp = g.index(xindex, yindex + 1, zindex);
			index_ym = g.index(xindex, yindex - 1, zindex);

			rhsone[index] = aaa * start[index_yp] + bbb * start[index]
					+ ccc * start[index_ym];
		};

		yindex = g.ngps_y() - 1;
		index = g.index(xindex, yindex, zindex);
		index_ym = g.index(xindex, yindex - 1, zindex);

		rhsone[index] = ccc * start[index_ym]
				+ ((4.0 + lambda) / 6.0 - llambda * vecpotwithprefactor)
						* start[index];
	}

	// The matrix  S_+y
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			rhsone_y[yindex] = rhsone[index];
		};

		aaa = OOS + vecpotwithprefactor;
		ccc = (OOS - vecpotwithprefactor);
		bbb = TOT;

		b_upperleft = ((4.0 + lambda) / 6.0 + lambda * vecpotwithprefactor);
		b_lowerright = ((4.0 + lambda) / 6.0 + llambda * vecpotwithprefactor);

		rhstwo_y.solve_toep(aaa, bbb, b_upperleft, b_lowerright, ccc, rhsone_y,
				g.ngps_y());

		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			index = g.index(xindex, yindex, zindex);
			start[index] = rhstwo_y[yindex];
		};
	};

}

void density::dump_to_file(grid g, FILE* os, int stepwidth) {

	long xindex, yindex, zindex, i, counter, counter_ii, counter_iii;
	double u, r, rho, z, legpolnew, legpolold, legpololder;
	dcomplex summingres;

	counter = 0;
	counter_ii = 0;
	counter_iii = 0;

	switch (g.dimens()) {
	case 2:
		counter = 0;
		for (xindex = 0; xindex < g.ngps_x(); xindex += stepwidth) {
			counter++;
			counter_ii = 0;
			for (yindex = 0; yindex < g.ngps_y(); yindex += stepwidth) {
				counter_ii++;
				i = g.index(xindex, yindex, 0);
				fprintf(os, "%e %e\n", real(start[i]), imag(start[i]));
			};
		}
		;

		break;

	default:
		for (long i = 0; i < dens_dim; i++) {
			fprintf(os, "%e %e\n", real(start[i]), imag(start[i]));
		}
		;
	};
}

void density::calculate_fixed_potential_array(grid g, hamop hamil,
		double time) {
	double x, y, z;
	long index, xindex, yindex, zindex;
	dcomplex imagi(0.0,1.0);

	switch (g.dimens()) {

	case 2:
		zindex = 0;
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			x = g.x(xindex);
			for (yindex = 0; yindex < g.ngps_y(); yindex++) {
				y = g.y(yindex);
				index = g.index(xindex, yindex, zindex);
				start[index] = hamil.scalarpot(x, y, z, time)
						- imagi
								* hamil.imagpot(xindex, yindex, zindex, time,
										g);
			};
		}
		;
		break;

	};

}

void density::calculate_fixed_potential_array_x(grid g, hamop hamil,
		double time, const density &densstart, const density &dens) {
	double x, y, z;
	long xindex, yindex, zindex, index;
	dcomplex imagi(0.0,1.0);

	density densx(g.ngps_x());
	densx.nullify();
	density densxstart(g.ngps_x());
	densxstart.nullify();
double k = 0.3;
	switch (g.dimens()) {
	case 2:
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			for (yindex = 0; yindex < g.ngps_y(); yindex++) {
				index = g.index(xindex, yindex, 0);
				densx[xindex] = densx[xindex]
						+ real(dens[index])* g.delt_y();

				densxstart[xindex] = densxstart[xindex]
										+ real(densstart[index])
												* g.delt_y();

			}

		}

		yindex = 0;
		zindex = 0;
		y = g.y(yindex);
		z = g.z(zindex);
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			x = g.x(xindex);
			start[xindex] = hamil.scalarpotx(x, y, z, time)

					- imagi * hamil.imagpotx(xindex, yindex, zindex, time, g);
		}
		;
		break;

	}
}

void density::calculate_fixed_potential_array_y(grid g, hamop hamil,
		double time, const density &densstart, const density &dens) {
	double x, y, z;
	long xindex, yindex, zindex, index;
	dcomplex imagi(0.0,1.0);
	density densy(g.ngps_x());
	densy.nullify();
	density densystart(g.ngps_x());
	densystart.nullify();
	double k = 0.4;
	switch (g.dimens()) {
	case 2:
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			for (yindex = 0; yindex < g.ngps_y(); yindex++) {
				index = g.index(xindex, yindex, 0);
				densy[yindex] = densy[yindex]
						+ real(dens[index]) * g.delt_x();

				densystart[yindex] = densystart[yindex]
										+ real(densstart[index])
												* g.delt_x();

			}

		}


		xindex = 0;
		yindex = 0;
		x = g.x(xindex);
		z = g.z(zindex);
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			y = g.y(yindex);
			start[yindex] = hamil.scalarpoty(x, y, z, time)

					- imagi * hamil.imagpoty(xindex, yindex, zindex, time, g);
		}
		;
		break;

	}
}

void density::calculate_fixed_potential_array_xy(grid g, hamop hamil,
		const density &dens, double time) {
	double x, y, z;
	long index, xindex, yindex, zindex;

	double k = 0.02;
	switch (g.dimens()) {
	case 2:

		zindex = 0;
		z = g.z(zindex);
		for (xindex = 0; xindex < g.ngps_x(); xindex++) {
			x = g.x(xindex);
			for (yindex = 0; yindex < g.ngps_y(); yindex++) {
				y = g.y(yindex);
				index = g.index(xindex, yindex, zindex);
				start[index] = hamil.interactionpotxy(x, y, z, time);

			};
		}
		;
		break;

	};
}

void density::calculate_nonlinearpot(grid g, hamop hamil, double startnorm,
		const density &dens, double time, double k) {
	double x, y, z;
	long index, xindex, yindex, zindex;
	double initialnorm;
	for (xindex = 0; xindex < g.ngps_x(); xindex++) {
		x = g.x(xindex);
		for (yindex = 0; yindex < g.ngps_y(); yindex++) {
			y = g.y(yindex);
			index = g.index(xindex, yindex, zindex);
			start[index] = -hamil.mitosis(time) * (k * startnorm - dens[index]);

		};
	};

}

double density::expect_x(grid g) {
	dcomplex result(0.0,0.0);
	long xindex, yindex, zindex;
	long index;
	double x;

	switch (g.dimens())
	{

		case 2 :

		for (xindex=0; xindex<g.ngps_x(); xindex++)
		{
			x=g.x(xindex);

			for (yindex=0; yindex<g.ngps_y(); yindex++)
			{
				index=g.index(xindex,yindex,0);

				result=result+x*x*real(start[index])*g.delt_x()*g.delt_y();
			};
		};

	};

	return real(result);

}

double density::expect_y(grid g) {
	dcomplex result(0.0,0.0);
	long xindex, yindex, zindex;
	long index;
	double y;

	switch (g.dimens())
	{

		case 2 :
		for (yindex=0; yindex<g.ngps_y(); yindex++)
		{
			y=g.y(yindex);
			for (xindex=0; xindex<g.ngps_x(); xindex++)
			{

				index=g.index(xindex,yindex,0);

				result=result+y*y*real(start[index])*g.delt_x()*g.delt_y();
			};
		};
	};

	return real(result);

}

double density::norm(grid g) {
	dcomplex result(0.0,0.0);
	long xindex, yindex, zindex;
	long index,indexm,indexmm;
	switch(g.dimens()) {
		case 2:
		for (xindex=0; xindex<g.ngps_x(); xindex++)
		{

			for (yindex=0; yindex<g.ngps_y(); yindex++)
			{
				index=g.index(xindex,yindex,0);

				result=result+real(start[index])*g.delt_x()*g.delt_y();
			};

		};
		break;

	}

	return real(result);

}
double density::normx(grid g) {
	dcomplex result(0.0,0.0);
	long xindex, yindex, zindex;
	long index,indexm,indexmm;

	for (xindex=0; xindex<g.ngps_x(); xindex++)
	{

		result=result+real(start[xindex])*g.delt_x();
	};

	return real(result);

}

double density::normy(grid g) {
	dcomplex result(0.0,0.0);
	long xindex, yindex, zindex;
	long index,indexm,indexmm;

	for (yindex=0; yindex<g.ngps_y(); yindex++)
	{

		result=result+real(start[yindex])*g.delt_y();
	};

	return real(result);

}

void density::solve(const density &a, const density &b, const density &c,
		const density &psi, long dimens) {

	density e(dimens), f(dimens);
	long i;

	e[0] = -b[0] / a[0];
	f[0] = psi[0] / a[0];
	for (i = 1; i < dimens; i++) {
		e[i] = (-c[i] / e[i - 1] - b[i]) / a[i];
		f[i] = (psi[i] + c[i] / e[i - 1] * f[i - 1]) / a[i];
	};

	start[dimens - 1] = -f[dimens - 1] / e[dimens - 1];
	for (i = dimens - 2; i >= 0; i--) {
		start[i] = (start[i + 1] - f[i]) / e[i];
	};

}

// this is for pseudo-Toeplitz structured matrices, i.e., a, c is always the sa
// b also, apart from the corner elements. Thus b_upperleft and b_lowerright are extra arguments.
// version with a, b, c, double !!!
void density::solve_toep(double a, double b, double b_upperleft,
		double b_lowerright, double c, const density &psi, long dimens) {

	density f(dimens);
	fluid e(dimens);
	long i;
	double covera, bovera;

	covera = c / a;
	bovera = b / a;

	e[0] = -b_upperleft / a;
	f[0] = psi[0] / a;
	for (i = 1; i < dimens - 1; i++) {
		e[i] = -covera / e[i - 1] - bovera;
		f[i] = psi[i] / a + covera / e[i - 1] * f[i - 1];
	};
	e[dimens - 1] = -covera / e[dimens - 2] - b_lowerright / a;
	f[dimens - 1] = psi[dimens - 1] / a
			+ covera / e[dimens - 2] * f[dimens - 2];

	start[dimens - 1] = -f[dimens - 1] / e[dimens - 1];
	for (i = dimens - 2; i >= 0; i--) {
		start[i] = (start[i + 1] - f[i]) / e[i];
	};

}

void density::apply(const density &a, const density &b, const density &c,
		const density &psi) {

	long dimens = psi.dens_size();
	long i;

	start[0] = a[0] * psi[1] + b[0] * psi[0];

	for (i = 1; i < dimens - 1; i++) {
		start[i] = a[i] * psi[i + 1] + b[i] * psi[i] + c[i] * psi[i - 1];
	};

	start[dimens - 1] = b[dimens - 1] * psi[dimens - 1]
			+ c[dimens - 1] * psi[dimens - 2];

}

