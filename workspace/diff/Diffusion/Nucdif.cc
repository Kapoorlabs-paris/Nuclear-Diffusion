#include"Nucdif.h"
#include<stdio.h>
#include"density.h"
#include"fluid.h"
#include"grid.h"
#include"hamop.h"
#include<cmath>

#define pi 3.1415926

int main(int argc, char **argv) {
	FILE *file_densdat;
	FILE *file_obser;
	FILE *file_reading;
	FILE *file_potential;

	// *** create some files with appropriate appendices

	char string_densdat[] = "res/densitytestkg.dat";
	char string_obser[] = "res/observtestkg.dat";
	char string_potential[] = "res/potentialtestkg.dat";
	char string_reading[] = "res/density_ini.dat";

	file_densdat = fopen(string_densdat, "w");
	file_obser = fopen(string_obser, "w");
	file_potential = fopen(string_potential, "w");
	file_reading = fopen(string_reading, "r");

	double diffconst = 0.1;
	double deltx = 0.5;
	double delty = 0.5;
	double deltz = 1;
	long ngpsx = 1500;
	long ngpsy = 1500;
	long ngpsz = 1;

	// *** declare grid ***
	grid g;
	g.set_dim(2);      // propagation mode is 2 for xy cartesian in 2 dimensions
	g.set_ngps(ngpsx, ngpsy, 1);    // N_x, N_y, 1
	g.set_delt(deltx, delty, deltz);  // delta_x, delta_y, 1.0
	g.set_offs(ngpsx / 2, ngpsy / 2, 0);    // origin at N_x/2, N_y/2

	// *** declare smaller grid for reading***
	grid g_small;
	g_small.set_dim(2);
	g_small.set_ngps(ngpsx / 2, ngpsy / 2, 1);
	g_small.set_delt(deltx, delty, deltz);
	g_small.set_offs(ngpsx / 4, ngpsy / 4, 0);

	// *** declare rest of variables ***
	double timestep = 0.1; //0.1*(deltx*deltx+delty*delty)/(2*diffconst); // Make sure you do not choose a very big time-step for small diffconst
	long no_of_timesteps = 5000;
	int obs_output_every = 1;
	long dens_output_every = 49;
	int dumpingstepwidth = 1;

	hamop hamilton(g, vecpot_x, vecpot_y, vecpot_z, scalarpotx, scalarpoty,
			scalarpotz, interactionpotxy, imagpotx, imagpoty, mitosis);

	density dens(g.ngps_x() * g.ngps_y() * g.ngps_z());
	density densold(g.ngps_x() * g.ngps_y() * g.ngps_z());
	density densstart(g.ngps_x() * g.ngps_y() * g.ngps_z());
	density densread(g_small.ngps_x() * g_small.ngps_y() * g_small.ngps_z());

	double time = 0.0;

	long counter_i = 0;
	long counter_ii = 0;
	long counter_iii = 0;

	// initialization
	long outputofinterest = 0;
	//	dens.nullify();
	//densread.init(g_small,99,0.1,0.0,0.0,file_reading,outputofinterest);
	//dens.regrid(g,g_small,densread);
	dens.init(g, 1, 20, 10, 0.0);
	dens *= 1.0 / (dens.norm(g));
	densold = dens;
	dens = densold;
	densstart = dens;
	fclose(file_reading);

	cout << "norm dens    : " << dens.norm(g) << "\n";

	density staticpot_x(g.ngps_x());
	density staticpot_xold(g.ngps_x());
	staticpot_x.calculate_fixed_potential_array_x(g, hamilton, 0.0, densstart,
			densstart);

	density staticpot_y(g.ngps_y());
	density staticpot_yold(g.ngps_y());
	staticpot_y.calculate_fixed_potential_array_y(g, hamilton, 0.0, densstart,
			densstart);

	density staticpot_xy(g.ngps_x() * g.ngps_y() * g.ngps_z());
	staticpot_xy.calculate_fixed_potential_array_xy(g, hamilton, dens, 0.0);


	dens.dump_to_file(g, file_densdat, dumpingstepwidth);

	long ts;


	// *************  timeprop

	for (ts = 0; ts < no_of_timesteps; ts++) {
		counter_i++;
		counter_ii++;
		time = (timestep * ts);

		cout << "Time : " << ts << "  " << "norm dens    : " << dens.norm(g)
				<< endl;
		staticpot_x.calculate_fixed_potential_array_x(g, hamilton, time, densstart,
					dens);
		staticpot_y.calculate_fixed_potential_array_y(g, hamilton, time, densstart,
					dens);
		staticpot_xy.calculate_fixed_potential_array_xy(g, hamilton, dens, time);

		 dens.propagator(timestep, time, diffconst,g, hamilton,staticpot_x, staticpot_y,staticpot_xy);

		//dens *= 1.0 / (dens.norm(g));


		// Non-linear part is included in staticpot_xy.
		if (counter_ii == obs_output_every) {

			fprintf(file_obser,
					"%.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14lex \n",
					(time + (timestep)), dens.norm(g), vecpot_x(time),
					vecpot_y(time), dens.expect_x(g), dens.expect_y(g),
					staticpot_x.normx(g), staticpot_y.normy(g));

			counter_ii = 0;
		};

		if (counter_i == dens_output_every) {
			dens.dump_to_file(g, file_densdat, dumpingstepwidth);

			counter_i = 0;
		};

	};

	fclose(file_obser);
	fclose(file_densdat);
	fclose(file_potential);

	cout << "Now plot the results." << endl;

}


// Declare some parameters for the extra-terms in the diffusion equation
	double frequ = 8;
	double alphahat = 5;
	double n = 60.0;
	double ww = 0.5 * frequ / n;
	double lambda = -1;


double vecpot_x(double time) {
	double result;

	result = alphahat *sin(ww*time)*sin(ww*time)*sin(frequ * time);

	//*sin(ww*time)*sin(ww*time)*
	return 0*result;
}

double vecpot_y(double time) {

	return vecpot_x(time);
}

double vecpot_z(double time)

{
	return 0.0;
}

double scalarpotx(double x, double y, double z, double time) {

	double result;
	double amp = alphahat *sin(frequ * time);

    result =   0.1/sqrt((x + amp)*(x + amp) + 5);




	return result;

}

double scalarpoty(double x, double y, double z, double time) {

	double result;

	double amp = alphahat *sin(frequ * time);
	double sum = 0;


	 result =   0.1/sqrt((y + amp) * (y + amp) + 5);
	return result;

}

double scalarpotz(double x, double y, double z, double time) {
	double result = 0.0;
	return result;
}

double interactionpotxy(double x, double y, double z, double time) {

	double result;

	result = 0.0001/sqrt((x-y)*(x-y) + 10) + 0.0001/sqrt((x+y)*(x+y) + 10)  ;




	return result;
}
// For the lambda part of the non-linear term, lambda*density*(1-density), rest is in the function calculate_staticpot_xy.
double mitosis(double time) {
	double lambda = 0.0;
	return lambda;

}

// For absorbing boundaries ampl>0 for reflecting boundaries ampl = 0
// Absorbing in X direction
double imagpotx(long xindex, long yindex, long zindex, double time, grid g) {
	double x, y, z;
	//     double ampl=0.0; // switch imaginary potential off
	double ampl = 50.0; // switch imaginary potential on

	if (ampl > 1.0) {
		x = ((double) xindex + 0.5 - 0.5 * g.ngps_x()) / (0.5 * g.ngps_x())
				* ((double) xindex + 0.5 - 0.5 * g.ngps_x())
				/ (0.5 * g.ngps_x());

		return ampl * x * x * x * x * x * x * x * x;
	} else {
		return 0.0;
	};

}

// Absorbing in Y direction
double imagpoty(long xindex, long yindex, long zindex, double time, grid g) {
	double y;
	//  double ampl=0.0; // switch imaginary potential off
	double ampl = 50.0; // switch imaginary potential on

	if (ampl > 1.0) {
		y = ((double) yindex + 0.5 - 0.5 * g.ngps_y()) / (0.5 * g.ngps_y())
				* ((double) yindex + 0.5 - 0.5 * g.ngps_y())
				/ (0.5 * g.ngps_y());
		return ampl * y * y * y * y * y * y * y * y;
	} else {
		return 0.0;
	};

}

