#ifndef density_h
#define density_h density_h
#include<assert.h>
#include<complex>
#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>

using namespace std;

class fluid;
class hamop;
class grid;

class density
{
  public:
    density(long x=0) : dens_dim(x), start(new complex<double>[x]) { }

    density(const density& v)
      {
	dens_dim  = v.dens_dim;
	start = new complex<double>[dens_dim];
	for (long i = 0; i < dens_dim; i++)
	  start[i] = v.start[i];
      }

    virtual ~density() { delete [] start;}

    long  dens_size()const {return dens_dim;}

    complex<double>&  operator[](long index)
    {
      //      assert(index >= 0  &&  index < dens_dim);
      return start[index];
    }

    const complex<double>&  operator[](long index) const
    {
      //      assert(index >= 0  &&  index < dens_dim);
      return start[index];
    }

    density& operator=(const density&);
    density& operator=(const fluid&);

    void init(grid g, int inittype, double width, int k, double time, FILE* filename, int output_of_interest);
    void init(grid g, int inittype, double width, int k, double time);
    complex<double>* begin() { return start;}
    complex<double>* end()   { return start + dens_dim;}
    const complex<double>* begin() const { return start;}
    const complex<double>* end()   const { return start + dens_dim;}

    density Static_KS(grid g);
    density& operator *= (double z);
    density& operator *= (complex<double> z);


    void regrid(grid g, grid g_small, const density &v);

    void propagator(double delta_time,
		   double time,
		   double diffconst,
		   grid g,
		   hamop hamil,
		   const density &staticpot_x,
		   const density &staticpot_y,
		   const density &staticpot_xy
		   );

   void propagatedft(double delta_time,
		   double time,
		   double diffconst,
		   grid g,
		   hamop hamil,
		   const density &staticpot,
		   const density &tddftpot
		   );

    complex<double> energy(double time, grid g, hamop hamil,  const double masses[], const density &staticpot_x,
    		const density &staticpot_y, const density &staticpot_xy );
    complex<double> energydft(double time, grid g, hamop hamil, const double masses[], const density &staticpot,
    		const density &tddftpot );

    double expect_x(grid g);
    double expect_y(grid g);



    double norm(grid g);
    double normx(grid g);
    double normy(grid g);

    double non_ionized(grid g, long box);
    double sing_ionized(grid g, long box);
    double doub_ionized(grid g, long box);
    double molecule_ionized(grid g, long box);


    void solve(const density &a, const density &b,
	       const density &c, const density &psi, long dimens);
    void solve_toep(double a, double b, double b_upperleft, double b_lowerright,
	       double c, const density &psi, long dimens);

    void apply(const density &a, const density &b,
	       const density &c, const density &psi);

    void dump_to_file(grid g, FILE* os, int stepwidth);
    void  calculate_fixed_potential_array(grid g, hamop hamil, double time);
    void  calculate_fixed_potential_arraydft(grid g, hamop hamil, double time);

    void nullify();

    void calculate_fixed_potential_array_x(grid g, hamop hamil, double time, const density &densstart, const density &dens);
    void calculate_fixed_potential_array_y(grid g, hamop hamil, double time, const density &densstart, const density &dens);
     void calculate_fixed_potential_array_z(grid g, hamop hamil, double time);
    void calculate_fixed_potential_array_xy(grid g, hamop hamil, const density &dens, double time);
    void calculate_nonlinearpot(grid g, hamop hamil, double startnorm, const density &dens, double time, double k);
    private:
    long  dens_dim;
    complex<double>   *start;

    complex<double> int_P1_lp(const long l, const long l_l, grid g);
    complex<double> int_P2_lp(const long l, const long l_l, grid g);



    void do_cn_step_x(double timestep, double time,double diffconst, grid g, hamop hamil,
    		long yindex, long zindex, const fluid &dens_one, const density &dens_two);

    void do_cn_step_x_muller(double timestep, double time, double diffconst, grid g, hamop hamil,
    		const density &staticpot, const density &tddftpot, long yindex, long zindex);

    void do_cn_step_xy_muller(double timestep, double time,double diffconst, grid g, hamop hamil,
    		const density &staticpot_x, const density &staticpot_y, const density &staticpot_xy, long yindex, long zindex);

    void do_muller(double timestep, double time, double diffconst, grid g, hamop hamil,
    		const density &staticpot,  long noofgridpoints[] );

    void do_cn_step_xy(double timestep, double time, double diffconst, grid g, hamop hamil,
    		long zindex, const fluid &dens_one, const density &dens_two);

};


ostream& operator<<(ostream& os, const density& v);
istream& operator>>(istream& is, density& v);

complex<double> operator * (const density &v, const density &w );
density operator * (double z, const density &v);
density operator * (complex<double> z, const density &v);
density operator * (const density &v, double z);
density operator * (const density &v, complex<double> z);
density operator / (const density &v, double z);
density operator / (const density &v, const density &w);
density operator | (const density &v, const density &w);

density operator - (const density &v, const density &w );
density operator + (const density &v, const density &w );
density operator + (const density &v, const fluid &w );
density operator + (const fluid &v, const density &w );
density operator + (const density &v, double z );

#endif // density_h






