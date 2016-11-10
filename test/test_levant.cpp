
// #define DEBUG 1  	

#include "levant_differentiator/levant_differentiator.hpp"

int main(int argc, char const *argv[])
{

	// double L = 1000.0;
	// levant< double > differtentiator;
	Levant< double > differtentiator;
	differtentiator.setParameters(2.0, 1000.0, { 3.0, 1.5, 1.1});


	double function;
	double functionp;
	double A = 1.0;
	double f = 100.0;
	double dt = 1.0 / f;

	// differtentiator.setFilter(true, 10, dt);

	// boost::numeric::odeint::runge_kutta4< std::vector<double> > stepper;
	std::vector< double > x0(3, 0);

	for (double t = 0.0; t < 10.0; t += dt)
	{

		function = A * std::cos(2 * M_PI * t);
		functionp = -1.0 * A * 2 * M_PI * std::sin(2 * M_PI * t);
		// stepper.do_step(differtentiator, x0, function, dt);
		differtentiator.step(function, dt);
		x0 = differtentiator.getX();
		std::cout << t << "\t" << function << "\t" << functionp << "\t" << x0[0] << "\t" << x0[1] << std::endl;
	}

	return 0;
}