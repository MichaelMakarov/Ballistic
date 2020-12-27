#include "PlahovIntegrator.h"
#include "RungeKuttaIntegrator.h"
#include "AdamsBashforthIntegrator.h"
#include "GeoPotential.h"
#include "EGM96.h"
#include "AstroValues.h"
#include "DateTime.h"
#include <iostream>
#include <iomanip>

using namespace ball;

void TestIntegrator();

int main()
{
	std::cout << "...Testing integrators...\n";
	TestIntegrator();
	return 0;
}

void TestIntegrator()
{
	std::cout << std::setprecision(16);
	types::PV x0(
		geometry::XYZ(-4688980.289, -11060428.914, 238914.750),
		geometry::XYZ(-1402.353, 1753.937, 5249.663));
	time::JD t0(time::DateTime(2015, 7, 11, 0, 0, 0, 0));
	auto f = [](const types::PV& x, const double tk) {
		tasks::GeoPotential geopot(
			std::move(std::unique_ptr<tasks::IPotential>{ new tasks::PotentialEGM96() }),
			tasks::egm96::R,
			tasks::egm96::Mu,
			16);
		auto rbl = tasks::GCS_OrthoToSpher(geometry::XYZ(x.P1, x.P2, x.P3));
		double pot = -geopot(rbl) / rbl.R;
		auto xyz = tasks::GCS_SpherToOrtho(geometry::RBL(pot, rbl.B, rbl.L));
		return types::PV(geometry::XYZ(x.V1, x.V2, x.V3), xyz);
	};

	std::cout << "Initial x0: " << x0 << std::endl;

	double step = 360.0 / time::SEC_PER_DAY;
	time::JD tk = t0;
	size_t n = 10;
	std::vector<std::pair<types::PV, double>> values(n);
	values[0] = { x0, tk };
	PlahovIntegrator plint;
	plint.Initialize(x0, tk);
	plint.Function(f);
	RungeKuttaIntegrator rkint;
	rkint.Initialize(x0, tk);
	rkint.Function(f);
	AdamsBashforthIntegrator abint;
	abint.Function(f);

	std::cout << "Plahov's integrator:\n";
	for (size_t i = 1; i < n; ++i)
	{
		auto xk = plint.Integrate(step);
		tk.AddDays(step);
		values[i] = { xk, tk };
		plint.Initialize(xk, tk);
		std::cout << i << ": " << xk << std::endl;
	}
	tk = t0;
	std::cout << "Runge-Kutta's integrator:\n";
	for (size_t i = 1; i < n; ++i)
	{
		auto xk = rkint.Integrate(step);
		tk.AddDays(step);
		rkint.Initialize(xk, tk);
		std::cout << i << ": " << xk << std::endl;
	}
	std::cout << "Adams-Bashforth's integrator:\n";
	for (size_t i = 8; i < n; ++i)
	{
		abint.Initialize({
			values[i - 8],
			values[i - 7],
			values[i - 6],
			values[i - 5],
			values[i - 4],
			values[i - 3],
			values[i - 2],
			values[i - 1]
		});
		auto xk = abint.Integrate(step);
		tk.AddDays(step);
		values[i] = { xk, tk };
		std::cout << i + 1 << ": " << xk << std::endl;
	}

}
