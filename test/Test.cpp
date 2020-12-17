#include "PlahovIntegrator.h"
#include "RungeKuttaIntegrator.h"
#include "GeoPotential.h"
#include "EGM96.h"
#include "AstroValues.h"
#include <iostream>

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
	types::PV x0(
		geometry::XYZ(-4688980.289, -11060428.914, 238914.750),
		geometry::XYZ(-1402.353, 1753.937, 5249.663));
	time::JD t0(time::DateTime(2015, 7, 11, 0, 0, 0, 0));
	auto f = [](const types::PV& x, const time::JD tk) {
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
	double step = 60.0 / time::SecPerDay;
	PlahovIntegrator plint;
	plint.Initialize(x0, t0);
	plint.Function(f);
	auto xk = plint.Integrate(step);
	std::cout << "Plahov's integrator: " << xk << std::endl;

	RungeKuttaIntegrator rkint;
	rkint.Initialize(x0, t0);
	rkint.Function(f);
	xk = rkint.Integrate(step);
	std::cout << "Runge-Kutta's integrator: " << xk << std::endl;
}
