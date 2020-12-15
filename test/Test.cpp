#include "PlahovIntegrator.h"
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
		double pot = geopot(rbl) / rbl.R;
		auto xyz = tasks::GCS_SpherToOrtho(geometry::RBL(pot, rbl.B, rbl.L));
		return types::PV(xyz, geometry::XYZ());
	};
	PlahovIntegrator integr;
	integr.Initialize(x0, t0);
	integr.Function(f);
	t0.AddSeconds(60);
	auto xk = integr.Integrate(t0);
	std::cout << xk << std::endl;
}
