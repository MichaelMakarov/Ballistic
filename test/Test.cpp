#include "PlahovIntegrator.h"
#include "RungeKuttaIntegrator.h"
#include "AdamsIntegrator.h"
#include "GeoPotential.h"
#include "EGM96.h"
#include "PZ90.h"
#include "StaticAtmosphere81.h"
#include "AstroValues.h"
#include "DateTime.h"
#include <Ballistic.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace ball;

void TestIntegrator();
void TestBallistic();

int main()
{
	TestBallistic();
	return 0;
}

void TestIntegrator()
{
	std::cout << "...Testing integrators...\n";

	std::cout << std::setprecision(16);
	types::PV x0(
		geometry::XYZ(-4688980.289, -11060428.914, 238914.750),
		geometry::XYZ(-1402.353, 1753.937, 5249.663));
	time::JD t0(time::DateTime(2015, 7, 11, 0, 0, 0, 0));
	auto f = [](const types::PV& x, const time::JD& tk) {
		tasks::GeoPotential geopot(
			std::move(std::unique_ptr<tasks::IGravity>{ new tasks::EGM96() }),
			16);
		auto rbl = tasks::GCS_OrthoToSpher(geometry::XYZ(x.P1, x.P2, x.P3));
		double pot = -geopot(rbl) / rbl.R;
		auto xyz = tasks::GCS_SpherToOrtho(geometry::RBL(pot, rbl.B, rbl.L));
		return types::PV(geometry::XYZ(x.V1, x.V2, x.V3), xyz);
	};

	std::cout << "Initial x0: " << x0 << std::endl;

	int step = 360;
	time::JD tk = t0;
	size_t n = 10;
	std::vector<std::pair<types::PV, time::JD>> values(n);
	values[0] = { x0, tk };
	PlahovIntegrator<> plint;
	plint.Initialize(x0, tk);
	plint.Func = f;
	RungeKuttaIntegrator<> rkint;
	rkint.Initialize(x0, tk);
	rkint.Func = f;
	AdamsIntegrator<> abint;
	abint.Func = f;

	std::cout << "Plahov's integrator:\n";
	for (size_t i = 1; i < n; ++i)
	{
		auto xk = plint.Integrate(step);
		tk.AddSeconds(step);
		values[i] = { xk, tk };
		plint.Initialize(xk, tk);
		std::cout << i << ": " << xk << std::endl;
	}
	tk = t0;
	std::cout << "Runge-Kutta's integrator:\n";
	for (size_t i = 1; i < n; ++i)
	{
		auto xk = rkint.Integrate(step);
		tk.AddSeconds(step);
		rkint.Initialize(xk, tk);
		std::cout << i << ": " << xk << std::endl;
	}
	std::cout << "Adams-Bashforth's integrator:\n";
	for (size_t i = 8; i < n; ++i)
	{
		abint.Initialize(values.begin());
		auto xk = abint.Integrate(step);
		tk.AddSeconds(step);
		values[i] = { xk, tk };
		std::cout << i + 1 << ": " << xk << std::endl;
	}

}
void TestBallistic()
{
	using namespace time;
	using namespace types;

	std::cout << "\n...Ballistic test...\n";

	auto x0 = StateParams(
		-5173447.1334, 4100505.6444, 0,
		-368.2895847, -575.7154581, 7689.0963483,
		0.018,
		JD(42162.430574375001015) + time::JD1899,
		1);
	auto x1 = StateParams(
		-3369338.2137, 5689531.5193, 0.0,
		-555.1752958, -393.8641492, 7680.8692571,
		0.018198167,
		JD(DateTime(2015, 6, 8, 11, 34, 11, 802)),
		1);
	auto x2 = StateParams(
		-5173447.1334, 4100505.6444, 0.0,
		-368.2895847, -575.7154581, 7689.0963483,
		0.0180,
		JD(DateTime(2015, 6, 7, 10, 20, 1, 626)),
		1);
	auto pGravity{ std::make_unique<PZ90>() };
	auto pAtmosphere{ std::make_unique<StaticAtmosphere81>(pGravity->R(), pGravity->Fl()) };
	auto ball = Ballistic(std::move(pGravity), std::move(pAtmosphere), 16);
	auto pBall = &ball;
	auto printFile = [pBall]() {
		auto fout = std::ofstream("ball test mma 16_1.txt");
		fout << std::setprecision(16);
		for (auto& x : pBall->Trajectory())
			fout << x.second.ToDateTime() << "; " << x.first << std::endl;
		fout.close();
	};
	try {
		ball.Run(x2, x2.T + 1.0);
		std::cout << "Succesfully calculated!\n";
		printFile();
	}
	catch (std::exception& ex)
	{
		std::cout << "An error occured! " << ex.what() << "\n";
		printFile();
	}
}
