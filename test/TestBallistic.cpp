#include "PZ90.h"
#include "Ballistic.h"
#include "StaticAtmosphere.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <future>

using namespace ball;
using namespace space;
using namespace geometry;

void TestIntegrator();
void TestBallistic();

int main()
{
	TestBallistic();
	return 0;
}

void TestIntegrator()
{
	std::cout << "\n...Test integrators...\n";

	std::cout << std::setprecision(16);
	/*types::PV x0(
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
	}*/

}

void TestBallistic1()
{
	std::cout << "\n...Test ballistic 1...\n";

	auto list =	{
		State(
			-5173447.1334, 4100505.6444, 0,
			-368.2895847, -575.7154581, 7689.0963483,
			0.018,
			JD(42162.430574375001015) + time::JD1899,
			1),
		State(
			-3369338.2137, 5689531.5193, 0.0,
			-555.1752958, -393.8641492, 7680.8692571,
			0.018198167,
			JD(DateTime(2015, 6, 8, 11, 34, 11, 802)),
			1),
		State(
			-5173447.1334, 4100505.6444, 0.0,
			-368.2895847, -575.7154581, 7689.0963483,
			0.0180,
			JD(DateTime(2015, 6, 7, 10, 20, 1, 626)),
			1)
	};
	auto pGravity{ std::make_shared<PZ90>() };
	auto pAtmosphere{ std::make_shared<StaticAtmosphere81>(pGravity->R(), pGravity->Fl()) };
	auto index{ 1 };
	auto calculate = [pGravity, pAtmosphere](const State& x, const double dt, const size_t index)
	{
		auto ball = Ballistic(pGravity, pAtmosphere, 16);
		auto pBall = &ball;
		auto printFile = [pBall](const std::string& filename) {
			auto fout = std::ofstream(filename);
			fout << std::setprecision(16);
			for (auto& x : pBall->Trajectory())
				fout << x.second.ToDateTime() << "; " << x.first << std::endl;
			fout.close();
		};
		auto filename = "ball test mma version " + std::to_string(index) + ".txt";
		try {
			ball.Run(x, x.T + dt);
			std::cout << index << " Successfully calculated!\n";
			printFile(filename);
		}
		catch (std::exception& ex)
		{
			std::cout << index << " An error occured! " << ex.what() << "\n";
			printFile(filename);
		}
	};
	/*auto ball = Ballistic(pGravity, pAtmosphere, 16);
	auto pBall = &ball;
	auto printFile = [pBall](const std::string& filename) {
		auto fout = std::ofstream(filename);
		fout << std::setprecision(16);
		for (auto& x : pBall->Trajectory())
			fout << x.second.ToDateTime() << "; " << x.first << std::endl;
		fout.close();
	};*/
	for (auto& x : list)
	{
		auto future = std::async(std::launch::async, calculate, x, 1.0, index++);
		//future.get();
		/*std::cout << index << std::ends;
		auto filename = "ball test mma version " + std::to_string(index++) + ".txt";
		try {
			ball.Run(x, x.T + 1.0);
			std::cout << " Successfully calculated!\n";
			printFile(filename);
		}
		catch (std::exception& ex)
		{
			std::cout << "An error occured !" << ex.what() << "\n";
			printFile(filename);
		}
		std::cout << std::endl;*/
	}
	
}

void TestBallistic2()
{
	std::cout << "\n...Test ballistic 2...\n";
	std::shared_ptr<State> ptr;
	std::cout << std::setprecision(15);

	auto x0{ State(
		1655697.1365, -3861686.2435, 5064524.7827,
		3705.8211976, -4786.9125982, -4812.0157646,
		0.018059454,
		JD(42187, 0.321120925924333) + JD1899,
		1) };
	auto pGravity{ std::make_shared<PZ90>() };
	auto pAtmosphere{ std::make_shared<StaticAtmosphere81>(pGravity->R(), pGravity->Fl()) };
	auto ball = Ballistic(pGravity, pAtmosphere, 16);
	std::cout << "initial point:\n";
	std::cout << "T: " << x0.T.ToDateTime() << "; x: " <<
		x0.Vec << "; s = " << x0.Sb << "; loop = " << x0.Loop << std::endl;
	try {
		ball.Run(x0, x0.T + 0.5);
		State x;
		JD list[] { x0.T, x0.T, x0.T, x0.T };
		int deltas[] { 17, 300, 3600, 86330 };
		for (size_t i = 0; i < sizeof(deltas) / sizeof(deltas[0]); ++i)
		{
			list[i].AddSeconds(deltas[i]);
			if (!ball.GetPoint(list[i], x))
			{
				std::cout << "Failed to calculate the point for time: " << list[i].ToDateTime() << std::endl;
			}
			else {
				std::cout << "T: " << list[i].ToDateTime() << "; x: " <<
					x.Vec << "; s = " << x.Sb << "; loop = " << x.Loop << std::endl;
			}
		}
		
	}
	catch (std::exception& ex)
	{
		std::cout << "An error occured during the calculation!\n";
	}
}

void TestBallistic()
{
	TestBallistic1();
	TestBallistic2();
}
