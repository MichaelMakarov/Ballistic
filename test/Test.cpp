#include "PZ90.h"
#include "Ballistic.h"
#include "StaticAtmosphere.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <future>
#include <AdamsIntegrator.h>
#include <RungeKuttaIntegrator.h>
#include <EGM96.h>
#include <GeoPotential.h>
#include "PAForecast.h"

using namespace ball;
using namespace space;

void test_conversions();
void test_geopotential();
void test_atmosphere();
void test_ballistic();

int main()
{
	test_conversions();
	test_geopotential();
	test_atmosphere();
	test_ballistic();
	return 0;
}

void test_conversions()
{
	std::cout << "\n...Test conversions...\n";

	auto xyzPosition{ XYZ(-4688980.289, -11060428.914, 238914.750) };
	std::cout << "Initial position in orthogonal: " << xyzPosition << std::endl;
	auto rblPosition = CS_ortho_to_spher(xyzPosition);
	std::cout << "Converted to spherical: " << rblPosition << std::endl;
	auto xyzCheck = CS_spher_to_ortho(rblPosition);
	std::cout << "Converted to orthogonal: " << xyzCheck << std::endl;

}

void test_geopotential()
{
	std::cout << "\n...Test geopotential...\n";

	std::cout << std::setprecision(16);
	auto gp{ GeoPotential(std::move(std::unique_ptr<EGM96>{ new EGM96() }), 16) };
	double delta = ball::math::deg_to_rad(5);
	double latitude{ 0 };
	double longitude{ ball::math::deg_to_rad(349.45) };
	double radius{ 6378137.000 + 5.63755000000E+06 };
	for (size_t i = 0; i < 60; ++i)
	{
		auto position{ RBL(radius, latitude, longitude) };
		std::cout << "U(" << position << ") = " << gp(position) << std::endl;
		latitude += 0.1 * delta;
		longitude += delta;
	}
}

void test_atmosphere()
{
	std::cout << "\n...Test atmosphere...\n";

	EGM96 gm;
	auto atmosphere{ StaticAtmosphere81(gm.R(), gm.Fl()) };
	double delta = 1000.0;
	double height{ 10 };
	auto rblPosition{ RBL(gm.R() + height, 0, ball::math::deg_to_rad(349.45)) };
	std::cout << "Position: " << rblPosition << std::endl;
	for (size_t i = 0; i < 60; ++i)
	{
		std::cout << "h = " << height << "; rho = " <<
			atmosphere.density(CS_spher_to_ortho(rblPosition), ball::time::JD2000) << std::endl;
		height += delta;
		rblPosition.R += delta;
	}
}

void TestBallistic1()
{
	std::cout << "\n...Test ballistic 1...\n";

	std::vector<State> list {
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
	auto pForecast{ std::make_shared<PAForecast<StaticAtmosphere81>>(pGravity, 16, pAtmosphere) };
	auto index{ 1 };
	auto calculate = [pForecast](const State& x, const double dt, const size_t index)
	{
		//std::this_thread::sleep_for(std::chrono::seconds(10));
		auto ball = Ballistic<PAForecast<StaticAtmosphere81>>(pForecast);
		auto pBall = &ball;
		auto printfile = [pBall](const std::string& filename) {
			auto fout = std::ofstream(filename);
			fout << std::setprecision(16);
			for (auto& x : pBall->trajectory())
				fout << x.second.to_datetime() << "; " << x.first << std::endl;
			fout.close();
		};
		auto filename = "ball test mma version " + std::to_string(index) + ".txt";
		try {
			ball.Run(x, x.T + dt);
			std::cout << index << " Successfully calculated!\n";
			printfile(filename);
		}
		catch (std::exception& ex) {
			std::cout << index << " An error occured! " << ex.what() << "\n";
			printfile(filename);
		}
	};
	auto tasks{ std::vector<std::future<void>>(list.size()) };
	for (size_t i = 0; i < list.size(); ++i)
	{
		tasks[i] = std::async(std::launch::async, calculate, list[i], 1.0, index++);
	}
	
}

void TestBallistic2()
{
	std::cout << "\n...Test ballistic 2...\n";
	std::cout << std::setprecision(15);

	auto x0{ State(
		1655697.1365, -3861686.2435, 5064524.7827,
		3705.8211976, -4786.9125982, -4812.0157646,
		0.018059454,
		JD(42187, 0.321120925924333) + JD1899,
		1) };
	auto pGravity{ std::make_shared<EGM96>() };
	auto pAtmosphere{ std::make_shared<StaticAtmosphere81>(pGravity->R(), pGravity->Fl()) };
	auto pForecast{ std::make_shared<PAForecast<StaticAtmosphere81>>(pGravity, 50, pAtmosphere) };
	auto ball = Ballistic<PAForecast<StaticAtmosphere81>>(pForecast);
	std::cout << "initial point:\n";
	std::cout << "T: " << x0.T.to_datetime() << "; x: " <<
		x0.Vec << "; s = " << x0.Sb << "; loop = " << x0.Loop << std::endl;
	try {
		ball.Run(x0, x0.T + 0.5);
		State x;
		JD list[] { x0.T, x0.T, x0.T, x0.T };
		int deltas[] { 17, 300, 3600, 86330 };
		for (size_t i = 0; i < sizeof(deltas) / sizeof(deltas[0]); ++i)
		{
			list[i].add_seconds(deltas[i]);
			if (!ball.get_point(list[i], x)) {
				std::cout << "Failed to calculate the point for time: " << list[i].to_datetime() << std::endl;
			}
			else {
				std::cout << "T: " << list[i].to_datetime() << "; x: " <<
					x.Vec << "; s = " << x.Sb << "; loop = " << x.Loop << std::endl;
			}
		}
		
	}
	catch (std::exception& ex) {
		std::cout << "An error occured during the calculation! " << ex.what() << std::endl;
	}
}

void test_ballistic()
{
	/*TestBallistic1();
	TestBallistic2();*/
	std::vector<double> v{ 1, 2, 3, 9 };
	std::allocator<double> a;
	const size_t n = 10;
	double* d = a.allocate(n);
	for (size_t i = 0; i < n; ++i) {
		a.construct(d + i, i + 1);
	}
	a.deallocate(d, n);
}
