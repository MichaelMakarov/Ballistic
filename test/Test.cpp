#include "PZ90.h"
#include "Ballistic.h"
#include "Atmosphere1981.h"
//#include "Atmosphere2004.h"
#include "SolarModel.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <future>
#include <AdamsIntegrator.h>
#include <RungeKuttaIntegrator.h>
#include <EGM96.h>
#include <GeoPotential.h>
#include "LowOrbModel.h"
#include "AdjustTrajectory.h"

using namespace ball;

void test_conversions();
void test_geopotential();
void test_atmosphere();
void test_ballistic();
void test_trajectoryproxy();

int main()
{
	//test_conversions();
	//test_geopotential();
	test_atmosphere();
	//test_ballistic();
	//test_trajectoryproxy();
	return 0;
}

void test_conversions()
{
	std::cout << "\n...Test conversions...\n";

	auto xyzPosition{ general::math::Vec3(-4688980.289, -11060428.914, 238914.750) };
	std::cout << "Initial position in orthogonal: " << xyzPosition << std::endl;
	auto rblPosition = ort_to_sph(xyzPosition);
	std::cout << "Converted to spherical: " << rblPosition << std::endl;
	auto xyzCheck = sph_to_ort(rblPosition);
	std::cout << "Converted to orthogonal: " << xyzCheck << std::endl;
}

void test_geopotential()
{
	using namespace general;
	std::cout << "\n...Test geopotential...\n";

	std::cout << std::setprecision(16);
	auto gp{ GeoPotential(EGM96::R(), EGM96::Mu(), EGM96(), 16) };
	double delta = math::deg_to_rad(5);
	double latitude{ 0 };
	double longitude{ math::deg_to_rad(349.45) };
	double radius{ 6378137.000 + 5.63755000000E+06 };
	for (size_t i = 0; i < 60; ++i)
	{
		auto position{ math::Vec3(radius, latitude, longitude) };
		std::cout << "U(" << position << ") = " << gp(position) << std::endl;
		latitude += 0.1 * delta;
		longitude += delta;
	}
}

void test_atmosphere()
{
	using namespace general;
	std::cout << "\n...Test atmosphere...\n";

	
}

void test_ballistic1()
{
	using namespace general::time;
	std::cout << "\n...Test ballistic 1...\n";

	std::vector<State> list {
		State(
			-5173447.1334, 4100505.6444, 0,
			-368.2895847, -575.7154581, 7689.0963483,
			0.018,
			JD(42162.430574375001015) + general::time::JD1899,
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
	auto index{ 1 };
	auto calculate = [](const State& x, const double dt, const size_t index)
	{
		//std::this_thread::sleep_for(std::chrono::seconds(10));
		auto ball = make_forecast(std::make_unique<StatAtmModel>(
			PZ90::Mu(), PZ90::R(), PZ90::W(), PZ90::Fl(), PZ90(), 16));
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
			ball.run(
				x, x.T + dt, 
				RKIntegrator<general::math::PV>(), 
				AdamsIntegrator<general::math::PV>());
			std::cout << index << " Successfully calculated!\n";
			printfile(filename);
		}
		catch (std::exception& ex) {
			std::cout << index << " An error occured! " << ex.what() << "\n";
			printfile(filename);
		}
	};
	auto tasks{ std::vector<std::future<void>>(list.size()) };
	for (size_t i = 0; i < list.size(); ++i) {
		tasks[i] = std::async(std::launch::async, calculate, list[i], 1.0, index++);
	}
	
}

void test_ballistic2()
{
	using namespace general::time;
	std::cout << "\n...Test ballistic 2...\n";
	std::cout << std::setprecision(15);

	auto x0{ State(
		1655697.1365, -3861686.2435, 5064524.7827,
		3705.8211976, -4786.9125982, -4812.0157646,
		0.018059454,
		JD(42187, 0.321120925924333) + JD1899,
		1) };
	auto ball = make_forecast(std::make_unique<StatAtmModel>(
		EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(), EGM96(), 50));
	std::cout << "initial point:\n";
	std::cout << "T: " << x0.T.to_datetime() << "; x: " <<
		x0.Vec << "; s = " << x0.Sb << "; loop = " << x0.Loop << std::endl;
	try {
		ball.run(
			x0, x0.T + 0.5,
			RKIntegrator<general::math::PV>(),
			AdamsIntegrator<general::math::PV>());
		State x;
		JD list[] { x0.T, x0.T, x0.T, x0.T };
		int deltas[] { 17, 300, 3600, 86330 };
		for (size_t i = 0; i < sizeof(deltas) / sizeof(deltas[0]); ++i) {
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
	test_ballistic1();
	test_ballistic2();
}

void test_trajectoryproxy1()
{
	std::cout << "\n...Trajectory proxy test...\n";
	using namespace general::time;
	using namespace general::math;

	std::cout << std::setprecision(15);

	auto read_measurements_from_txt = [](std::string& filepath) -> std::list<std::pair<PV, JD>> {
		auto is = std::ifstream(filepath);
		if (!is.is_open())
			throw std::invalid_argument("Invalid istream!");
		auto measurements = std::list<std::pair<PV, JD>>();
		unsigned short date[3]{}, time[3]{};
		PV coords;
		char buf{ ' ' };
		while (!is.eof()) {
			is >> date[0] >> date[1] >> date[2] >>
				time[0] >> time[1] >> time[2] >>
				coords.pos.x >> coords.pos.Y >> coords.pos.Z >>
				coords.Vel.x >> coords.Vel.Y >> coords.Vel.Z;
			coords *= 1e3;
			measurements.push_back(std::make_pair(coords, JD(DateTime(date[2], date[1], date[0], time[0], time[1], time[2]))));
			while (!is.eof() && buf != '\n')
				is.read(&buf, 1);
		}
		is.close();
		return measurements;
	};

	auto x0{ State(
		-5173447.1334, 4100505.6444, 0,
		-368.2895847, -575.7154581, 7689.0963483,
		0.018,
		JD(42162 + static_cast<size_t>(JD1899), 0.430574375001015),
		1
	) };
	std::string directory = "D:\\User\\Desktop\\disser\\progs\\Ballistic\\resources\\ExampleDataTxt\\";
	auto filepath = directory + "meas1.txt";
	auto list = read_measurements_from_txt(filepath);
	auto measurements{ std::vector<std::pair<PV, JD>>(list.size()) };
	size_t index{ 0 };
	for (auto& v : list)
		measurements[index++] = v;
	std::function<std::unique_ptr<StatAtmModel>()> make_statmodel = 
		[]() -> std::unique_ptr<StatAtmModel> {
		return std::make_unique<StatAtmModel>(PZ90::Mu(), PZ90::R(), PZ90::W(), PZ90::Fl(), PZ90(), 16);
	};
	std::cout << "initial x0: " << x0 << std::endl;
	auto iter = refine_initpoint(x0, measurements, make_statmodel);
	std::cout << "refined x0: " << x0 << std::endl;
	//OUT_LINE(std::cout, x0);
}

void test_trajectoryproxy2()
{
	using namespace general::math;
	using namespace general::time;

	struct Difference
	{
		DateTime T;
		double R{ 0 }, V{ 0 };
	};

	auto read_measurements_from_txt = [](const std::string& filepath) -> std::list<std::pair<PV, JD>> {
		auto reader = std::ifstream(filepath);
		if (!reader.is_open()) 
			throw std::runtime_error("Failed to open file: " + filepath);
		std::string date, time;
		PV vec;
		DateTime dt;
		char buf{ 0 };
		auto list = std::list<std::pair<PV, JD>>();

		while (!reader.eof()) {
			reader >> date >> time >>
				vec.pos.x >> vec.pos.Y >> vec.pos.Z >>
				vec.Vel.x >> vec.Vel.Y >> vec.Vel.Z;
			if (!datetime_from_str(date + " " + time, dt, "y-M-d h:m:s"))
				throw std::runtime_error("Failed to parse datetime from " + date + " " + time);
			list.push_back(std::make_pair(vec, JD(dt).add_seconds(-18)));
			while (buf != '\n' && !reader.eof()) reader.read(&buf, 1);
		}
		reader.close();
		return list;
	};
	auto save_to_file = [](const std::string& filepath, const State& x0, const std::vector<Difference>& arr) {
		auto writer = std::ofstream(filepath);
		if (!writer.is_open()) 
			throw std::runtime_error("Failed to open file: " + filepath);
		writer << x0 << "\n";
		for (const auto& p : arr) writer << p.T << "; dR = " << p.R << "; dV = " << p.V << "\n";
		writer.close();
	};
	try {
		const std::string directory = "D:\\User\\Desktop\\disser\\progs\\Ballistic\\resources\\ExampleDataTxt\\";
		const std::string filename = "25_09_2019.txt";
		const std::string filepath = directory + filename;
		auto list = read_measurements_from_txt(filepath);
		auto& measurements = list;// filter_measurements(list);
		JD tk{ (--measurements.cend())->second };
		auto x0{ State(measurements.cbegin()->first, 0, measurements.cbegin()->second, 1) };
		const size_t harmonics{ 32 };
		std::function<std::unique_ptr<StatAtmModel>()> make_statmodel =
			[harmonics]() -> std::unique_ptr<StatAtmModel> {
			return std::make_unique<StatAtmModel>(EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(), EGM96(), harmonics);
		};
		std::cout << "initial x0: " << x0 << std::endl;
		auto iter = refine_initpoint(x0, measurements, make_statmodel);
		std::cout << "refined x0: " << x0 << std::endl;

		auto ball = make_forecast(make_statmodel());
		ball.run(x0, tk, RKIntegrator<PV>(), AdamsIntegrator<PV>());
		auto diversities = std::vector<Difference>(measurements.size());
		State point;
		size_t index = 0;
		PV diff;
		for (const auto& m : measurements) {
			ball.get_point(m.second, point);
			diff = m.first - point.Vec;
			diversities[index++] = Difference{ m.second.to_datetime(), diff.pos.length(), diff.Vel.length() };
		}
		save_to_file(std::to_string(harmonics) + " harm data of " + filename, x0, diversities);
	}
	catch (std::exception& ex) {
		std::cout << ex.what() << std::endl;
	}
}

void test_trajectoryproxy()
{
	//test_trajectoryproxy1();
	test_trajectoryproxy2();
}
