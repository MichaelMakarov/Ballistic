#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <utility>
#include <list>
#include <vector>
#include <functional>
#include <filesystem>
#include "general/Times.h"
#include "general/Geometry.h"
#include "Structures.h"
#include "LowOrbModel.h"
#include "EGM96.h"
#include "AdjustTrajectory.h"
#include "RungeKuttaIntegrator.h"
#include "AdamsIntegrator.h"
#include "DataProviders.h"

using namespace general::time;
using namespace general::math;
using namespace ball;

struct Difference
{
	DateTime T;
	double R{ 0 }, V{ 0 };
};

template<class Model> void adjust_trajectory(
	State& x0, const std::list<std::pair<PV, JD>>& measurements,
	std::function<std::unique_ptr<Model>()>& make_model)
{
	std::cout << "init point: " << x0 << std::endl;
	std::cout << "iterations number " << refine_initpoint(x0, measurements, make_model)	<< std::endl;
	std::cout << "calc point: " << x0 << std::endl;
}
std::list<std::pair<PV, JD>> read_measurements(const std::string& filepath);
std::list<std::pair<PV, JD>> read_measurements(const std::string& filepath, const std::string& dtformat);

template<class Iterator>
void write_differences(const std::string& filepath, const State& x0, Iterator first, Iterator last)
{
	auto fout{ std::ofstream(filepath) };
	if (!fout.is_open())
		throw std::runtime_error("Failed to open file: " + filepath);
	fout << x0 << "\n";
	for (auto iter = first; iter != last; ++iter)
		fout << iter->T << "; dR = " << iter->R << "; dV = " << iter->V << "\n";
	fout.close();
}
template<class Model> std::vector<Difference> calc_differences(
	const State& x0, const std::list<std::pair<PV, JD>>& measurements,
	std::function<std::unique_ptr<Model>()>& make_model)
{
	auto diversities{ std::vector<Difference>(measurements.size()) };
	auto forecast{ make_forecast(make_model()) };
	forecast.run(x0, (--measurements.end())->second, RKIntegrator<PV>(), AdamsIntegrator<PV>());
	size_t index{ 0 };
	for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter) {
		const auto& diff = iter->first - forecast.get_point(iter->second).Vec;
		diversities[index++] = Difference{ iter->second.to_datetime(), diff.Pos.length(), diff.Vel.length() };
	}
	return diversities;
}

int main()
{
	std::cout << std::setprecision(15);

	auto path = std::filesystem::current_path().parent_path().parent_path().parent_path();
	std::string celestrackfile = path.string() + "/resources/celestrack/spaceweather.txt";

	const size_t harmonics{ 32 };

	auto spaceweather = load_spaceweather_data(celestrackfile);
	const auto& data = spaceweather[Date(2019, 7, 20)];

	std::function<std::unique_ptr<StatAtmModel>()> make_statmodel = [harmonics]() {
		return std::make_unique<StatAtmModel>(EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(), EGM96(), harmonics);
	};
	std::function<std::unique_ptr<DynAtmModel>()> make_dynmodel = [harmonics, data]() {
		return std::make_unique<DynAtmModel>(EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(), data.F10_7, data.F81, data.Kpsum / 8, EGM96(), harmonics);
	};
	
	std::string directory = path.string() + "/resources/exampledata/";
	std::string filename;
	
	filename = "20_07_2019.txt";
	std::cout << "testing " << filename << std::endl;
	auto measurements = read_measurements(directory + filename, "y-M-d h:m:s.f");

	auto x0{ State(measurements.begin()->first, 0, measurements.cbegin()->second, 1) };
	adjust_trajectory(x0, measurements, make_statmodel);
	auto differences = calc_differences(x0, measurements, make_statmodel);
	write_differences(std::to_string(harmonics) + " harm stat atm test data diff " + filename, x0, differences.cbegin(), differences.cend());

	x0 = State(measurements.begin()->first, 0, measurements.cbegin()->second, 1);
	adjust_trajectory(x0, measurements, make_dynmodel);
	differences = calc_differences(x0, measurements, make_dynmodel);
	write_differences(std::to_string(harmonics) + " harm dynm atm test data diff " + filename, x0, differences.cbegin(), differences.cend());


	return 0;
}


std::list<std::pair<PV, JD>> read_measurements(const std::string& filepath)
{
	auto fin = std::ifstream(filepath);
	if (!fin.is_open())
		throw std::invalid_argument("Invalid istream!");
	auto measurements = std::list<std::pair<PV, JD>>();
	unsigned short date[3]{}, time[3]{};
	PV coords;
	char buf{ ' ' };
	while (!fin.eof()) {
		fin >> date[0] >> date[1] >> date[2] >>
			time[0] >> time[1] >> time[2] >>
			coords.Pos.X >> coords.Pos.Y >> coords.Pos.Z >>
			coords.Vel.X >> coords.Vel.Y >> coords.Vel.Z;
		coords *= 1e3;
		measurements.push_back(std::make_pair(coords, JD(DateTime(date[2], date[1], date[0], time[0], time[1], time[2]))));
		while (!fin.eof() && buf != '\n')
			fin.read(&buf, 1);
	}
	fin.close();
	return measurements;
}

std::list<std::pair<PV, JD>> read_measurements(const std::string& filepath, const std::string& dtformat)
{
	auto fin{ std::ifstream(filepath) };
	if (!fin.is_open())
		throw std::runtime_error("Failed to open file: " + filepath);
	std::string date, time;
	PV vec;
	DateTime dt;
	char buf{ 0 };
	auto list = std::list<std::pair<PV, JD>>();
	while (!fin.eof()) {
		fin >> date >> time >>
			vec.Pos.X >> vec.Pos.Y >> vec.Pos.Z >>
			vec.Vel.X >> vec.Vel.Y >> vec.Vel.Z;
		if (!try_parse(date + " " + time, dt, dtformat))
			throw std::runtime_error("Failed to parse datetime from " + date + " " + time);
		list.push_back(std::make_pair(vec, JD(dt).add_seconds(-18)));
		while (buf != '\n' && !fin.eof()) fin.read(&buf, 1);
	}
	fin.close();
	return list;
}
