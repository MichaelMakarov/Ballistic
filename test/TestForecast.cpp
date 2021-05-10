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
#include "EverhartIntegrator.h"
#include "DataProviders.h"
#include "serialization/Serialization.h"

using namespace general::time;
using namespace general::math;
using namespace ball;

struct Difference
{
	DateTime T;
	double R{ 0 }, V{ 0 };
};

std::list<std::pair<PV, JD>> read_measurements_from_txt(const std::string& filepath);
std::list<std::pair<PV, JD>> read_measurements_from_txt(const std::string& filepath, const std::string& dtformat);
std::list<std::pair<PV, JD>> read_measurements_from_xml(const std::string& filepath);
void write_residuals(const std::string& filepath, const State& x0, const std::vector<Difference>& list);
std::vector<Difference> calc_residuals(
	const State& x0,
	const std::list<std::pair<PV, JD>>& measurements,
	const std::function<std::vector<State>(const State&, const std::list<std::pair<PV, JD>>&)>& calc_traj);

int main()
{
	std::cout << std::setprecision(15);

	auto path = std::filesystem::current_path().parent_path().parent_path().parent_path();
	std::string celestrackfile = path.string() + "/resources/celestrack/spaceweather.txt";
	std::string directory = path.string() + "/resources/exampledata/";

	const size_t harmonics{ 32 };

	auto spaceweather = load_spaceweather_data(celestrackfile);
	const auto& data = spaceweather[Date(2019, 7, 20)];
	EGM96("resources/gravity/egm96.txt");
	Stopwatch sw;
	
	std::string files[]{
		"08_06_2018.txt",
		"18_11_2018.txt",
		"25_09_2019.txt",
		"20_07_2019.txt"
	};
	std::string filename = files[1];
	std::cout << "testing " << filename << std::endl;
	auto measurements = read_measurements_from_txt(directory + filename, "y-M-d h:m:s.f");
	auto x0{ State(measurements.cbegin()->first, 2e-6, measurements.cbegin()->second, 1) };
	JD tk = x0.T;
	for (const auto& m : measurements) {
		if (tk < m.second) tk = m.second;
	}
	std::function<std::vector<State>(const State&, const std::list<std::pair<PV, JD>>&)> calc_trajectory =
		[harmonics, data, tk](const State& x0, const std::list<std::pair<PV, JD>>& list) {
		auto forecast = make_forecast<DynAtmModel>(
			std::make_unique<DynAtmModel>(
				EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(),
				data.F10_7, data.F81, data.Kpsum / 8,
				EGM96(), harmonics));
		forecast.run(x0, tk, RKIntegrator<PV, JD>(), AdamsIntegrator<PV, JD>());
		auto traj = std::vector<State>(list.size());
		size_t index{ 0 };
		for (const auto& m : list) {
			traj[index++] = forecast.get_point(m.second);
		}
		return traj;
	};

	std::cout << "init point: " << x0 << std::endl;
	sw.start();
	std::cout << "iterations number " << refine_initpoint(x0, measurements, calc_trajectory) << std::endl;
	sw.finish();
	std::cout << "calc point: " << x0 << std::endl;
	std::cout << sw.duration() << std::endl;
	sw.start();
	auto differences = calc_residuals(x0, measurements, calc_trajectory);
	sw.finish();
	std::cout << sw.duration() << std::endl;
	write_residuals(std::to_string(harmonics) + " harm dynm atm test data diff " + filename, x0, differences);

	return 0;
}

std::vector<Difference> calc_residuals(
	const State& x0,
	const std::list<std::pair<PV, JD>>& measurements,
	const std::function<std::vector<State>(const State&, const std::list<std::pair<PV, JD>>&)>& calc_traj)
{
	auto diversities{ std::vector<Difference>(measurements.size()) };
	auto trajectory = calc_traj(x0, measurements);
	size_t index{ 0 };
	for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter) {
		const auto& diff = iter->first - trajectory[index].Vec;
		diversities[index++] = Difference{ iter->second.to_datetime(), diff.pos.length(), diff.vel.length() };
	}
	return diversities;
}

void write_residuals(const std::string& filepath, const State& x0, const std::vector<Difference>& list)
{
	auto fout{ std::ofstream(filepath) };
	if (!fout.is_open())
		throw std::runtime_error("Failed to open file: " + filepath);
	fout << x0 << "\n";
	for (size_t i = 0; i < list.size(); ++i)
		fout << list[i].T << "; dR = " << list[i].R << "; dV = " << list[i].V << "\n";
	fout.close();
}

std::list<std::pair<PV, JD>> read_measurements_from_txt(const std::string& filepath)
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
			coords.pos.x() >> coords.pos.y() >> coords.pos.z() >>
			coords.vel.x() >> coords.vel.y() >> coords.vel.z();
		coords *= 1e3;
		measurements.push_back(std::make_pair(coords, JD(DateTime(date[2], date[1], date[0], time[0], time[1], time[2]))));
		while (!fin.eof() && buf != '\n')
			fin.read(&buf, 1);
	}
	fin.close();
	return measurements;
}

std::list<std::pair<PV, JD>> read_measurements_from_txt(const std::string& filepath, const std::string& dtformat)
{
	auto fin{ std::ifstream(filepath) };
	if (!fin.is_open())
		throw std::runtime_error("Failed to open file: " + filepath);
	std::string date, time;
	PV vec;
	char buf{ 0 };
	auto list = std::list<std::pair<PV, JD>>();
	while (!fin.eof()) {
		fin >> date >> time >>
			vec.pos.x() >> vec.pos.y() >> vec.pos.z() >>
			vec.vel.x() >> vec.vel.y() >> vec.vel.z();
		list.push_back(std::make_pair(vec, JD(datetime_from_str(date + " " + time, dtformat)).add_seconds(-18)));
		while (buf != '\n' && !fin.eof()) fin.read(&buf, 1);
	}
	fin.close();
	return list;
}

std::list<std::pair<PV, JD>> read_measurements_from_xml(const std::string& filepath)
{
	using namespace serialization;

	struct Nu : public IXmlSerializable
	{
		std::string params[7];

		std::unique_ptr<XmlNode> serialize() const override { return nullptr; }
		void deserialize(const XmlNode& node) override
		{
			size_t index{ 0 };
			for (const auto& pNode : node.nodes) {
				params[index++] = pNode->value;
			}
		}
	};

	auto parse_nu = [](const Nu& nu) -> std::pair<PV, JD> {
		return std::make_pair(
			PV(
				std::strtod(nu.params[1].c_str(), nullptr),
				std::strtod(nu.params[2].c_str(), nullptr),
				std::strtod(nu.params[3].c_str(), nullptr),
				std::strtod(nu.params[4].c_str(), nullptr),
				std::strtod(nu.params[5].c_str(), nullptr),
				std::strtod(nu.params[6].c_str(), nullptr)
			),
			JD(general::time::datetime_from_str(nu.params[0]))
		);
	};

	std::list<Nu> list;
	std::function<IXmlSerializable&(std::list<Nu>&)> func = [](std::list<Nu>& list) -> IXmlSerializable& {
		list.push_back(Nu());
		return *(--list.end());
	};

	read_xml(func, list, filepath);
	std::list<std::pair<PV, JD>> data;
	for (const auto& nu : list) {
		data.push_back(parse_nu(nu));
	}
	return data;
}