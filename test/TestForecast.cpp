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
#include "Structures.h"
#include "LowOrbModel.h"
#include "EGM96.h"
#include "OrbitAdjustment.h"
#include "RungeKuttaIntegrator.h"
#include "AdamsIntegrator.h"
#include "EverhartIntegrator.h"
#include "DataProviders.h"
#include "serialization/Serialization.h"
#include <tuple>

using namespace general::time;
using namespace general::math;
using namespace ball;

#define NAMEOF(x) #x

struct Difference
{
	DateTime T;
	double R{ 0 }, V{ 0 };

	Difference() = default;
	Difference(const DateTime& t, const double r, const double v) : T{ t }, R{ r }, V{ v } {}
	Difference(const Difference& d) = default;
	Difference(Difference&& d) noexcept = default;
	Difference& operator = (const Difference& d) = default;
	Difference& operator = (Difference&& d) noexcept = default;
};


std::list<std::pair<Vec6, JD>> read_measurements_from_txt(const std::string& filepath);
std::list<std::pair<Vec6, JD>> read_measurements_from_txt(const std::string& filepath, const std::string& dtformat);
std::list<std::pair<Vec6, JD>> read_measurements_from_xml(const std::string& filepath);
void write_residuals(const std::string& filepath, const Params<>& x0, const std::vector<Difference>& list);
std::vector<Difference> calc_residuals(
	const Params<>& x0,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const std::vector<Params<>>& trajectory);
void write_trajectory(const std::vector<Params<>>& list, const std::string& filepath);

void run_evaluation_mdasm_depr(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics
);
void run_evaluation_mdasm(
	const SpaceWeatherData& data, 
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics);
void run_evaluation_mdasm_vars(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics);
void run_evaluation_mdasmvb(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics);

template<size_t n> 
Params<6 + n * n> make_params_with_variations(const Params<6>& params, const double(&vars)[n])
{
	size_t index{ 0 };
	Params<6 + n * n> res;
	res.loop = params.loop;
	res.T = params.T;
	res.sb = params.sb;
	for (index = 0; index < 6; ++index) res.vec[index] = params.vec[index];
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			if (j == i) res.vec[index] = vars[i];
			++index;
		}
	}
	return res;
}

std::string g_filename;

int main()
{
	std::cout << std::setprecision(10);
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
	g_filename = files[1];
	std::cout << "testing " << g_filename << std::endl;
	auto measurements = read_measurements_from_txt(directory + g_filename, "y-M-d h:m:s.f");
	
	run_evaluation_mdasm(data, measurements, harmonics);
	//run_evaluation_mdasmvb(data, measurements, harmonics);

	return 0;
}

void run_evaluation_mdasm_depr(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics)
{
	Stopwatch sw;
	Params<> x0{ measurements.cbegin()->first, measurements.cbegin()->second, 2e-6, 1 };
	auto tk = std::max_element(
		measurements.cbegin(), measurements.cend(),
		[](const auto& l, const auto& r) { return l.second < r.second; })->second;
	auto calc_trajectory = [harmonics, &data, tk](const Params<>& x0, const std::list<std::pair<Vec6, JD>>& list) {
		auto model = MDASM(
			EGM96(), harmonics,
			EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(),
			data.F10_7, data.F81, data.Kpsum / 8,
			x0.sb);
		Forecast<6> forecast;
		forecast.run(
			x0, tk,
			make_memberfunc(&MDASM::function, &model),
			rk_integrator<Vec6, JD>(), adams_integrator<Vec6, JD>());
		auto trajectory = std::vector<Params<>>(list.size());
		size_t index{ 0 };
		for (const auto& m : list) trajectory[index++] = forecast.get_point(m.second);
		return trajectory;
	};
	std::cout << "init point: " << x0 << std::endl;
	sw.start();
	std::cout << "iterations number " << refine_initpoint(x0, measurements, std::function(calc_trajectory)) << std::endl;
	sw.finish();
	std::cout << "calc point: " << x0 << std::endl;
	std::cout << sw.duration() << std::endl;
	sw.start();
	auto differences = calc_residuals(x0, measurements, calc_trajectory(x0, measurements));
	sw.finish();
	std::cout << sw.duration() << std::endl;
	write_residuals(std::to_string(harmonics) + " harm " + __func__ + " result " + g_filename, x0, differences);
}

template<class R, class T, class Inv> void integrate(
	const R& x0,
	const T& t0,
	const double step,
	R& xk,
	T& tk,
	const invoker<Inv, R, const R&, const T&>& func) {
	func(x0, t0);
}

void run_evaluation_mdasm(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics)
{
	Params<> x0{ measurements.cbegin()->first, measurements.cbegin()->second, 2e-6, 1 };
	auto tk = std::max_element(
		measurements.cbegin(), measurements.cend(),
		[](const auto& l, const auto& r) { return l.second < r.second; })->second;
	const double vars[]{ 25, 25, 25, 0.25, 0.25, 0.25, 0.00016 };

	auto calc_forecast = [harmonics, &data](
		const Params<>& x0, const JD& tk, 
		const std::list<std::pair<Vec6, JD>>& measurements) {
		auto model = MDASM(
			EGM96(), harmonics,
			EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(),
			data.F10_7, data.F81, data.Kpsum / 8,
			x0.sb);
		Forecast<6> forecast;
		forecast.run(
			x0, tk,
			make_memberfunc(&MDASM::function, &model),
			rk_integrator<Vec6, JD>(), adams_integrator<Vec6, JD>());
		auto points = std::vector<Params<>>(measurements.size());
		size_t index{ 0 };
		for (const auto& m : measurements) points[index++] = forecast.get_point(m.second);
		return points;
	};
	auto calc_traj = [&calc_forecast](
		const Params<>& x0, const JD& tk, 
		const std::list<std::pair<Vec6, JD>>& measurements,
		const double(&vars)[7]) {
		Params<> xlist[7]{ x0, x0, x0, x0, x0, x0, x0 };
		for (size_t i = 0; i < 6; ++i) xlist[i].vec[i] += vars[i];
		xlist[6].sb += vars[6];
		std::future<std::vector<Params<>>> futures[7];
		
		for (size_t i = 0; i < 7; ++i) futures[i] = std::async(std::launch::async, calc_forecast, xlist[i], tk, measurements);
		std::array<std::vector<Params<>>, 8> result;
		result[7] = calc_forecast(x0, tk, measurements);
		for (size_t i = 0; i < 7; ++i) result[i] = futures[i].get();
		return result;
	};
	auto correct_params = [](const Vector& dx0, Params<>& x0) {
		for (size_t i = 0; i < 6; ++i) x0.vec[i] += dx0[i];
		x0.sb += dx0[6];
	};
	auto fout = std::ofstream("log data.txt");
	Stopwatch sw;
	std::cout << "init point: " << x0 << std::endl;
	sw.start();
	std::cout << "iterations = " << evaluate_params<7>(x0, measurements, std::function(calc_traj), std::function(correct_params), vars, fout.rdbuf()) << std::endl;
	sw.finish();
	fout.close();
	std::cout << "calc point: " << x0 << std::endl;
	std::cout << sw.duration() << std::endl;;
	auto differences = calc_residuals(x0, measurements, calc_forecast(x0, tk, measurements));
	write_residuals(std::to_string(harmonics) + " harm " + __func__ + " result " + g_filename, x0, differences);

}

void run_evaluation_mdasm_vars(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics)
{
	Params<> x0{ measurements.cbegin()->first, measurements.cbegin()->second, 2e-6, 1 };
	const auto tk = std::max_element(
		measurements.cbegin(), measurements.cend(),
		[](const auto& l, const auto& r) { return l.second < r.second; })->second;
	const auto tn = x0.T;
	const double vars[7]{ 5, 5, 5, 0.05, 0.05, 0.05, 0.00016 };
	auto xv = make_params_with_variations(x0, vars);
	auto make_model = [&data, harmonics](const double sb) {
		return MDASM(
			EGM96(), harmonics,
			EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(),
			data.F10_7, data.F81, data.Kpsum / 8,
			sb);
	};
	auto calc_equation = [harmonics, &data, tn, tk, &measurements, &make_model](const Params<55>& x0, Matrix& mx, Vector& v) {
		auto model = make_model(x0.sb);
		Forecast<55> forecast;
		forecast.run(
			x0, tk,
			make_memberfunc(&MDASM::extfunction, &model),
			rk_integrator<Vec<55>, JD>(), adams_integrator<Vec<55>, JD>());
		size_t index{ 0 };
		for (const auto& m : measurements) {
			const auto& point = forecast.get_point(m.second);
			for (size_t i = 0; i < 3; ++i) {
				v[index] = m.first[i] - point.vec[i];
				for (size_t j = 0; j < 7; ++j)
					mx(index, j) = point.vec[6 + j * 7 + i] / x0.vec[6 + j * 7 + j];
				++index;
			}
		}
	};
	auto correct = [](const Vector& dx0, Params<55>& x0) {
		for (size_t i = 0; i < 6; ++i) x0.vec[i] += dx0[i];
		x0.sb += dx0[6];
	};

	Stopwatch sw;

	std::cout << "init point: " << xv << std::endl;
	sw.start();
	std::cout << "iterations = " << evaluate_params<55>(xv, std::function(calc_equation), std::function(correct), measurements.size() * 3) << std::endl;
	sw.finish();
	std::cout << "calc point: " << xv << std::endl;
	std::cout << sw.duration();
	for (size_t i = 0; i < 6; ++i) x0.vec[i] = xv.vec[i];
	x0.sb = xv.sb;
	sw.start();
	Forecast<6> forecast;
	auto model = make_model(x0.sb);
	forecast.run(x0, tk, make_memberfunc(&MDASM::function, &model), rk_integrator<Vec6, JD>(), adams_integrator<Vec6, JD>());
	auto trajectory = std::vector<Params<>>(measurements.size());
	size_t index{ 0 };
	for (const auto& m : measurements) trajectory[index++] = forecast.get_point(m.second);
	auto differences = calc_residuals(x0, measurements, trajectory);
	sw.finish();
	write_residuals(std::to_string(harmonics) + " harm " + __func__ + " result " + g_filename, x0, differences);
}

void run_evaluation_mdasmvb(
	const SpaceWeatherData& data,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const size_t harmonics)
{
	Params<> x0{ measurements.cbegin()->first, measurements.cbegin()->second, 2e-6, 1 };
	const size_t count{ 5 };
	auto ballparams = std::array<double, count>();
	const double vars[]{ 25, 25, 25, 0.25, 0.25, 0.25, 0.00016, 1e-6, 1e-8, 1e-8, 1e-8 };
	const size_t num = sizeof(vars) / sizeof(double);
	const auto tk = std::max_element(
		measurements.cbegin(), measurements.cend(),
		[](const auto& l, const auto& r) { return l.second < r.second; })->second;

	auto make_model = [harmonics, &data](const std::array<double, count>& sb) {
		return MDASMVb(
			EGM96(), harmonics,
			EGM96::Mu(), EGM96::R(), EGM96::W(), EGM96::Fl(),
			data.F10_7, data.F81, data.Kpsum / 8,
			sb);
	};
	auto calc_forecast = [&make_model](
		const Params<>& x0, const JD& tk,
		const std::list<std::pair<Vec6, JD>>& measurements,
		const std::array<double, count>& par) {
		auto model = make_model(par);
		Forecast<6> forecast;
		forecast.run(
			x0, tk,
			make_memberfunc(&MDASMVb::function, &model),
			rk_integrator<Vec6, JD>(), adams_integrator<Vec6, JD>());
		auto points = std::vector<Params<>>(measurements.size());
		size_t index{ 0 };
		for (const auto& m : measurements) points[index++] = forecast.get_point(m.second);
		return points;
	};
	auto calc_traj = [&calc_forecast, &ballparams](
		const Params<>& x0, const JD& tk, 
		const std::list<std::pair<Vec6, JD>>& measurements, 
		const double(&vars)[num]) {
		Params<> xlist[num]{ x0, x0, x0, x0, x0, x0, x0, x0, x0, x0, x0 };
		std::array<std::array<double, count>, num> params{ 
			ballparams, ballparams, ballparams, ballparams, ballparams, ballparams, 
			ballparams, ballparams, ballparams, ballparams, ballparams
		};
		for (size_t i = 0; i < 6; ++i) xlist[i].vec[i] += vars[i];
		std::get<0>(params[6]) += vars[6];
		std::get<1>(params[7]) += vars[7];
		std::get<2>(params[8]) += vars[8];
		std::get<3>(params[9]) += vars[9];
		std::get<4>(params[10]) += vars[10];
		std::future<std::vector<Params<>>> futures[num];
		
		for (size_t i = 0; i < num; ++i)
			futures[i] = std::async(std::launch::async, calc_forecast, xlist[i], tk, measurements, params[i]);
		std::array<std::vector<Params<>>, num + 1> result;
		result[num] = calc_forecast(x0, tk, measurements, ballparams);
		for (size_t i = 0; i < num; ++i) result[i] = futures[i].get();
		return result;
	};
	std::function<void(const Vector&, Params<>&)> correct_params =
		[&ballparams](const Vector& dx0, Params<>& x0) {
		for (size_t i = 0; i < 6; ++i) x0.vec[i] += dx0[i];
		x0.sb += dx0[6];
		std::get<0>(ballparams) += dx0[6];
		std::get<1>(ballparams) += dx0[7];
		std::get<2>(ballparams) += dx0[8];
		std::get<3>(ballparams) += dx0[9];
		std::get<4>(ballparams) += dx0[10];
	};

	Stopwatch sw;
	auto fout = std::ofstream("log oscul.txt");
	auto traj = calc_forecast(x0, tk, measurements, ballparams);
	for (const auto& p : traj) {
		double sidt = sidereal_time_mean(p.T);
		const auto oscul = oscul_from_ACS(GCS_to_ACS(slice<0, 2>(p.vec), sidt), GCS_to_ACS(slice<3, 5>(p.vec), sidt), EGM96::Mu());
		fout << rad_to_deg(oscul.ascendnode) << " " << rad_to_deg(oscul.latitudearg) << std::endl;
	}
	fout.close();
	return;

	std::cout << "init point: " << x0 << std::endl;
	sw.start();
	std::cout << "iterations = " << evaluate_params<11>(x0, measurements, calc_traj, correct_params, vars, nullptr) << std::endl;
	sw.finish();
	std::cout << "calc point: " << x0 << std::endl;
	std::cout << sw.duration();
	sw.start();
	auto differences = calc_residuals(x0, measurements, calc_forecast(x0, tk, measurements, ballparams));
	sw.finish();
	write_residuals(std::to_string(harmonics) + " harm " + __func__ + " result " + g_filename, x0, differences);
}

std::vector<Difference> calc_residuals(
	const Params<>& x0,
	const std::list<std::pair<Vec6, JD>>& measurements,
	const std::vector<Params<>>& trajectory)
{
	auto diversities{ std::vector<Difference>(measurements.size()) };
	size_t index{ 0 };
	for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter) {
		const auto& diff = iter->first - trajectory[index].vec;
		diversities[index++] = Difference{ 
			jd_to_datetime(iter->second),
			Vec3{ diff[0], diff[1], diff[2] }.length(),
			Vec3{ diff[3], diff[4], diff[5] }.length() 
		};
	}
	std::sort(diversities.begin(), diversities.end(), [](const Difference& l, const Difference& r) { return l.R < r.R; });
	std::cout << "mediane: " << diversities[diversities.size() / 2].R << std::endl;
	return diversities;
}

void write_residuals(const std::string& filepath, const Params<>& x0, const std::vector<Difference>& list)
{
	auto fout{ std::ofstream(filepath) };
	if (!fout.is_open())
		throw std::runtime_error("Failed to open file: " + filepath);
	fout << std::setprecision(6) << x0 << "\n";
	for (size_t i = 0; i < list.size(); ++i) {
		fout << list[i].T << " dR = ";
		fout.width(7);
		fout << std::left << list[i].R << " dV = ";
		fout.width(7);
		fout << std::left << list[i].V << "\n";
	}
	fout.close();
}

std::list<std::pair<Vec6, JD>> read_measurements_from_txt(const std::string& filepath)
{
	auto fin = std::ifstream(filepath);
	if (!fin.is_open())
		throw std::invalid_argument("Invalid istream!");
	auto measurements = std::list<std::pair<Vec6, JD>>();
	unsigned short date[3]{}, time[3]{};
	Vec6 coords;
	char buf{ ' ' };
	while (!fin.eof()) {
		fin >> date[0] >> date[1] >> date[2] >>
			time[0] >> time[1] >> time[2] >> coords;
		coords *= 1e3;
		measurements.push_back(std::make_pair(coords, JD(DateTime(date[2], date[1], date[0], time[0], time[1], time[2]))));
		while (!fin.eof() && buf != '\n')
			fin.read(&buf, 1);
	}
	fin.close();
	return measurements;
}

std::list<std::pair<Vec6, JD>> read_measurements_from_txt(const std::string& filepath, const std::string& dtformat)
{
	auto fin{ std::ifstream(filepath) };
	if (!fin.is_open())
		throw std::runtime_error("Failed to open file: " + filepath);
	std::string date, time;
	Vec6 vec;
	char buf{ 0 };
	auto list = std::list<std::pair<Vec6, JD>>();
	while (!fin.eof()) {
		fin >> date >> time >> vec;
		list.push_back(std::make_pair(vec, JD(datetime_from_str(date + " " + time, dtformat)).add_seconds(-18)));
		while (buf != '\n' && !fin.eof()) fin.read(&buf, 1);
	}
	fin.close();
	return list;
}

std::list<std::pair<Vec6, JD>> read_measurements_from_xml(const std::string& filepath)
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

	auto parse_nu = [](const Nu& nu) -> std::pair<Vec6, JD> {
		return std::make_pair(
			Vec6{
				std::strtod(nu.params[1].c_str(), nullptr),
				std::strtod(nu.params[2].c_str(), nullptr),
				std::strtod(nu.params[3].c_str(), nullptr),
				std::strtod(nu.params[4].c_str(), nullptr),
				std::strtod(nu.params[5].c_str(), nullptr),
				std::strtod(nu.params[6].c_str(), nullptr)
			},
			JD(general::time::datetime_from_str(nu.params[0]))
		);
	};

	std::list<Nu> list;
	std::function<IXmlSerializable&(std::list<Nu>&)> func = [](std::list<Nu>& list) -> IXmlSerializable& {
		list.push_back(Nu());
		return *(--list.end());
	};

	read_xml(func, list, filepath);
	std::list<std::pair<Vec6, JD>> data;
	for (const auto& nu : list) {
		data.push_back(parse_nu(nu));
	}
	return data;
}


void write_trajectory(const std::vector<Params<>>& list, const std::string& filepath)
{
	auto fout{ std::ofstream(filepath) };
	if (!fout.is_open()) throw std::runtime_error("Failed to create file: " + filepath);
	for (size_t i = 0; i < list.size(); ++i) {
		fout << list[i] << std::endl;
	}
	fout.close();
}