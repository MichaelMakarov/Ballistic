#include "AdjustTrajectory.h"
#include "general/Polynomial.h"

namespace ball
{
	std::vector<general::math::Vector> create_array_of_vectors(const size_t count, const size_t size)
	{
		auto result{ std::vector<general::math::Vector>(count) };
		for (auto& vec : result) vec.resize(size);
		return result;
	}

	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		std::list<std::pair<general::math::PV, general::time::JD>>& measurements,
		const double dt,
		const size_t polydegree)
	{
		using namespace general::math;
		using namespace general::time;
		if (measurements.size() < 2) return measurements;
		measurements.sort([](std::pair<PV, JD>& f, std::pair<PV, JD>& s) { return f.second < s.second; });
		const size_t vars{ 6 };
		const size_t degree{ polydegree };
		const size_t sko{ 3 };
		JD	t0{ measurements.cbegin()->second },
			tk{ (--measurements.cend())->second };
		std::list<std::pair<PV, JD>> list, result;
		tk = t0 + 2 * dt;
		auto filter = [degree, vars, sko](const std::list<std::pair<PV, JD>>& list, std::list<std::pair<PV, JD>>& filtered) {
			if (list.size() < degree) {
				for (const auto& m : list) filtered.push_back(m);
			}
			else {
				auto buffer{ std::list<std::pair<PV, JD>>() };
				auto x{ Vector(list.size()) };
				auto Y = create_array_of_vectors(vars, list.size());
				size_t index{ 0 };
				for (auto iter = list.cbegin(); iter != list.cend(); ++iter) {
					x[index] = iter->second - list.cbegin()->second;
					for (size_t n = 0; n < vars; ++n) Y[n][index] = iter->first[n];
					++index;
				}
				auto results{ std::vector<std::vector<bool>>(vars) };
				for (size_t n = 0; n < vars; ++n) {
					results[n] = std::async(
						std::launch::async,
						[sko, degree](const Vector x, const Vector& y) -> std::vector<bool> {
							auto poly = create_polynom(x, y, degree);
							double currsko{ 0 };
							auto diff{ std::vector<double>(x.size()) };
							auto result{ std::vector<bool>(x.size()) };
							for (size_t k = 0; k < x.size(); ++k) {
								diff[k] = std::abs(y[k] - poly(x[k]));
								currsko += diff[k] * diff[k];
							}
							currsko = std::sqrt(currsko / (x.size() - 1));
							for (size_t k = 0; k < x.size(); ++k)
								result[k] = (diff[k] - sko * currsko) < 0.0;
							return result;
						}, x, Y[n]
							).get();
				}
				bool flag;
				index = 0;
				for (auto iter = list.cbegin(); iter != list.cend(); ++iter) {
					flag = true;
					for (size_t n = 0; n < vars; ++n) flag &= results[n][index];
					if (flag) buffer.push_back(*iter);
				}
				if (buffer.size() > 0) {
					auto iterator{ buffer.cbegin() };
					for (size_t k = 0; k < buffer.size() / 2; ++k)
						++iterator;
					filtered.push_back(*iterator);
				}
			}
		};

		for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter) {
			if (iter->second <= tk)
				list.push_back(*iter);
			else {
				if (list.size() > 0) {
					filter(list, result);
					tk -= dt;
					auto iterator{ list.cbegin() };
					for (; iterator != list.cend(); ++iterator) {
						if (iterator->second >= tk) break;
					}
					list.erase(list.cbegin(), iterator);
					tk += dt;
				}
				tk += dt;
			}
		}
		if (list.size() > 0) filter(list, result);
		return result;
	}

	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		std::list<std::pair<general::math::PV, general::time::JD>>&& measurements,
		const double dt,
		const size_t polydegree)
	{
		return filter_measurements(measurements, dt, polydegree);
	}
}