#pragma once
#include "Ballistic.h"
#include "general/Matrix.h"
#include <list>
#include <utility>
#include <string>
#include <future>


namespace ball
{
	enum class SolveType
	{
		POSITION_ONLY = 3,
		ALL_COORDINATES = 6
	};

	std::vector<general::math::Vector> create_array_of_vectors(const size_t count, const size_t size);

	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		std::list<std::pair<general::math::PV, general::time::JD>>& measurements,
		const double dt = 300 / 86400.0,
		const size_t polydegree = 4);

	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		std::list<std::pair<general::math::PV, general::time::JD>>&& measurements,
		const double dt = 300 / 86400.0,
		const size_t polydegree = 4);

	template<class Iterator>
	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		Iterator first, Iterator last,
		const double dt = 300 / 86400.0,
		const size_t polydegree = 4)
	{
		return filter_measurements(std::list<std::pair<general::math::PV, general::time::JD>>(first, last), dt, polydegree);
	}

	/// <summary>
	/// A concept of the container from standard library
	/// </summary>
	template<class T, class M> concept StdContainer = requires (T ex) { 
		{ ex.size() }; 
		{ T::iterator }; 
		{ ex.cbegin() }; 
		{ ex.cend() }; 
	};
	/// <summary>
	/// Calculating the residuals
	/// </summary>
	/// <param name="trajectory"> - a calculated trajectory</param>
	/// <param name="measurements"> - an amount of measuring data</param>
	/// <param name="vars"> - a number of coordinates to consider</param>
	/// <param name="vec"> - a vector of residuals</param>
	template<StdContainer<std::pair<general::math::PV, general::time::JD>> Container>
	void _calc_residuals(
		const std::vector<State>& trajectory,
		const Container& measurements,
		const size_t vars,
		general::math::Vector& vec)
	{
		size_t index{ 0 }, count{ 0 };
		for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter, ++count)
			for (size_t n = 0; n < vars; ++n)
				vec[index++] = iter->first[n] - trajectory[count].Vec[n];
	}
	/// <summary>
	/// Calculating the matrix of partial derivatives
	/// </summary>
	/// <param name="x0"> - an initial point</param>
	/// <param name="measurements"> - an amount of measuring data</param>
	/// <param name="calc_trajectory"> - a function for calculating the trajectory of motion</param>
	/// <param name="vars"> - a number of coordinates to consider</param>
	/// <param name="matrix"> - a matrix of partial derivatives</param>
	/// <returns>a calculated trajectory</returns>
	template<StdContainer<std::pair<general::math::PV, general::time::JD>> Container>
	std::vector<State> _calc_partial_derivatives(
		const State& x0,
		const Container& measurements,
		const std::function<std::vector<State>(const State&, const Container&)>& calc_trajectory,
		const size_t vars,
		general::math::Matrix& matrix)
	{
		const double deltas[7]{ 25, 25, 25, 0.025, 0.025, 0.025, 0.2 * 0.0007855809 };
		State xlist[7]{ x0, x0, x0, x0, x0, x0, x0 };
		auto traj = calc_trajectory(x0, measurements);
		std::future<std::vector<State>> futures[6];
		for (size_t i = 0; i < 6; ++i) {
			xlist[i].Vec[i] += deltas[i];
			futures[i] = std::async(std::launch::async, calc_trajectory, xlist[i], measurements);
		}
		xlist[6].Sb += deltas[6];
		auto varied = calc_trajectory(xlist[6], measurements);
		size_t index{ 0 };
		for (size_t n = 0; n < traj.size(); ++n)
			for (size_t k = 0; k < vars; ++k)
				matrix(index++, 6) = (varied[n].Vec[k] - traj[n].Vec[k]) / deltas[6];
		for (size_t i = 0; i < 6; ++i) {
			index = 0;
			varied = futures[i].get();
			for (size_t n = 0; n < traj.size(); ++n)
				for (size_t k = 0; k < vars; ++k)
					matrix(index++, i) = (varied[n].Vec[k] - traj[n].Vec[k]) / deltas[i];
		}
		return traj;
	}
	/// <summary>
	/// Recalculating an initial point to adjust the trajectory of motion to measurements
	/// </summary>
	/// <param name="x0"> - initial point</param>
	/// <param name="measurements"> - an amount of measurements</param>
	/// <param name="calc_trajectory"> - a function for trajectory calculating</param>
	/// <param name="type"> - a solution type</param>
	/// <param name="iterations"> - a max number of iterations for optimization</param>
	/// <param name="pLogstrbuf"> - a pointer to the logging stream</param>
	/// <returns>a number of iterations required for solution</returns>
	template<StdContainer<std::pair<general::math::PV, general::time::JD>> Container>
	size_t refine_initpoint(
		State& x0,
		const Container& measurements,
		const std::function<std::vector<State>(const State&, const Container&)>& calc_trajectory,
		const SolveType type = SolveType::POSITION_ONLY,
		const size_t iterations = 10,
		std::streambuf* pLogstrbuf = nullptr
	)
	{
		using namespace general::math;
		using namespace general::time;

		auto logout = std::ostream(pLogstrbuf);
		const size_t vars = static_cast<size_t>(type);
		size_t iteration{ 0 };
		logout << "Initial: " << x0 << std::endl;
		logout << "Measurements:\n";
		for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter) {
			logout << iter->first << " " << iter->second << std::endl;
		}
		auto B{ Matrix(vars * measurements.size(), 7) };
		auto L{ Vector(vars * measurements.size()) };

		while (iteration++ < iterations) {
			auto trajectory = _calc_partial_derivatives(x0, measurements, calc_trajectory, vars, B);
			_calc_residuals(trajectory, measurements, vars, L);
			auto BT = transpose(B);
			auto M = BT * B;
			auto D{ Matrix(M.rows(), M.columns()) };
			for (size_t i = 0; i < 7; ++i) D(i, i) = 1 / std::sqrt(M(i, i));
			M = AxD(DxA(D, M), D);
			auto dx0 = AxD(DxA(D, inverse(M)), D) * BT * L;
			x0.Vec.pos.x() += dx0[0];
			x0.Vec.pos.y() += dx0[1];
			x0.Vec.pos.z() += dx0[2];
			x0.Vec.vel.x() += dx0[3];
			x0.Vec.vel.y() += dx0[4];
			x0.Vec.vel.z() += dx0[5];
			x0.Sb += dx0[6];
			logout << "Iteration " << iteration << "\nResiduals: " << L << std::endl;
			auto V = B * dx0 - L;
			double mV = V.length();
			double mL = L.length();
			if (std::abs(mV - mL) / mL < 1e-3) break;
		}
		return iteration;
	}
}