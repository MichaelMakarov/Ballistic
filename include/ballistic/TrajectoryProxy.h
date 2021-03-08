#pragma once
#include "Ballistic.h"
#include "general/Matrix.h"
#include "general/Polynom.h"
#include <vector>
#include <list>
#include <utility>
#include <string>
#include <future>
#include <iostream>


namespace ball
{
	enum class SolveType
	{
		POSITION_ONLY = 3,
		ALL_COORDINATES = 6
	};

	std::vector<general::math::VectorDyn> create_array_of_vectors(const size_t count, const size_t size);

	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		std::list<std::pair<general::math::PV, general::time::JD>>& measurements,
		const double dt = 300 / 86400.0);

	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		std::list<std::pair<general::math::PV, general::time::JD>>&& measurements,
		const double dt = 300 / 86400.0);

	template<class Iterator>
	std::list<std::pair<general::math::PV, general::time::JD>> filter_measurements(
		Iterator first, Iterator last,
		const double dt = 300 / 86400.0)
	{
		return filter_measurements(std::list<std::pair<general::math::PV, general::time::JD>>(first, last), dt);
	}

	template<class T> concept StdContainer = requires (T ex) { 
		{ ex.size() }; 
		{ T::iterator }; 
		{ ex.cbegin() }; 
		{ ex.cend() }; 
	};

	template<TransModel Model, StdContainer Container, class ... Args>
	size_t refine_initpoint(
		State& x0,
		const Container& measurements,
		std::function<std::unique_ptr<Model>()> create_model,
		const SolveType type = SolveType::POSITION_ONLY,
		const size_t iterations = 10
	)
	{
		using namespace general::math;
		using namespace general::time;
#define DEBUG_PROXY_STREAM std::cout

#ifdef DEBUG_PROXY_STREAM
		DEBUG_PROXY_STREAM << x0 << std::endl;
		DEBUG_PROXY_STREAM << "measurements:\n";
#endif
		const size_t vars = static_cast<size_t>(type);
		size_t iteration{ 0 };
		const double deltas[7]{ 25, 25, 25, 0.025, 0.025, 0.025, 0.2 * 0.0007855809 };
		general::time::JD tk;
		for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter) {
			if (iter->second > tk) tk = iter->second;
#ifdef DEBUG_PROXY_STREAM
			DEBUG_PROXY_STREAM << iter->second.to_datetime() << "; " << iter->first << std::endl;
#endif
		}
		auto B{ MatrixDyn(vars * measurements.size(), 7) };
		auto L{ VectorDyn(vars * measurements.size()) };

		auto calc_trajectory = [create_model, measurements](const State& x0, const JD& tk) -> std::vector<State> {
			auto ball = create_forecast<Model>(create_model());
			ball.run<AdamsIntegrator<PV>, RKIntegrator<PV>>(
				x0, tk,
				RKIntegrator<PV>(),
				AdamsIntegrator<PV>());
			auto traj{ std::vector<State>(measurements.size()) };
			size_t index{ 0 };
			for (auto& m : measurements)
				ball.get_point(m.second, traj[index++]);
			return traj;
		};
		auto calc_diversities = [measurements, vars](const std::vector<State>& trajectory, VectorDyn& vec) {
			size_t index{ 0 }, count{ 0 };
			for (auto iter = measurements.cbegin(); iter != measurements.cend(); ++iter, ++count)
				for (size_t n = 0; n < vars; ++n)
					vec[index++] = iter->first[n] - trajectory[count].Vec[n];
		};
		auto calc_partial_derivatives = [deltas, calc_trajectory, vars](const State& x0, const JD& tk, MatrixDyn& matrix)
			-> std::vector<State> {
			State xlist[7]{ x0, x0, x0, x0, x0, x0, x0 };
			auto traj = calc_trajectory(x0, tk);
			std::future<std::vector<State>> futures[6];
			for (size_t i = 0; i < 6; ++i) {
				xlist[i].Vec[i] += deltas[i];
				futures[i] = std::async(std::launch::async, calc_trajectory, xlist[i], tk);
			}
			xlist[6].Sb += deltas[6];
			auto vared = calc_trajectory(xlist[6], tk);
			size_t index{ 0 };
			for (size_t n = 0; n < traj.size(); ++n)
				for (size_t k = 0; k < vars; ++k)
					matrix(index++, 6) = (vared[n].Vec[k] - traj[n].Vec[k]) / deltas[6];
			for (size_t i = 0; i < 6; ++i) {
				index = 0;
				vared = futures[i].get();
				for (size_t n = 0; n < traj.size(); ++n)
					for (size_t k = 0; k < vars; ++k)
						matrix(index++, i) = (vared[n].Vec[k] - traj[n].Vec[k]) / deltas[i];
			}
			return traj;
		};

		while (iteration++ < iterations) {
			try {
				auto trajectory = calc_partial_derivatives(x0, tk, B);
				calc_diversities(trajectory, L);
				auto BT = transpose(B);
				auto M = BT * B;
				auto D{ MatrixDyn(M.rows(), M.columns()) };
				for (size_t i = 0; i < 7; ++i) D(i, i) = 1 / std::sqrt(M(i, i));
				M = AxD(DxA(D, M), D);
				auto dx0 = AxD(DxA(D, inverse(M)), D) * BT * L;
				x0.Vec.Pos.X += dx0[0];
				x0.Vec.Pos.Y += dx0[1];
				x0.Vec.Pos.Z += dx0[2];
				x0.Vec.Vel.X += dx0[3];
				x0.Vec.Vel.Y += dx0[4];
				x0.Vec.Vel.Z += dx0[5];
				x0.Sb += dx0[6];
				auto V = B * dx0 - L;
				double mV = (V * V) / V.length();
				double mL = (L * L) / L.length();
#ifdef DEBUG_PROXY_STREAM
				DEBUG_PROXY_STREAM << "iteration: " << iteration << std::endl;
				DEBUG_PROXY_STREAM << "diversities: " << L << std::endl;
				DEBUG_PROXY_STREAM << "dx0: " << dx0 << std::endl;
#endif
				if (std::abs((mV - mL) / mV) < 1e-2) break;
			}
			catch (const std::exception) { return 0; }
		}
		return iteration;
	}
}