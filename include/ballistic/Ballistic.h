#pragma once
#include "TranslationModel.h"
#include "Structures.h"
#include "general/GeneralConstants.h"
#include "Integrators.h"
#include <memory>
#include <algorithm>

namespace ball
{
	//using namespace space;

	template<class T>
	concept TransModel = std::is_base_of<TranslationModel<T>, T>::value;
	template<class Int>
	concept PVSinglestepInt = std::is_base_of<SinglestepIntegrator<Int, general::math::PV, general::time::JD>, Int>::value;
	template<class Int>
	concept PVMultistepInt = std::is_base_of<MultistepIntegrator<Int, general::math::PV, general::time::JD>, Int>::value;
	
	/// <summary>
	/// Object provides the functionality for ballistic calculation.
	/// The trajectory of the center of mass for the space vehicle can be calculated.
	/// </summary>
	template<TransModel Model>
	class Forecast
	{
	private:
		std::unique_ptr<TranslationModel<Model>> _pModel;
		std::vector<general::math::PV> _trajectory;
		std::vector<general::time::JD> _times;
		std::vector<size_t> _loops;
		double _deltatime;

	private:
		template<PVSinglestepInt SinglestepInt>
		void start_run(
			const double start_step, 
			const double continue_step, 
			const size_t index, 
			const SinglestepIntegrator<SinglestepInt, general::math::PV, general::time::JD>& integrator)
		{
			auto x0{ _trajectory[0] };
			auto t0{ _times[0] };
			auto xk{ general::math::PV() };
			auto tk{ general::time::JD() };
			auto func{ make_memberfunc(&TranslationModel<Model>::function, _pModel.get()) };
			bool intersection;
			const size_t n{ static_cast<size_t>(continue_step / start_step) };
			const double step{ continue_step / n };
			for (size_t i = 0; i < index; ++i) {
				_loops[i + 1] = _loops[i];
				for (size_t k = 0; k < n; ++k) {
					integrator.integrate(x0, t0, step, xk, tk, func);
					intersection = std::signbit(x0.pos.z()) == true && std::signbit(xk.pos.z()) == false;
					x0 = xk;
					t0 = tk;
					if (intersection) _loops[i + 1] += 1;
				}
				_trajectory[i + 1] = xk;
				_times[i + 1] = tk;
			}
		}

		template<PVMultistepInt MultistepInt>
		void continue_run(
			const double step,
			const size_t index,
			const MultistepIntegrator<MultistepInt, general::math::PV, general::time::JD>& integrator)
		{
			bool intersection;
			auto x_iter{ _trajectory.cbegin() };
			auto t_iter{ _times.cbegin() };
			auto func{ make_memberfunc(&TranslationModel<Model>::function, _pModel.get()) };
			for (size_t i = index; i < _trajectory.size() - 1; ++i) {
				integrator.integrate(x_iter._Ptr, t_iter._Ptr, step, _trajectory[i + 1], _times[i + 1], func);
				intersection = std::signbit(_trajectory[i].pos.z()) == true && std::signbit(_trajectory[i + 1].pos.z()) == false;
				_loops[i + 1] = _loops[i] + (intersection ? 1 : 0);
				++x_iter;
				++t_iter;
			}
		}

		template<PVSinglestepInt SinglestepInt>
		void single_run(
			const double step,
			const SinglestepIntegrator<SinglestepInt, general::math::PV, general::time::JD>& integrator)
		{
			bool intersection;
			auto func{ make_memberfunc(&TranslationModel<Model>::function, _pModel.get()) };
			for (size_t i = 0; i < _trajectory.size() - 1; ++i) {
				const auto& x0{ _trajectory[i] };
				const auto& t0{ _times[i] };
				auto& xk{ _trajectory[i + 1] };
				auto& tk{ _times[i + 1] };
				integrator.integrate(x0, t0, step, xk, tk, func);
				intersection = std::signbit(x0.pos.z()) == true && std::signbit(xk.pos.z()) == false;
				_loops[i + 1] = _loops[i] + (intersection ? 1 : 0);
			}
		}

	public:
		/// <summary>
		/// Initializing a ballistic object 
		/// </summary>
		/// <param name="pModel">A pointer to the model of the linear movement</param>
		explicit Forecast(std::unique_ptr<TranslationModel<Model>>&& pModel) :
			_pModel{ std::move(pModel) },
			_deltatime{ 0 }
		{}
		Forecast(const Forecast& ball) = delete;
		Forecast(Forecast&& b) noexcept : 
			_pModel{ std::move(b._pModel) },
			_trajectory{ std::move(b._trajectory)},
			_loops{ std::move(b._loops) },
			_deltatime{ b._deltatime }
		{}
		~Forecast() = default;

		Forecast& operator = (const Forecast& ball) = delete;
		Forecast& operator = (Forecast&& b) noexcept
		{
			_pModel = std::move(b._pModel);
			_trajectory = std::move(b._trajectory);
			_loops = std::move(b._loops);
			_deltatime = b._deltatime;
			_deltatime = 0;
			return *this;
		}

		/// <summary>
		/// Run the trajectory calculation. 
		/// Throws an exception when height is out of bounds (MinHeight, MaxHeight) or input arguments are invalid.
		/// </summary>
		/// <typeparam name="MultistepIntType">type of the multistep integrator that inherits the MultistepIntegrator base class</typeparam>
		/// <typeparam name="SinglestepIntType">type of the multistep integrator that inherits the SinglestepIntegrator base class</typeparam>
		/// <param name="x0">an initial point</param>
		/// <param name="tk">final time of the calculation</param>
		/// <param name="singlestep_int">a multistep intergrator that will be used</param>
		/// <param name="multistep_int">a single step integrator that will be used to calculate initial points</param>
		/// <param name="step">time step for the multistep integrator</param>
		template<PVMultistepInt MultistepInt, PVSinglestepInt SinglestepInt>
		void run(
			const State& x0,
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepInt, general::math::PV, general::time::JD>& singlestep_int,
			MultistepIntegrator<MultistepInt, general::math::PV, general::time::JD>& multistep_int,
			const double step = 30.0)
		{
			if (step < 1.0)
				throw std::invalid_argument("Invalid time step (should be > 1.0 sec)!");
			if (tk <= x0.T)
				throw std::invalid_argument("Invalid tk < tn!");
			_deltatime = step / general::time::SEC_PER_DAY;
			size_t index = multistep_int.degree() - 1;
			size_t count = static_cast<size_t>((tk - x0.T) / _deltatime) + 1;
			_trajectory.resize(count);
			_times.resize(count);
			_loops.resize(count);

			_trajectory[0] = x0.Vec;
			_times[0] = x0.T;
			_loops[0] = x0.Loop;
			_pModel->sBall = x0.Sb;
			// performing initial calculations
			start_run(step / 6, step, index, singlestep_int);
			// performing remaining calculations
			continue_run(step, index, multistep_int);
		}

		/// <summary>
		/// Run the trajectory calculation. 
		/// Throws an exception when height is out of bounds (MinHeight, MaxHeight) or input arguments are invalid.
		/// </summary>
		/// <typeparam name="MultistepInt">type of the multistep integrator that inherits the MultistepIntegrator base class</typeparam>
		/// <typeparam name="SinglestepInt">type of the multistep integrator that inherits the SinglestepIntegrator base class</typeparam>
		/// <param name="x0">an initial point</param>
		/// <param name="tk">final time of the calculation</param>
		/// <param name="singlestep_int">a multistep intergrator that will be used</param>
		/// <param name="multistep_int">a single step integrator that will be used to calculate initial points</param>
		/// <param name="step">time step for the multistep integrator</param>
		template<PVMultistepInt MultistepInt, PVSinglestepInt SinglestepInt>
		void run(
			const State& x0, 
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepInt, general::math::PV, general::time::JD>&& singlestep_int,
			MultistepIntegrator<MultistepInt, general::math::PV, general::time::JD>&& multistep_int,
			const double step = 30.0)
		{
			run<MultistepInt, SinglestepInt>(x0, tk, singlestep_int, multistep_int, step);
		}

		template<PVSinglestepInt SinglestepInt>
		void run(
			const State& x0,
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepInt, general::math::PV, general::time::JD>& integrator,
			const double step = 30.0)
		{
			if (step < 1.0)
				throw std::invalid_argument("Invalid time step (should be > 1.0 sec)!");
			if (tk <= x0.T)
				throw std::invalid_argument("Invalid tk < tn!");
			_deltatime = step / general::time::SEC_PER_DAY;
			size_t count = static_cast<size_t>((tk - x0.T) / _deltatime) + 1;
			_trajectory.resize(count);
			_times.resize(count);
			_loops.resize(count);

			_trajectory[0] = x0.Vec;
			_times[0] = x0.T;
			_loops[0] = x0.Loop;
			_pModel->sBall = x0.Sb;

			single_run(step, integrator);
		}

		template<PVSinglestepInt SinglestepInt>
		void run(
			const State& x0,
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepInt, general::math::PV, general::time::JD>&& integrator,
			const double step = 30.0)
		{
			run<SinglestepInt>(x0, tk, integrator, step);
		}

		// Calculating the state parameters.
		// time - time of the trajectory point (should be between the start and final time of the integration);
		// x - will be filled by the content.
		State get_point(const general::time::JD& time) const
		{
			using long_t = long long;
			auto count{ _trajectory.size() };
			if (count < 4) throw std::length_error("Not enough points of trajectory!");
			auto index = static_cast<size_t>((time - _times[0]) / _deltatime);
			if (index < count) {
				// checking the value sign of the z coordinate of the nearest point
				bool intersection = std::signbit(_trajectory[index].pos.z()) == true;
				const size_t loop{ _loops[index] };
				// get the index of the first point for the interpolation
				index = static_cast<size_t>(
					std::max(
						long_t(0),
						std::min(static_cast<long_t>(count) - 4, static_cast<long_t>(index) - 2)
					));
				// hence the interpolation performing
				// P(t) = sum{n = 0..4} (mult{m = 0..4, m != n} (t - t_m)/(t_n - t_m)) x_n
				general::math::PV result;
				double mult;
				size_t indexn{ index };
				for (size_t n = 0; n < 4; ++n) {
					mult = 1.0;
					for (size_t m = 0; m < 4; ++m) {
						if (m != n)	{
							mult *= (time - _times[index + m]) / (_times[indexn] - _times[index + m]);
						}
					}
					result += mult * _trajectory[indexn++];
				}
				intersection = intersection && std::signbit(result.pos.z()) == false;
				return State(result, _pModel->sBall, time, loop + intersection ? 1 : 0);
			}
			throw std::invalid_argument("Time is out of bounds!");
		}

		/// <summary>
		/// The reference to current trajectory
		/// </summary>
		/// <returns>an array of points by reference</returns>
		const std::vector<general::math::PV>& trajectory() const { return _trajectory; }
		/// <summary>
		/// The reference to current list of times corresponding to trajectory points
		/// </summary>
		/// <returns>a list of time moments by reference</returns>
		const std::vector<general::time::JD>& times() const { return _times; }
	};
	/// <summary>
	/// Fast creating ball class example
	/// </summary>
	/// <param name="pModel">a model of the movement</param>
	template<TransModel Model>
	auto make_forecast(std::unique_ptr<Model>&& pModel) -> Forecast<Model>
	{
		return Forecast<Model>(std::forward<std::unique_ptr<Model>>(pModel));
	}
}
