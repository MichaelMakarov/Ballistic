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
	template<class T>
	concept TransSinglestepInt = std::is_base_of<SinglestepIntegrator<T, general::math::PV>, T>::value;
	template<class T>
	concept TransMultistepInt = std::is_base_of<MultistepIntegrator<T, general::math::PV>, T>::value;
		
	// Object provides the functionality for ballistic calculation.
	// The trajectory of the center of mass for the space vehicle can be calculated.
	template<TransModel Model>
	class Forecast
	{
	private:
		std::unique_ptr<TranslationModel<Model>> _pModel;
		std::vector<std::pair<general::math::PV, general::time::JD>> _trajectory;
		std::vector<size_t> _loops;
		double _deltatime;

	private:
		template<TransSinglestepInt SinglestepIntType>
		void start_run(
			const double start_step, 
			const double continue_step, 
			const size_t index, 
			const SinglestepIntegrator<SinglestepIntType, general::math::PV>& integrator)
		{
			auto x0{ _trajectory[0] }; 
			std::pair<general::math::PV, general::time::JD> xk;
			//size_t loop{ _loops[0] };
			bool intersection;
			const size_t n{ static_cast<size_t>(continue_step / start_step) };
			const double step{ continue_step / n };
			const double dt{ step / general::time::SEC_PER_DAY };
			for (size_t i = 0; i < index; ++i) {
				for (size_t k = 0; k < n; ++k) {
					integrator.integrate(x0.first, x0.second, start_step, xk.first, xk.second);
					intersection = std::signbit(x0.first.Pos.Z) == true && std::signbit(xk.first.Pos.Z) == false;
					x0 = xk;
					_loops[i + 1] = _loops[i] + (intersection ? 1 : 0);
				}
				_trajectory[i + 1] = xk;
			}
		}

		template<TransMultistepInt MultistepIntType>
		void continue_run(
			const double step,
			const size_t index,
			const MultistepIntegrator<MultistepIntType, general::math::PV>& integrator)
		{
			bool intersection;
			auto iter{ _trajectory.begin() };
			for (size_t i = index; i < _trajectory.size() - 1; ++i) {
				integrator.integrate(iter._Ptr, step, _trajectory[i + 1].first, _trajectory[i + 1].second);
				intersection = std::signbit(_trajectory[i].first.Pos.Z) == true && std::signbit(_trajectory[i + 1].first.Pos.Z) == false;
				_loops[i + 1] = _loops[i] + (intersection ? 1 : 0);
				iter++;
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
		/// <param name="startstep">time step for the single step integrator</param>
		/// <param name="continuestep">time step for the multistep integrator</param>
		template<TransMultistepInt MultistepIntType, TransSinglestepInt SinglestepIntType>
		void run(
			const State& x0,
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepIntType, general::math::PV>& singlestep_int,
			MultistepIntegrator<MultistepIntType, general::math::PV>& multistep_int,
			const double startstep = 5.0,
			const double continuestep = 30.0)
		{
			auto second{ 1.0 / general::time::SEC_PER_DAY };
			if (startstep < second || continuestep < second)
				throw std::invalid_argument("Invalid time step (should be > 1.0 sec)!");
			if (startstep > continuestep)
				throw std::invalid_argument("Invalid start step > continue step!");
			if (tk <= x0.T)
				throw std::invalid_argument("Invalid tk < tn!");
			_deltatime = continuestep / general::time::SEC_PER_DAY;
			_pModel->sBall = x0.Sb;
			size_t index = multistep_int.degree() - 1;
			size_t count = static_cast<size_t>((tk - x0.T) / _deltatime) + 1;
			_trajectory.resize(count);
			_loops.resize(count);

			auto func = &TranslationModel<Model>::function;
			auto ptr = _pModel.get();
			auto function = [func, ptr](const general::math::PV& vec, const general::time::JD& time) {
				return (ptr->*func)(vec, time);
			};
			multistep_int.func = function;
			singlestep_int.func = function;

			_trajectory[0] = { x0.Vec, x0.T };
			_loops[0] = x0.Loop;
			// performing initial calculations
			start_run(startstep, continuestep, index, singlestep_int);
			// performing remaining calculations
			continue_run(continuestep, index, multistep_int);
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
		/// <param name="startstep">time step for the single step integrator</param>
		/// <param name="continuestep">time step for the multistep integrator</param>
		template<TransMultistepInt MultistepIntType, TransSinglestepInt SinglestepIntType>
		void run(
			const State& x0, 
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepIntType, general::math::PV>&& singlestep_int,
			MultistepIntegrator<MultistepIntType, general::math::PV>&& multistep_int,
			const double startstep = 5.0,
			const double continuestep = 30.0)
		{
			run<MultistepIntType, SinglestepIntType>(x0, tk, singlestep_int, multistep_int, startstep, continuestep);
		}

		// Calculating the state parameters.
		// time - time of the trajectory point (should be between the start and final time of the integration);
		// x - will be filled by the content.
		State get_point(const general::time::JD& time) const
		{
			using long_t = long long;
			auto count{ _trajectory.size() };
			if (count < 4) throw std::length_error("Not enough points of trajectory!");
			auto index = static_cast<size_t>((time - _trajectory[0].second) / _deltatime);
			if (index < count)
			{
				// checking the value sign of the z coordinate of the nearest point
				bool intersection = std::signbit(_trajectory[index].first.Pos.Z) == true;
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
							mult *= (time - _trajectory[index + m].second) /
								(_trajectory[indexn].second - _trajectory[index + m].second);
						}
					}
					result += mult * _trajectory[indexn++].first;
				}
				intersection = intersection && std::signbit(result.Pos.Z) == false;
				return State(result, _pModel->sBall, time, loop + intersection ? 1 : 0);
			}
			throw std::invalid_argument("Time out of bounds of integrated trajectory!");
		}
		// Returns the reference to current trajectory.
		const std::vector<std::pair<general::math::PV, general::time::JD>>& trajectory() const
		{
			return _trajectory;
		}
	};
	/// <summary>
	/// Fast creating ball class example
	/// </summary>
	/// <param name="pModel">a model of the movement</param>
	template<TransModel Model>
	constexpr auto make_forecast(std::unique_ptr<Model>&& pModel) -> Forecast<Model>
	{
		return Forecast<Model>(std::forward<std::unique_ptr<Model>>(pModel));
	}
}
