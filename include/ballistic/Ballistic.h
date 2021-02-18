#pragma once
#include "TranslationModel.h"
#include "Conversions.h"
#include "Structures.h"
#include "general/GeneralConstants.h"
#include "Integrators.h"
#include "AdamsIntegrator.h"
#include "RungeKuttaIntegrator.h"
#include <memory>
#include <algorithm>

namespace ball
{
	using namespace space;
		
	// Object provides the functionality for ballistic calculation.
	// The trajectory of the center of mass for the space vehicle can be calculated.
	template<class ModelType>
	class Ballistic
	{
	private:
		std::shared_ptr<TranslationModel<ModelType>> _pModel;
		std::vector<std::pair<general::math::PV, general::time::JD>> _trajectory;
		std::vector<size_t> _loops;
		double _deltatime;

	private:

		template<class SinglestepIntType>
		void start_run(
			const double startStep, 
			const double continueStep, 
			const size_t index, 
			const SinglestepIntegrator<SinglestepIntType, general::math::PV>& integrator)
		{
			auto x0{ _trajectory[0] }; 
			std::pair<general::math::PV, general::time::JD> xk;
			size_t loop{ _loops[0] };
			bool intersection;
			const size_t n{ static_cast<size_t>(continueStep / startStep) };
			const double step{ continueStep / n };
			const double dt{ step / general::time::SEC_PER_DAY };
			for (size_t i = 0; i < index; ++i)
			{
				for (size_t k = 0; k < n; ++k)
				{
					integrator.integrate(x0.first, x0.second, startStep, xk.first, xk.second);
					intersection = std::signbit(x0.first.Pos.Z) == true && std::signbit(xk.first.Pos.Z) == false;
					x0 = xk;
					loop += intersection ? 1 : 0;
				}
				_trajectory[i + 1] = xk;
				_loops[i + 1] = loop;
			}
		}

		template<class MultistepIntType>
		void continue_run(
			const double step,
			const size_t index,
			const MultistepIntegrator<MultistepIntType, general::math::PV>& integrator)
		{
			bool intersection;
			auto iter{ _trajectory.begin() };
			for (size_t i = index; i < _trajectory.size() - 1; ++i)
			{
				integrator.integrate(iter._Ptr, step, _trajectory[i + 1].first, _trajectory[i + 1].second);
				intersection = std::signbit(_trajectory[i].first.Pos.Z) == true && std::signbit(_trajectory[i + 1].first.Pos.Z) == false;
				_loops[i + 1] = _loops[i] + intersection ? 1 : 0;
				iter++;
			}
		}

	public:
		explicit Ballistic(const std::shared_ptr<TranslationModel<ModelType>>& pModel) :
			_pModel{ pModel },
			_deltatime{ 0 }
		{}
		Ballistic(const Ballistic& ball) = delete;
		Ballistic(Ballistic&& b) noexcept : 
			_pModel{ std::move(b._pModel) },
			_trajectory{ std::move(b._trajectory)},
			_loops{ std::move(b._loops) },
			_deltatime{ b._deltatime }
		{}
		~Ballistic() = default;

		Ballistic& operator = (const Ballistic& ball) = delete;
		Ballistic& operator = (Ballistic&& b) noexcept
		{
			_pModel = std::move(b._pModel);
			_trajectory = std::move(b._trajectory);
			_loops = std::move(b._loops);
			_deltatime = b._deltatime;
			_deltatime = 0;
			return *this;
		}

		template<class MultistepIntType, class SinglestepIntType>
		// Calculating the trajectory.
		// x0 - initial conditions;
		// tk - final time;
		// startStep - period for the single step integrator in seconds;
		// continueStep - period for the multi step integrator in seconds.
		// Throws runtime error when height is out of bounds (MinHeight, MaxHeight).
		void Run(
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

			auto func = &TranslationModel<ModelType>::function;
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

		template<class MultistepIntType, class SinglestepIntType>
		void Run(
			const State& x0, 
			const general::time::JD& tk,
			SinglestepIntegrator<SinglestepIntType, general::math::PV>&& singlestep_int,
			MultistepIntegrator<MultistepIntType, general::math::PV>&& multistep_int,
			const double startstep = 5.0,
			const double continuestep = 30.0)
		{
			Run<MultistepIntType, SinglestepIntType>(x0, tk, singlestep_int, multistep_int, startstep, continuestep);
			//auto second{ 1.0 / general::time::SEC_PER_DAY };
			//if (startstep < second || continuestep < second)
			//	throw std::invalid_argument("Invalid time step (should be > 1.0 sec)!");
			//if (startstep > continuestep)
			//	throw std::invalid_argument("Invalid start step > continue step!");
			//if (tk <= x0.T)
			//	throw std::invalid_argument("Invalid tk < tn!");
			//_deltatime = continuestep / general::time::SEC_PER_DAY;
			//_pForecast->sBall = x0.Sb;
			//size_t index = multistep_int.degree() - 1;
			//size_t count = static_cast<size_t>((tk - x0.T) / _deltatime) + 1;
			//_trajectory.resize(count);
			//_loops.resize(count);

			//auto func = &TranslationModel<ModelType>::function;
			//auto ptr = _pForecast.get();
			//auto function = [func, ptr](const general::math::PV& vec, const general::time::JD& time) {
			//	return (ptr->*func)(vec, time);
			//};
			//multistep_int.func = function;
			//singlestep_int.func = function;

			//_trajectory[0] = { x0.Vec, x0.T };
			//_loops[0] = x0.Loop;
			//// performing initial calculations
			//start_run(startstep, continuestep, index, singlestep_int);
			//// performing remaining calculations
			//continue_run(continuestep, index, multistep_int);
		}


		template<class MultistepIntType = AdamsIntegrator<>, class SinglestepIntType = RKIntegrator<>>
		// Calculating the trajectory.
		// x0 - initial conditions;
		// tk - final time;
		// startStep - period for the single step integrator in seconds;
		// continueStep - period for the multi step integrator in seconds.
		// Throws runtime error when height is out of bounds (MinHeight, MaxHeight).
		void Run(
			const State& x0,
			const general::time::JD& tk,
			const double startStep = 5.0,
			const double continueStep = 30.0)
		{
			auto second{ 1.0 / general::time::SEC_PER_DAY };
			if (startStep < second || continueStep < second)
				throw std::invalid_argument("Invalid time step (should be > 1.0 sec)!");
			if (startStep > continueStep)
				throw std::invalid_argument("Invalid start step > continue step!");
			if (tk <= x0.T)
				throw std::invalid_argument("Invalid tk < tn!");
			_deltatime = continueStep / general::time::SEC_PER_DAY;
			_pModel->sBall = x0.Sb;
			SinglestepIntType singlestep_int;
			MultistepIntType multistep_int;
			size_t index = multistep_int.degree() - 1;
			size_t count = static_cast<size_t>((tk - x0.T) / _deltatime) + 1;
			_trajectory.resize(count);
			_loops.resize(count);

			auto func = &TranslationModel<ModelType>::function;
			auto ptr = _pModel.get();
			auto function = [func, ptr](const general::math::PV& vec, const general::time::JD& time) {
				return (ptr->*func)(vec, time);
			};
			multistep_int.func = function;
			singlestep_int.func = function;

			_trajectory[0] = { x0.Vec, x0.T };
			_loops[0] = x0.Loop;
			// performing initial calculations
			start_run(startStep, continueStep, index, singlestep_int);
			// performing remaining calculations
			continue_run(continueStep, index, multistep_int);
		}

		// Calculating the state parameters.
		// time - time of the trajectory point (should be between the start and final time of the integration);
		// x - will be filled by the content.
		bool get_point(const general::time::JD& time, State& x) const
		{
			using long_t = long long;
			auto count{ _trajectory.size() };
			if (count < 4) return false;
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
				for (size_t n = 0; n < 4; ++n)
				{
					mult = 1.0;
					for (size_t m = 0; m < 4; ++m)
					{
						if (m != n)
						{
							mult *= (time - _trajectory[index + m].second) /
								(_trajectory[indexn].second - _trajectory[index + m].second);
						}
					}
					result += mult * _trajectory[indexn++].first;
				}
				intersection = intersection && std::signbit(result.Pos.Z) == false;
				x = State(result, _pModel->sBall, time, loop + intersection ? 1 : 0);
				return true;
			}
			return false;
		}
		// Returns the reference to current trajectory.
		const std::vector<std::pair<general::math::PV, general::time::JD>>& trajectory() const
		{
			return _trajectory;
		}


	};
}
