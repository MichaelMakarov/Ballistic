#pragma once
#include "TranslationModel.h"
#include "Structures.h"
#include "general/GeneralConstants.h"
#include "Integrators.h"
#include <memory>
#include <algorithm>

namespace ball
{
	/// <summary>
	/// Implements calculating, store and retrieving of orbit data functionality
	/// </summary>
	template<size_t dim> class Forecast
	{
		using R = general::math::Vec<dim>;
	private:
		std::vector<R> _trajectory;
		std::vector<general::time::JD> _times;
		std::vector<size_t> _loops;
		double _deltatime;

	public:
		Forecast() : _deltatime{ 0 } {}
		Forecast(const Forecast& other) = default;
		Forecast(Forecast&& other) noexcept : 
			_trajectory{ std::move(other._trajectory) },
			_times{ std::move(other._times) },
			_loops{ std::move(other._loops) },
			_deltatime{ std::move(other._deltatime) }
		{}
		~Forecast() = default;

		Forecast& operator= (const Forecast& other) = default;
		Forecast& operator= (Forecast&& other) noexcept
		{
			_trajectory = std::move(other._trajectory);
			_times = std::move( other._times);
			_loops = std::move(other._loops);
			_deltatime = std::move(other._deltatime);
			return *this;
		}
		/// <summary>
		/// An orbit trajectory
		/// </summary>
		/// <returns>an array of vectors</returns>
		const std::vector<R>& get_points() const { return _trajectory; }
		/// <summary>
		/// Time moments related to trajectory
		/// </summary>
		/// <returns>an arrayt of time moments</returns>
		const std::vector<general::time::JD>& get_times() const { return _times; }
		/// <summary>
		/// Loop numbers related to trajectory
		/// </summary>
		/// <returns>an array of loop numbers</returns>
		const std::vector<size_t>& get_loops() const { return _loops; }
		/// <summary>
		/// Retreiving an orbit point
		/// </summary>
		/// <param name="time">is a julian date</param>
		/// <returns>corresponding orbit data</returns>
		template<size_t degree = 4>	Params<dim> get_point(const general::time::JD& time) const
		{
			using llong = long long;
			auto count{ _trajectory.size() };
			if (count < degree) throw std::length_error("Not enough points of trajectory!");
			auto index = static_cast<size_t>((time - _times[0]) / _deltatime);
			if (index < count) {
				// checking the value sign of the z coordinate of the nearest point
				bool intersection = std::signbit(_trajectory[index][2]);
				const size_t loop{ _loops[index] };
				// get the index of the first point for the interpolation
				index = static_cast<size_t>(
					std::max(
						llong(0),
						std::min(llong(count) - llong(degree), llong(index) - llong(degree / 2))
					));
				// hence the interpolation performing
				// P(t) = sum{n = 0..4} (mult{m = 0..4, m != n} (t - t_m)/(t_n - t_m)) x_n
				R result;
				double mult;
				size_t indexn{ index };
				for (size_t n = 0; n < degree; ++n) {
					mult = 1.0;
					for (size_t m = 0; m < degree; ++m) {
						if (m != n) {
							mult *= (time - _times[index + m]) / (_times[indexn] - _times[index + m]);
						}
					}
					result += mult * _trajectory[indexn++];
				}
				intersection &= !std::signbit(result[2]);
				return Params<dim>{ result, time, 0.0, loop + static_cast<size_t>(intersection) };
			}
			throw std::invalid_argument("Time is out of bounds!");
		}
		/// <summary>
		/// Calculating an orbit
		/// </summary>
		/// <typeparam name="Intss">is a single step integrator type</typeparam>
		/// <typeparam name="Intms">is a multiple step integrator type</typeparam>
		/// <param name="x0">is initial orbit point</param>
		/// <param name="tk">is final julian date for calculation</param>
		/// <param name="func">is right part of </param>
		/// <param name="singlestep_int"></param>
		/// <param name="multistep_int"></param>
		/// <param name="deltatime"></param>
		template<class Inv, class Intss, class Intms> void run(
			const Params<dim>& x0, const general::time::JD& tk,
			const invoker<Inv, R, const R&, const general::time::JD&>& func,
			const singlestep_integrator<Intss, R, general::time::JD>& singlestep_int,
			const multistep_integrator<Intms, R, general::time::JD>& multistep_int,
			const double deltatime = 30)
		{
			if (deltatime < 0.0)
				throw std::invalid_argument("Invalid deltatime " + std::to_string(deltatime) + "!");
			if (tk < x0.T + deltatime)
				throw std::invalid_argument("Invalid tk < tn + deltatime!");
			size_t index = multistep_int.degree() - 1;
			size_t count = static_cast<size_t>((tk - x0.T) / deltatime) + 1;
			_trajectory.resize(count);
			_times.resize(count);
			_loops.resize(count);
			_trajectory[0] = x0.vec;
			_times[0] = x0.T;
			_loops[0] = x0.loop;
			_deltatime = deltatime;
			// performing initial calculations
			start_run(index, singlestep_int, func);
			// performing remaining calculations
			continue_run(index, multistep_int, func);
		}

		private:
			template<class Inv, class Intss>
			void start_run(
				const size_t index,
				const singlestep_integrator<Intss, R, general::time::JD>& integrator,
				const invoker<Inv, R, const R&, const general::time::JD&>& func)
			{
				auto x0{ _trajectory[0] };
				auto t0{ _times[0] };
				R xk;
				general::time::JD tk;
				bool intersection;
				const size_t n{ 6 };
				const double step{ _deltatime / n };
				for (size_t i = 0; i < index; ++i) {
					_loops[i + 1] = _loops[i];
					for (size_t k = 0; k < n; ++k) {
						integrator.integrate(x0, t0, step, xk, tk, func);
						intersection = std::signbit(x0[2]) && !std::signbit(xk[2]);
						x0 = xk;
						t0 = tk;
						_loops[i + 1] += static_cast<size_t>(intersection);
					}
					_trajectory[i + 1] = xk;
					_times[i + 1] = tk;
				}
			}

			template<class Inv, class Intms>
			void continue_run(
				const size_t index,
				const multistep_integrator<Intms, R, general::time::JD>& integrator,
				const invoker<Inv, R, const R&, const general::time::JD&>& func)
			{
				bool intersection;
				auto x_iter{ _trajectory.cbegin() };
				auto t_iter{ _times.cbegin() };
				for (size_t i = index; i < _trajectory.size() - 1; ++i) {
					integrator.integrate(x_iter._Ptr, t_iter._Ptr, _deltatime, _trajectory[i + 1], _times[i + 1], func);
					intersection = std::signbit(_trajectory[i][2]) && !std::signbit(_trajectory[i + 1][2]);
					_loops[i + 1] = _loops[i] + static_cast<size_t>(intersection);
					++x_iter;
					++t_iter;
				}
			}

	};

	
}
