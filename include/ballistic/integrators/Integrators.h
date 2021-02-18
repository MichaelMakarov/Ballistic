#pragma once
#include "general/Geometry.h"
#include "general/Times.h"
#include <functional>
#include <type_traits>
#include <concepts>

namespace ball
{
	template<class T> concept Arithmetic = 
		std::is_default_constructible<T>::value &&
		requires (T a, double b) {
			{ a * b };
			{ a / b };
			{ a *= b };
			{ a /= b };
			{ b * a };
		} && requires (T a, T b) { 
			{ a += b };
			{ a -= b };
			{ a + b };
			{ a - b };
		};

	template<class R>
	class Integrator
	{
	public:
		Integrator() = default;
		~Integrator() = default;

		std::function<R(const R&, const general::time::JD& t)> func;
	};

	template<class IntType, Arithmetic R>
	class SinglestepIntegrator : public Integrator<R>
	{
	public:
		SinglestepIntegrator() : Integrator<R>() {}
		~SinglestepIntegrator() = default;

		void integrate(
			const R& x0,
			const general::time::JD& t0,
			const double step,
			R& xk,
			general::time::JD& tk) const
		{
			static_cast<const IntType*>(this)->integrate(x0, t0, step, xk, tk);
		}
	};

	template<class IntType, Arithmetic R>
	class MultistepIntegrator : public Integrator<R>
	{
	private:
		const size_t _degree;

	public:
		explicit MultistepIntegrator(const size_t degree) : Integrator<R>(), _degree(degree) {}
		~MultistepIntegrator() = default;

		void integrate(
			std::pair<R, general::time::JD>* pData,
			const double step,
			R& xk,
			general::time::JD& tk) const
		{
			static_cast<const IntType*>(this)->integrate(pData, step, xk, tk);
		}

		size_t degree() const {	return _degree;	}
	};
}