#pragma once
#include "general/Geometry.h"
#include "general/Times.h"
#include <functional>

namespace ball
{
	template<class ... ArgsType>
	class Integrator
	{
	public:
		Integrator() {}
		~Integrator() {}

		std::function<general::math::PV(const general::math::PV&, const general::time::JD& t, const ArgsType& ...)> func;
	};

	template<class InvType, class ... ArgsType>
	class SinglestepIntegrator : public Integrator<ArgsType ...>
	{
	public:
		SinglestepIntegrator() : Integrator<ArgsType ...>() {}
		~SinglestepIntegrator() {}

		void integrate(
			const general::math::PV& x0,
			const general::time::JD& t0,
			const double step,
			general::math::PV& xk,
			general::time::JD& tk,
			const ArgsType& ... args) const
		{
			static_cast<const InvType*>(this)->integrate(x0, t0, step, xk, tk, args ...);
		}
	};

	template<class InvType, class ... ArgsType>
	class MultistepIntegrator : public Integrator<ArgsType ...>
	{
	private:
		const size_t _degree;

	public:
		explicit MultistepIntegrator(const size_t degree = 1) : Integrator<ArgsType ...>(), _degree(degree) {}
		~MultistepIntegrator() = default;

		void integrate(
			std::pair<general::math::PV, general::time::JD>* pData,
			const double step,
			general::math::PV& xk,
			general::time::JD& tk,
			const ArgsType& ... args) const
		{
			static_cast<const InvType*>(this)->integrate(pData, step, xk, tk, args ...);
		}

		size_t degree() const {	return _degree;	}
	};
}