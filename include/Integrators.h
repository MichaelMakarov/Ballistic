#pragma once
#include "Geometry.h"
#include "Times.h"
#include <functional>

namespace ball
{
	class Integrator
	{
	public:
		Integrator() {}
		~Integrator() {}

		std::function<geometry::PV(const geometry::PV&, const time::JD& t)> func;
	};

	template<class InvType>
	class SinglestepIntegrator : public Integrator
	{
	public:
		SinglestepIntegrator() : Integrator() {}
		~SinglestepIntegrator() {}

		void integrate(
			const geometry::PV& x0,
			const time::JD& t0,
			const double step,
			geometry::PV& xk,
			time::JD& tk) const
		{
			static_cast<const InvType*>(this)->integrate(x0, t0, step, xk, tk);
		}
	};

	template<class InvType>
	class MultistepIntegrator : public Integrator
	{
	private:
		const size_t _degree;

	public:
		explicit MultistepIntegrator(const size_t degree = 1) : Integrator(), _degree(degree) {}
		~MultistepIntegrator() {}

		void integrate(
			std::pair<geometry::PV, time::JD>* pData,
			const double step,
			geometry::PV& xk,
			time::JD& tk) const
		{
			static_cast<const InvType*>(this)->integrate(pData, step, xk, tk);
		}

		size_t degree() const
		{
			return _degree;
		}
	};
}