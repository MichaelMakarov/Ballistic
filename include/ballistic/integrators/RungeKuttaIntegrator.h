#pragma once
#include "Integrators.h"

namespace ball
{
	template<Arithmetic R, Time T>
	class rk_integrator : public singlestep_integrator<rk_integrator<R, T>, R, T>
	{
	public:
		template<class Inv> void integrate(
			const R& x0,
			const T& t0,
			const double step,
			R& xk,
			T& tk,
			const invoker<Inv, R, const R&, const T&>& func) const
		{
			const double step_2 = 0.5 * step, step_6 = step / 6;
			auto t = t0 + step_2;
			tk = t0 + step;
			auto k1{ func(x0, t0) };
			auto k2{ func(x0 + k1 * step_2, t) };
			auto k3{ func(x0 + k2 * step_2, t) };
			auto k4{ func(x0 + k3 * step,  tk) };
			xk = x0 + (k1 + (k2 + k3) * 2 + k4) * step_6;
		}
	};
}