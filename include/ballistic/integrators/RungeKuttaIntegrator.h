#pragma once
#include "Integrators.h"

namespace ball
{
	template<class ... ArgsType>
	class RKIntegrator : public SinglestepIntegrator<RKIntegrator<ArgsType ...>, ArgsType ...>
	{
	public:
		RKIntegrator() : SinglestepIntegrator<RKIntegrator<ArgsType ...>, ArgsType ...>() {}
		~RKIntegrator() = default;

		void integrate(
			const general::math::PV& x0,
			const general::time::JD& t0,
			const double step,
			general::math::PV& xk,
			general::time::JD& tk,
			const ArgsType& ... args) const
		{
			const double 	step_2 = 0.5 * step,
							step_6 = step / 6;
			general::time::JD	t{ t0 };
			t.add_seconds(static_cast<int>(step_2));
			tk = t0;
			tk.add_seconds(static_cast<int>(step));
			auto k1{ func(x0, t0, args ...) };
			auto k2{ func(x0 + step_2 * k1, t, args ...) };
			auto k3{ func(x0 + step_2 * k2, t, args ...) };
			auto k4{ func(x0 + step * k3, tk, args ...) };
			xk = x0 + step_6 * (k1 + 2 * (k2 + k3) + k4);
		}
	};
}