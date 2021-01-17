#pragma once
#include "Integrators.h"

namespace ball
{
	class RKIntegrator : public SinglestepIntegrator<RKIntegrator>
	{
	public:
		RKIntegrator() : SinglestepIntegrator() {}
		~RKIntegrator() {}

		void integrate(
			const geometry::PV& x0,
			const time::JD& t0,
			const double step,
			geometry::PV& xk,
			time::JD& tk) const
		{
			const double 	step_2 = 0.5 * step,
							step_6 = step / 6;
			time::JD		t{ t0 };
			t.add_seconds(static_cast<int>(step_2));
			tk = t0;
			tk.add_seconds(static_cast<int>(step));
			auto k1{ func(x0, t0) };
			auto k2{ func(x0 + step_2 * k1, t) };
			auto k3{ func(x0 + step_2 * k2, t) };
			auto k4{ func(x0 + step * k3, tk) };
			xk = x0 + step_6 * (k1 + 2 * (k2 + k3) + k4);
		}
	};
}