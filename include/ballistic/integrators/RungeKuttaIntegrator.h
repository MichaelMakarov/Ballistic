#pragma once
#include "Integrators.h"

namespace ball
{
	template<Arithmetic R>
	class RKIntegrator : public SinglestepIntegrator<RKIntegrator<R>, R>
	{
	public:
		RKIntegrator() : SinglestepIntegrator<RKIntegrator<R>, R>() {}
		~RKIntegrator() = default;

		void integrate(
			const R& x0,
			const general::time::JD& t0,
			const double step,
			R& xk,
			general::time::JD& tk) const
		{
			const double 	step_2 = 0.5 * step,
							step_6 = step / 6;
			auto			t = tk = t0;
			t.add_seconds(static_cast<int>(step_2));
			tk.add_seconds(static_cast<int>(step));
			auto k1{ this->func(x0, t0) };
			auto k2{ this->func(x0 + step_2 * k1, t) };
			auto k3{ this->func(x0 + step_2 * k2, t) };
			auto k4{ this->func(x0 + step * k3, tk) };
			xk = x0 + step_6 * (k1 + 2 * (k2 + k3) + k4);
		}
	};
}