#pragma once
#include "SingleStepIntegrator.h"

namespace ball
{
	template<class ... ArgType>
	class RungeKuttaIntegrator : public SingleStepIntegrator<ArgType ...>
	{
	public:
		RungeKuttaIntegrator() {}
		~RungeKuttaIntegrator() {}

		types::PV Integrate(const double step, const ArgType ... args) const override
		{
			const double 	step_2 = 0.5 * step,
							step_6 = step / 6;
			time::JD		t1{ _t0 },
							t2{ _t0 };
			t1.AddSeconds(static_cast<int>(step_2));
			t2.AddSeconds(static_cast<int>(step));
			auto k1{ Func(_x0, _t0, args ...) };
			auto k2{ Func(_x0 + step_2 * k1, t1, args ...) };
			auto k3{ Func(_x0 + step_2 * k2, t1, args ...) };
			auto k4{ Func(_x0 + step * k3, t2, args ...) };
			return _x0 + step_6 * (k1 + 2 * (k2 + k3) + k4);
		}
	};
}