#pragma once
#include "SingleStepIntegrator.h"

namespace ball
{
	class RungeKuttaIntegrator : public SingleStepIntegrator
	{
	public:
//		using SingleStepIntegrator::Integrator;
		~RungeKuttaIntegrator() {}

		types::PV Integrate(const double step) const override;
	};
}